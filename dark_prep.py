#! /usr/bin/env python

'''
Module for preparing a given dark current exposure for
integration with a ramp of simulated sources created
using (e.g.) catalog_seed_image.py. This is part of the 
refactored ramp_simulator.py

For the moment, just keep the same yaml input file
as is used for catalog_seed_image.py and ramp_simulator.py.
I'm still not sure if it is worth breaking out into 3 separate
input files for the refactored simulator...
'''

import argparse
import sys, os
import yaml

class DarkPrep:

    def __init__(self):
        pass

    def run(self):
        # Read in the yaml parameter file
        self.readParameterFile()

        # Expand locations to be full path names
        self.fullPaths()

        # Base name for output files
        self.basename = os.path.join(self.params['Output']['directory'],
                                     self.params['Output']['file'][0:-5])
        
        # Check the entered read pattern info
        self.readPatternCheck()

        # Check input parameters for any bad values
        self.checkParams()

        # Read in the subarray definition file
        self.readSubarrayDefinitionFile()

        # Define the subarray bounds from param file
        self.getSubarrayBounds()
        
        # Read in the input dark current frame
        if self.runStep['linearized_darkfile'] == False:
            self.getBaseDark()
            self.linDark = None
        else:
            print(self.params['Reffiles']['linearized_darkfile'])
            self.readLinearDark()
            self.dark = self.linDark

        # Make sure there is enough data (frames/groups)
        # in the input integration to produce
        # the proposed output integration
        self.dataVolumeCheck(self.dark)

        # Compare the requested number of integrations
        # to the number of integrations in the input dark
        print("Dark shape as read in: {}".format(self.dark.data.shape))
        self.darkints()
        print("Dark shape after copying integrations to match request: {}".format(self.dark.data.shape))

        # Put the input dark (or linearized dark) into the
        # requested readout pattern
        self.dark,sbzeroframe = self.reorderDark(self.dark)
        print('DARK has been reordered to {} to match the input readpattern of {}'.format(self.dark.data.shape,self.dark.header['READPATT']))

        # If a raw dark was read in, create linearized version
        # here using the SSB pipeline. Better to do this 
        # on the full dark before cropping, so that reference
        # pixels can be used in the processing.
        if ((self.params['Inst']['use_JWST_pipeline']) & (self.runStep['linearized_darkfile'] == False)):

            #Linear ize the dark ramp via the SSB pipeline.
            # Also save a diff image of the original dark minus
            # the superbias and refpix subtracted dark, to use later.

            # In order to linearize the dark, the JWST pipeline must
            # be present, and self.dark will have to be translated back
            # into a RampModel instance
            self.linDark = self.linearizeDark(self.dark)
            print("Linearized dark shape: {}".format(self.linDark.data.shape))

            if self.params['Readout']['readpatt'].upper() == 'RAPID':
                print("Output is rapid, grabbing zero frame from linearized dark")
                self.zeroModel = read_fits.Read_fits()
                #can't grab all zeroframes (for multiple ints)'
                self.zeroModel.data = self.linDark.data[:,0,:,:] 
                self.zeroModel.sbAndRefpix = self.linDark.sbAndRefpix[:,0,:,:]
            elif ((self.params['Readout']['readpatt'].upper() != 'RAPID') & (self.dark.zeroframe is not None)):
                print("Now we need to linearize the zeroframe because the")
                print("output readpattern is not RAPID")
                # Now we need to linearize the zeroframe. Place it
                # into a RampModel instance before running the
                # pipeline steps
                self.zeroModel = read_fits.Read_fits()
                self.zeroModel.data = self.dark.zeroframe
                self.zeroModel.header = self.linDark.header
                self.zeroModel.header['NGROUPS'] = 1
                self.zeroModel = self.linearizeDark(self.zeroModel)
            else:
                self.zeroModel = None

            # Now crop self.linDark, self.dark, and zeroModel
            # to requested subarray
            self.dark = self.cropDark(self.dark)
            self.linDark= self.cropDark(self.linDark)
            if self.zeroModel is not None:
                self.zeroModel = self.cropDark(self.zeroModel)
        else:

            # If no pipeline is run
            self.zeroModel = read_fits.Read_fits()
            self.zeroModel.data = self.dark.zeroframe
            self.zeroModel.sbAndRefpix = sbzeroframe

            # Crop the linearized dark to the requested
            # subarray size
            # THIS WILL CROP self.dark AS WELL SINCE
            # self.linDark IS JUST A REFERENCE IN THE NON
            # PIPELINE CASE!!
            self.linDark = self.cropDark(self.linDark)
            if self.zeroModel.data is not None:
                self.zeroModel = self.cropDark(self.zeroModel)

        #save the linearized dark for testing
        if self.params['Output']['save_intermediates']:
            h0=fits.PrimaryHDU()
            h1 = fits.ImageHDU(self.linDark.data)
            h2 = fits.ImageHDU(self.linDark.sbAndRefPix)
            h3 = fits.ImageHDU(self.zeroModel.data)
            h4 = fits.ImageHDU(self.zeroModel.sbAndRefPix)
            hl=fits.HDUList([h0,h1,h2,h3,h4])
            hl.writeto(self.basename + '_linearizedDark.fits',overwrite=True)

        #important variables
        #self.linDark
        #self.zeroModel


    def fullPaths(self):
        # Expand all input paths to be full paths
        # This is to allow for easier Condor-ization of
        # many runs
        pathdict = {'Reffiles':['dark','linearized_darkfile','hotpixmask','superbias',
                                'subarray_defs','readpattdefs','linearity',
                                'saturation','gain','phot','pixelflat',
                                'illumflat','astrometric','distortion_coeffs','ipc',
                                'crosstalk','occult','filtpupilcombo','pixelAreaMap',
                                'flux_cal'],
                    'cosmicRay':['path'],
                    'simSignals':['pointsource','psfpath','galaxyListFile','extended',
                                  'movingTargetList','movingTargetSersic',
                                  'movingTargetExtended','movingTargetToTrack'],
                    'newRamp':['dq_configfile','sat_configfile','superbias_configfile',
                               'refpix_configfile','linear_configfile'],
                    'Output':['file','directory']}

        for key1 in pathdict:
            for key2 in pathdict[key1]:
                if self.params[key1][key2].lower() != 'none':
                    self.params[key1][key2] = os.path.abspath(self.params[key1][key2])

        
    def cropDark(self,model):
        # Cut the dark current array down to the size dictated
        # by the subarray bounds
        modshape = model.data.shape
        yd = modshape[-2]
        xd = modshape[-1]
        
        if ((self.subarray_bounds[0] != 0) or (self.subarray_bounds[2] != (xd-1))
            or (self.subarray_bounds[1] != 0) or (self.subarray_bounds[3] != (yd-1))):

            if len(modshape) == 4:
                model.data = model.data[:,:,self.subarray_bounds[1]:self.subarray_bounds[3]+1,
                                        self.subarray_bounds[0]:self.subarray_bounds[2]+1]
                try:
                    model.sbAndRefpix = model.sbAndRefpix[:,:,self.subarray_bounds[1]:
                                                          self.subarray_bounds[3]+1,
                                                          self.subarray_bounds[0]:self.subarray_bounds[2]+1]
                except:
                    pass
                try:
                    model.zeroframe = model.zeroframe[:,self.subarray_bounds[1]:
                                                          self.subarray_bounds[3]+1,
                                                          self.subarray_bounds[0]:self.subarray_bounds[2]+1]
                except:
                    pass
            if len(modshape) == 3:
                model.data = model.data[:,self.subarray_bounds[1]:self.subarray_bounds[3]+1,
                                        self.subarray_bounds[0]:self.subarray_bounds[2]+1]
                try:
                    model.sbAndRefpix = model.sbAndRefpix[:,self.subarray_bounds[1]:
                                                          self.subarray_bounds[3]+1,
                                                          self.subarray_bounds[0]:self.subarray_bounds[2]+1]
                except:
                    pass
                try:
                    model.zeroframe = model.zeroframe[:,self.subarray_bounds[1]:
                                                      self.subarray_bounds[3]+1,
                                                      self.subarray_bounds[0]:self.subarray_bounds[2]+1]
                except:
                    pass
            elif len(modshape) == 2:
                model.data = model.data[self.subarray_bounds[1]:self.subarray_bounds[3]+1,
                                        self.subarray_bounds[0]:self.subarray_bounds[2]+1]
                try:
                    model.sbAndRefpix = model.sbAndRefpix[self.subarray_bounds[1]:self.subarray_bounds[3]+1,
                                                          self.subarray_bounds[0]:self.subarray_bounds[2]+1]
                except:
                    pass
                try:
                    model.zeroframe = model.zeroframe[self.subarray_bounds[1]:
                                                      self.subarray_bounds[3]+1,
                                                      self.subarray_bounds[0]:self.subarray_bounds[2]+1]
                except:
                    pass

        # Make sure that if the output is supposedly a
        # 4-amplifier file, that the number of
        # pixels in the x direction is a multiple of 4.
        nfast=self.subarray_bounds[2]-self.subarray_bounds[0]+1
        nslow=self.subarray_bounds[3]-self.subarray_bounds[1]+1
        nramp=(self.params['Readout']['nframe']+self.params['Readout']['nskip'])*self.params['Readout']['ngroup']

        if self.params['Readout']['namp'] != 4 and self.params['Readout']['namp'] != 1:
            print('ERROR: amplifier mode specified ({}) is not allowed'.format(self.params['Readout']['namp']))
            sys.exit()

        if self.params['Readout']['namp'] == 4:
            n=int(nfast/4)*4
            if n != nfast:
                print('ERROR: 4 amplifier mode specified but the number of pixels in the fast\nread direction ({}) is not a multiple of 4.'.format(nfast))
                sys.exit()

        return model

             
    def reorderDark(self,dark):
        # Reorder the input dark ramp using the requested
        # readout pattern (nframe,nskip). If the initial
        # dark ramp is RAPID, then save and return the 0th frame.

        if self.params['Reffiles']['linearized_darkfile']:
            datatype = np.float
        else:
            datatype = np.int32
            
        # Get the info for the dark integration
        darkpatt = dark.header['READPATT']
        dark_nframe = dark.header['NFRAMES']
        mtch = self.readpatterns['name'].data == darkpatt
        dark_nskip = self.readpatterns['nskip'].data[mtch][0]

        nint,ngroup,yd,xd = dark.data.shape
        outdark = np.zeros((self.params['Readout']['nint'],self.params['Readout']['ngroup'],yd,xd))

        if dark.sbAndRefpix is not None:
            outsb = np.zeros((self.params['Readout']['nint'],self.params['Readout']['ngroup'],yd,xd))
            
        # We can only keep a zero frame around if the input dark
        # is RAPID. Otherwise that information is lost.
        darkzero = None
        sbzero = None
        if ((darkpatt == 'RAPID') & (dark.zeroframe is None)):
            dark.zeroframe = dark.data[:,0,:,:]
            print("Saving 0th frame from data to the zeroframe extension")
            if dark.sbAndRefpix is not None:
                sbzero = dark.sbAndRefpix[:,0,:,:]
        elif ((darkpatt != 'RAPID') & (dark.zeroframe is None)):
            print("Unable to save the zeroth frame because the input dark current ramp is not RAPID.")
            sbzero = None            
        elif ((darkpatt == 'RAPID') & (dark.zeroframe is not None)):
            # In this case we already have the zeroframe
            if dark.sbAndRefpix is not None:
                sbzero = dark.sbAndRefpix[:,0,:,:]
        elif ((darkpatt != 'RAPID') & (dark.zeroframe is not None)):
            # Non RAPID dark, zeroframe is present, but
            # we can't get sbAndRefpix for the zeroth frame
            # because the pattern is not RAPID.
            sbzero = None
                
        # We have already guaranteed that either the readpatterns match
        # or the dark is RAPID, so no need to worry about checking for 
        # other cases here.

        if ((darkpatt == 'RAPID') and (self.params['Readout']['readpatt'] != 'RAPID')): 

            deltaframe = self.params['Readout']['nskip']+self.params['Readout']['nframe']
            frames = np.arange(0,self.params['Readout']['nframe'])
            accumimage = np.zeros_like(outdark[0,0,:,:],dtype=datatype)

            if dark.sbAndRefpix is not None:
                zeroaccumimage = np.zeros_like(outdark[0,0,:,:],dtype=np.float)
            
            # Loop over integrations
            for integ in range(self.params['Readout']['nint']):
                
                # Loop over groups
                for i in range(self.params['Readout']['ngroup']):
                    # Average together the appropriate frames,
                    # skip the appropriate frames
                    print('Averaging dark current ramp. Frames {}, to become group {}'.format(frames,i))

                    # If averaging needs to be done
                    if self.params['Readout']['nframe'] > 1:
                        accumimage = np.mean(dark.data[integ,frames,:,:],axis=0)
                        if dark.sbAndRefpix is not None:
                            zeroaccumimage = np.mean(dark.sbAndRefpix[integ,frames,:,:],axis=0)
                        
                        # If no averaging needs to be done
                    else:
                        accumimage = dark.data[integ,frames[0],:,:]
                        if dark.sbAndRefpix is not None:
                            zeroaccumimage = dark.sbAndRefpix[integ,frames[0],:,:]
                            
                    outdark[integ,i,:,:] += accumimage
                    if dark.sbAndRefpix is not None:
                        outsb[integ,i,:,:] += zeroaccumimage
                    
                    # Increment the frame indexes
                    frames = frames + deltaframe

        elif (self.params['Readout']['readpatt'] == darkpatt):
            # If the input dark is not RAPID, or if the readout
            # pattern of the input dark and
            # the output ramp match, then no averaging needs to
            # be done 
            outdark = dark.data[:,0:self.params['Readout']['ngroup'],:,:]
            if dark.sbAndRefpix is not None:
                outsb = dark.sbAndRefpix[:,0:self.params['Readout']['ngroup'],:,:]
        else:
            # This check should already have been done,
            # but just to be sure...
            print("WARNING: dark current readout pattern is {} and requested output is {}.".format(darkpatt,self.params['Readout']['readpatt']))
            print("Cannot convert between the two.")
            sys.exit()

        # Now place the reorganized dark into the data object and
        # update the appropriate metadata
        dark.data = outdark
        if dark.sbAndRefpix is not None:
            dark.sbAndRefpix = outsb
        dark.header['READPATT'] = self.params['Readout']['readpatt'] 
        dark.header['NFRAMES'] = self.params['Readout']['nframe']
        dark.header['NSKIP'] = self.params['Readout']['nskip']
        dark.header['NGROUPS'] = self.params['Readout']['ngroup']
        dark.header['NINTS'] = self.params['Readout']['nint']

        return dark,sbzero


        
    def darkints(self):
        # Check the number of integrations in the dark
        # current file and compare with the requested
        # number of integrations. Add/remove integrations
        # if necessary
        ndarkints,ngroups,ny,nx = self.dark.data.shape
        reqints = self.params['Readout']['nint']

        if reqints <= ndarkints:
            # Fewer requested integrations than are in the
            # input dark.
            print("{} output integrations requested.".format(reqints))
            self.dark.data = self.dark.data[0:reqints,:,:,:]
            if self.dark.sbAndRefpix is not None:
                self.dark.sbAndRefpix = self.dark.sbAndRefpix[0:reqints,:,:,:]
            if self.dark.zeroframe is not None:
                self.dark.zeroframe = self.dark.zeroframe[0:reqints,:,:]

        elif reqints > ndarkints:
            # More requested integrations than are in input dark.
            print('Requested output has {} integrations, while input dark has only {}.'.format(reqints,ndarkints))
            print('Adding integrations to input dark by making copies of the input.')
            self.integration_copy(reqints,ndarkints)


    def dataVolumeCheck(self,obj):
        # Make sure that the input integration has
        # enough frames/groups to create the requested
        # number of frames/groups of the output
        ngroup = int(self.params['Readout']['ngroup'])
        nframe = int(self.params['Readout']['nframe'])
        nskip = int(self.params['Readout']['nskip'])

        inputframes = obj.data.shape[1]
        if ngroup*(nskip+nframe) > inputframes:
            print("WARNING: Not enough frames in the input integration to create the requested number of output groups. Input has {} frames. Requested output is {} groups each created from {} frames plus skipping {} frames between groups.".format(inputframes,ngroup,nframe,nskip))
            print("for a total of {} frames.".format(ngroup*(nframe+nskip)))
            print("Making copies of {} dark current frames and adding them to the end of the dark current integration.".format(ngroup*(nskip+nframe) - inputframes))

            # Figure out how many more frames we need,
            # in terms of how many copies of the original dark
            div = (ngroup*(nskip+nframe)) / inputframes
            mod = (ngroup*(nskip+nframe)) % inputframes
            
            # If more frames are needed than there are frames
            # in the original dark, then make copies of the
            # entire thing as many times as necessary, adding
            # the signal from the previous final frame to each.
            for ints in xrange(div-1):
                extra_frames = np.copy(obj.data)
                obj.data = np.hstack((obj.data,extra_frames+obj.data[:,-1,:,:]))
                if obj.sbAndRefpix is not None:
                    extra_sb = np.copy(obj.sbAndRefpix)
                    obj.sbAndRefpix = np.hstack((obj.sbAndRefpix,extra_sb+obj.sbAndRefpix[:,-1,:,:]))
                
            # At this point, if more frames are needed, but fewer
            # than an entire copy of self.dark.data, then add the
            # appropriate number of frames here.
            extra_frames = np.copy(obj.data[:,1:mod+1,:,:]) - obj.data[:,0,:,:]
            obj.data = np.hstack((obj.data,extra_frames+obj.data[:,-1,:,:]))
            if obj.sbAndRefpix is not None:
                extra_sb = np.copy(obj.sbAndRefpix[:,1:mod+1,:,:]) - obj.sbAndRefpix[:,0,:,:]
                obj.sbAndRefpix = np.hstack((obj.sbAndRefpix,extra_sb+obj.sbAndRefpix[:,-1,:,:]))
            
        elif ngroup*(nskip+nframe) < inputframes:
            # If there are more frames in the dark than we'll need,
            # crop the extras in order to reduce memory use
            obj.data = obj.data[:,0:ngroup*(nskip+nframe),:,:]
            if obj.sbAndRefpix is not None:
                obj.sbAndRefpix = obj.sbAndRefpix[:,0:ngroup*(nskip+nframe),:,:]
        obj.header['NGROUPS'] = ngroup*(nskip+nframe)

        
    def getSubarrayBounds(self):
        #find the bounds of the requested subarray
        if self.params['Readout']['array_name'] in self.subdict['AperName']:
            mtch = self.params['Readout']['array_name'] == self.subdict['AperName']
            self.subarray_bounds = [self.subdict['xstart'].data[mtch][0],self.subdict['ystart'].data[mtch][0],self.subdict['xend'].data[mtch][0],self.subdict['yend'].data[mtch][0]]
            self.refpix_pos = {'x':self.subdict['refpix_x'].data[mtch][0],'y':self.subdict['refpix_y'][mtch][0],'v2':self.subdict['refpix_v2'].data[mtch][0],'v3':self.subdict['refpix_v3'].data[mtch][0]}

            namps = self.subdict['num_amps'].data[mtch][0]
            if namps != 0:
                self.params['Readout']['namp'] = namps
            else:
                if ((self.params['Readout']['namp'] == 1) or (self.params['Readout']['namp'] == 4)):
                    print("CAUTION: Aperture {} can be used with either a 1-amp".format(self.subdict['AperName'].data[mtch][0]))
                    print("or a 4-amp readout. The difference is a factor of 4 in")
                    print("readout time. You have requested {} amps.".format(self.params['Readout']['namp']))
                else:
                    print("WARNING: {} requires the number of amps to be 1 or 4. You have requested {}.".format(self.params['Readout']['array_name'],self.params['Readout']['namp']))
                    sys.exit()


        else:
            print("WARNING: subarray name {} not found in the subarray dictionary {}.".format(self.params['Readout']['array_name'],self.params['Reffiles']['subarray_defs']))
            sys.exit()
        

    def readSubarrayDefinitionFile(self):
        #read in the file that contains a list of subarray names and positions on the detector
        try:
            self.subdict = ascii.read(self.params['Reffiles']['subarray_defs'],data_start=1,header_start=0)
        except:
            print("Error: could not read in subarray definitions file.")
            sys.exit()


    def readParameterFile(self):
        #read in the parameter file
        try:
            with open(self.paramfile,'r') as infile:
                self.params = yaml.load(infile)
        except:
            print("WARNING: unable to open {}".format(self.paramfile))
            sys.exit()


        
    def readpatternCheck(self):
        '''check the readout pattern that's entered and set nframe and nskip 
           accordingly'''
        self.params['Readout']['readpatt'] = self.params['Readout']['readpatt'].upper()

        #read in readout pattern definition file 
        #and make sure the possible readout patterns are in upper case
        self.readpatterns = ascii.read(self.params['Reffiles']['readpattdefs'])
        self.readpatterns['name'] = [s.upper() for s in self.readpatterns['name']]

        #if the requested readout pattern is in the table of options,
        #then adopt the appropriate nframe and nskip
        if self.params['Readout']['readpatt'] in self.readpatterns['name']:
            mtch = self.params['Readout']['readpatt'] == self.readpatterns['name']
            self.params['Readout']['nframe'] = self.readpatterns['nframe'][mtch].data[0]
            self.params['Readout']['nskip'] = self.readpatterns['nskip'][mtch].data[0]
            print('Requested readout pattern {} found in {}. Using the nframe and nskip'.format(self.params['Readout']['readpatt'],self.params['Reffiles']['readpattdefs']))
            print('values from the definition file. nframe = {} and nskip = {}.'.format(self.params['Readout']['nframe'],self.params['Readout']['nskip']))
        else:
            #if readpatt is not present in the definition file but the nframe/nskip combo is, then reset 
            #readpatt to the appropriate value from the definition file
            readpatt_nframe = self.readpatterns['nframe'].data
            readpatt_nskip = self.readpatterns['nskip'].data
            readpatt_name = self.readpatterns['name'].data
            if self.params['Readout']['nframe'] in readpatt_nframe:
                nfmtch = self.params['Readout']['nframe'] == readpatt_nframe
                nskip_subset = readpatt_nskip[nfmtch]
                name_subset = readpatt_name[nfmtch]
                if self.params['Readout']['nskip'] in nskip_subset:
                    finalmtch = self.params['Readout']['nskip'] == nskip_subset
                    finalname = name_subset[finalmtch][0]
                    print("CAUTION: requested readout pattern {} not recognized.".format(self.params['Readout']['readpatt']))
                    print('but the requested nframe/nskip combination ({},{}), matches those values for'.format(self.params['Readout']['nframe'],self.params['Readout']['nskip']))
                    print("the {} readout pattern, as listed in {}.".format(finalname,self.params['Reffiles']['readpattdefs']))
                    print('Continuing on using that as the readout pattern.')
                    self.params['Readout']['readpatt'] = finalname
                else:
                    #case where readpatt is not recognized, nframe is present in the definition file, but the nframe/nskip combination is not
                    print('Unrecognized readout pattern {}, and the input nframe/nskip combination {},{} does not'.format(self.params['Readout']['readpatt'],self.params['Readout']['nframe'],self.params['Readout']['nskip']))
                    print('match any present in {}. Continuing simulation with the input nframe/nskip, and'.format(self.params['Reffiles']['readpattdefs']))
                    print("setting the readout pattern to 'ANY' in order to allow the output file to be saved via RampModel without error.")
                    self.params['readpatt'] = 'ANY'
            else:
                #case where readpatt is not recognized, and nframe is not present in the definition file
                print('Unrecognized readout pattern {}, and the input nframe/nskip combination {},{} does not'.format(self.params['Readout']['readpatt'],self.params['Readout']['nframe'],self.params['Readout']['nskip']))
                print('match any present in {}. '.format(self.params['Reffiles']['readpattdefs']))
                print('Continuing simulation with the input nframe/nskip values, and')
                print("setting the readout pattern to 'ANY' in order to allow the output file to be saved via RampModel without error.")
                self.params['readpatt'] = 'ANY'

        
    def checkParams(self):
        # Check instrument name
        if self.params['Inst']['instrument'].lower() not in inst_list:
            print("WARNING: {} instrument not implemented within ramp simulator")
            sys.exit()

        # If user requests not to use the pipeline,
        # make sure the input is a linearized dark, and not
        # a raw dark
        if self.params['Inst']['use_JWST_pipeline'] == False:
            if self.params['Reffiles']['linearized_darkfile'] is None:
                print("WARNING: You specified no use of the JWST pipeline, but have")
                print("not provided a linearized dark file to use. Without the")
                print("pipeline, a raw dark cannot be used.")
                sys.exit()
            
        # Make sure nframe,nskip,ngroup are all integers
        try:
            self.params['Readout']['nframe'] = int(self.params['Readout']['nframe'])
        except:
            print("WARNING: Input value of nframe is not an integer.")
            sys.exit

        try:
            self.params['Readout']['nskip'] = int(self.params['Readout']['nskip'])
        except:
            print("WARNING: Input value of nskip is not an integer.")
            sys.exit

        try:
            self.params['Readout']['ngroup'] = int(self.params['Readout']['ngroup'])
        except:
            print("WARNING: Input value of ngroup is not an integer.")
            sys.exit

        try:
            self.params['Readout']['nint'] = int(self.params['Readout']['nint'])
        except:
            print("WARNING: Input value of nint is not an integer.")
            sys.exit


        # Make sure that the requested number of groups is
        # less than or equal to the maximum allowed. If you're
        # continuing on with an unknown readout pattern (not
        # recommended) then assume a max of 10 groups.
        # For science operations, ngroup is going to be limited
        # to 10 for all readout patterns except for the DEEP
        # patterns, which can go to 20. (FOR FULL FRAME!!)
        # Subarray exposures can use more!!
        match = self.readpatterns['name'] == self.params['Readout']['readpatt'].upper()
        if sum(match) == 1:
            maxgroups = self.readpatterns['maxgroups'].data[match][0]
        if sum(match) == 0:
            print("Unrecognized readout pattern {}. Assuming a maximum allowed number of groups of 10.".format(self.params['Readout']['readpatt']))
            maxgroups = 10

        if (self.params['Readout']['ngroup'] > maxgroups):
            print("WARNING: {} is limited to a maximum of {} groups. Proceeding with ngroup = {}.".format(self.params['Readout']['readpatt'],maxgroups,maxgroups))
            self.params['Readout']['readpatt'] = maxgroups

        # Check for entries in the parameter file that are None or blank,
        # indicating the step should be skipped. Create a dictionary of steps 
        # and populate with True or False
        self.runStep = {}
        #self.runStep['linearity'] = self.checkRunStep(self.params['Reffiles']['linearity'])
        self.runStep['linearized_darkfile'] = self.checkRunStep(self.params['Reffiles']['linearized_darkfile'])

        # NON-LINEARITY
        #make sure the input accuracy is a float with reasonable bounds
        #self.params['nonlin']['accuracy'] = self.checkParamVal(self.params['nonlin']['accuracy'],'nlin accuracy',1e-12,1e-6,1e-6)
        #self.params['nonlin']['maxiter'] = self.checkParamVal(self.params['nonlin']['maxiter'],'nonlin max iterations',5,40,10)
        #self.params['nonlin']['limit'] = self.checkParamVal(self.params['nonlin']['limit'],'nonlin max value',30000.,1.e6,66000.)
    
        # If the pipeline is going to be used to create
        # the linearized dark current ramp, make sure the
        # specified configuration files are present.
        if self.runStep['linearized_darkfile'] == False:
            dqcheck = self.checkRunStep(self.params['newRamp']['dq_configfile'])
            if dqcheck == False:
                print("WARNING: DQ pipeline step configuration file not provided. This file is needed to run the pipeline.")
                sys.exit()
            satcheck = self.checkRunStep(self.params['newRamp']['sat_configfile'])
            if satcheck == False:
                print("WARNING: Saturation pipeline step configuration file not provided. This file is needed to run the pipeline.")
                sys.exit()
            sbcheck = self.checkRunStep(self.params['newRamp']['superbias_configfile'])
            if sbcheck == False:
                print("WARNING: Superbias pipeline step configuration file not provided. This file is needed to run the pipeline.")
                sys.exit()
            refpixcheck = self.checkRunStep(self.params['newRamp']['refpix_configfile'])
            if refpixcheck == False:
                print("WARNING: Refpix pipeline step configuration file not provided. This file is needed to run the pipeline.")
                sys.exit()
            lincheck = self.checkRunStep(self.params['newRamp']['linear_configfile'])
            if lincheck == False:
                print("WARNING: Linearity pipeline step configuration file not provided. This file is needed to run the pipeline.")
                sys.exit()


