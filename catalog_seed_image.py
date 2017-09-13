#! /usr/bin env python

'''
Reorganization of ramp simulator code. This class is used to construct
a "seed image" from a point source catalog. This seed image is a 
noiseless countrate image containing only sources. No noise, no 
cosmic rays.
'''


import argparse, sys, glob, os
import yaml
import scipy.signal as s1
import numpy as np
import math
from astropy.io import fits, ascii
from astropy.table import Table
from astropy.modeling.models import Sersic2D
import rotations  #this is Colin's set of rotation functions
from asdf import AsdfFile
import polynomial #more of Colin's functions
from astropy.convolution import convolve
import read_siaf_table
import set_telescope_pointing_separated as set_telescope_pointing
import moving_targets
import segmentation_map as segmap

inst_list = ['nircam']
modes = {'nircam':['imaging','moving_target','wfss','coron']}
inst_abbrev = {'nircam':'NRC'}
pixelScale = {'nircam':{'sw':0.031,'lw':0.063}}
full_array_size = {'nircam':2048}
allowedOutputFormats = ['DMS']



class Catalog_seed():#RampSim):
    def __init__(self):
        #RampSim.__init__(self) #creates grism_direct_factor, coord_adjust

        #if a grism signal rate image is requested, expand
        #the width and height of the signal rate image by this 
        #factor, so that the grism simulation software can 
        #track sources that are outside the requested subarray
        #in order to calculate contamination.
        self.grism_direct_factor = np.sqrt(2.)

        #self.coord_adjust contains the factor by which the nominal output array size
        #needs to be increased (used for TSO and WFSS modes), as well as the coordinate
        #offset between the nominal output array coordinates, and those of the expanded 
        #array. These are needed mostly for TSO observations, where the nominal output  
        #array will not sit centered in the expanded output image.
        self.coord_adjust = {'x':1.,'xoffset':0.,'y':1.,'yoffset':0.}

        # NIRCam rough noise values. Used to make educated guesses when
        # creating segmentation maps
        self.single_ron = 6. #e-/read
        self.grism_background = 0.25 #e-/sec
        

    def run(self):
        # Read in input parameters and quality check
        self.readParameterFile()
        self.fullPaths()
        self.basename = os.path.join(self.params['Output']['directory'],
                                     self.params['Output']['file'][0:-5])
        self.checkParams()
        self.readSubarrayDefinitionFile()
        self.getSubarrayBounds()
        self.instrument_specific_dicts(self.params['Inst']['instrument'].lower())
        
        #If the output is a TSO ramp, use the slew rate and angle (and whether a grism
        #direct image is requested) to determine how much larger than nominal the
        #signal rate image should be
        if self.params['Inst']['mode'] == 'moving_target' or self.params['Output']['grism_source_image']:
            self.calcCoordAdjust()

        #image dimensions
        self.nominal_dims = np.array([self.subarray_bounds[3]-self.subarray_bounds[1]+1,self.subarray_bounds[2]-self.subarray_bounds[0]+1])
        self.output_dims = (self.nominal_dims * np.array([self.coord_adjust['y'],self.coord_adjust['x']])).astype(np.int)

        #calculate the exposure time of a single frame, based on the size of the subarray
        self.calcFrameTime()

        # For imaging mode, generate the countrate image using the catalogs
        if self.params['Inst']['mode'].lower() != 'moving_target':
            print('Creating signal rate image of synthetic inputs.')
            self.seedimage, self.seed_segmap = self.addedSignals()
            outapp = ''
            
        # If we are tracking a non-sidereal target, then
        # everything in the catalogs needs to be streaked across
        # the detector
        if self.params['Inst']['mode'].lower() == 'moving_target':
            print('Creating signal ramp of synthetic inputs')
            print("need to adjust moving target work for multiple integrations! everything above has been modified")
            self.seedimage = self.non_sidereal_seed()
            outapp = '_nonsidereal_target'

        #if moving targets are requested (KBOs, asteroids, etc, NOT moving_target mode
        #where the telescope slews), then create a RAPID integration which 
        #includes those targets
        mov_targs_ramps = []
        if (self.runStep['movingTargets'] | self.runStep['movingTargetsSersic'] | self.runStep['movingTargetsExtended']):
            print('Creating signal ramp of sources that are moving with respect to telescope tracking.')
            trailed_ramp = self.make_trailed_ramp()
            outapp += '_trailed_sources'

            # Now we need to expand frameimage into a ramp
            # so we can add the trailed objects
            print('Combining trailed object ramp with that containing tracked targets')
            if self.params['Inst']['mode'].lower() != 'moving_target':
                self.seedimage = self.combineSimulatedDataSources('countrate',self.seedimage,trailed_ramp)
            else:
                self.seedimage = self.combineSimulatedDataSources('ramp',self.seedimage,trailed_ramp)
                
        #save the combined static+moving targets ramp
        self.saveSeedImage()


    def saveSeedImage(self):
        #Create the grism direct image or ramp to be saved 

        #Get the photflambda and photfnu values that go with
        #the filter
        module = self.params['Readout']['array_name'][3]
            
        zps = ascii.read(self.params['Reffiles']['flux_cal'])
        if self.params['Readout']['pupil'][0] == 'F':
            usephot = 'pupil'
        else:
            usephot = 'filter'
        mtch = ((zps['Filter'] == self.params['Readout'][usephot]) & (zps['Module'] == module))
        photflam = zps['PHOTFLAM'][mtch][0]
        photfnu = zps['PHOTFNU'][mtch][0]
        pivot = zps['Pivot_wave'][mtch][0]
        
        arrayshape = self.seedimage.shape
        if len(arrayshape) == 2:
            units = 'e-/sec'
            yd,xd = arrayshape
            tgroup = 0.
        elif len(arrayshape) == 3:
            units = 'e-'
            g,yd,xd = arrayshape
            tgroup = self.frametime*(self.params['Readout']['nframe']+self.params['Readout']['nskip'])
            
        seedName = os.path.join(self.basename + '_' + self.params['Readout']['filter'] + '_seed_image.fits')
        xcent_fov = xd / 2
        ycent_fov = yd / 2
        kw = {}
        kw['xcenter'] = xcent_fov
        kw['ycenter'] = ycent_fov
        kw['units'] = units
        kw['TGROUP'] = tgroup
        if self.params['Readout']['pupil'][0].upper() == 'F':
            usefilt = 'pupil'
        else:
            usefilt = 'filter'
        kw['filter'] = self.params['Readout'][usefilt]
        kw['PHOTFLAM'] = photflam
        kw['PHOTFNU'] = photfnu
        kw['PHOTPLAM'] = pivot * 1.e4 #put into angstroms
        self.seedinfo = kw
        self.saveSingleFits(self.seedimage,seedName,key_dict=kw,image2=self.seed_segmap,image2type='SEGMAP')
        print("Seed image and segmentation map saved as {}".format(seedName))
        print("Seed image, segmentation map, and metadata available as:")
        print("self.seedimage, self.seed_segmap, self.seedinfo.")

        
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

        
    def combineSimulatedDataSources(self,inputtype,input1,mov_tar_ramp):
        #inputtype can be 'countrate' in which case input needs to be made
        #into a ramp before combining with mov_tar_ramp, or 'ramp' in which
        #case you can combine directly. Use 'ramp' with
        #moving_target MODE data, and 'countrate' with imaging MODE data

        #Combine the countrate image with the moving target ramp.

        if inputtype == 'countrate':
            #First change the countrate image into a ramp
            yd,xd = input1.shape
            num_frames = self.params['Readout']['ngroup'] * (self.params['Readout']['nframe'] + self.params['Readout']['nskip'])
            print("Countrate image of synthetic signals being converted to RAPID integration with {} frames.".format(num_frames))
            input1_ramp = np.zeros((num_frames,yd,xd))
            for i in xrange(num_frames):
                input1_ramp[i,:,:] = input1 * self.frametime * (i+1)

        else:
            #if input1 is a ramp rather than a countrate image
            input1_ramp = input1

        #combine the input1 ramp and the moving target ramp, which are
        #now both RAPID mode
        totalinput = input1_ramp + mov_tar_ramp
        return totalinput

            

    def make_trailed_ramp(self):
        # Create a ramp for objects that are trailing through
        # the field of view during the integration
        mov_targs_ramps = []

        if self.params['Inst']['mode'].lower() != 'moving_target':
            tracking = False
            ra_vel = None
            dec_vel = None
        else:
            tracking = True
            ra_vel = self.ra_vel
            dec_vel = self.dec_vel
            #print("Moving target mode, creating trailed object with")
            #print("RA velocity of {} and dec_val of {}".format(ra_vel,dec_vel))
        
        if self.runStep['movingTargets']:
            #print('Starting moving targets for point sources!')
            mov_targs_ptsrc = self.movingTargetInputs(self.params['simSignals']['movingTargetList'],'pointSource',MT_tracking=tracking,tracking_ra_vel=ra_vel,tracking_dec_vel=dec_vel)
            mov_targs_ramps.append(mov_targs_ptsrc)

        #moving target using a sersic object
        if self.runStep['movingTargetsSersic']:
            #print("Moving targets, sersic!")
            mov_targs_sersic = self.movingTargetInputs(self.params['simSignals']['movingTargetSersic'],'galaxies',MT_tracking=tracking,tracking_ra_vel=ra_vel,tracking_dec_vel=dec_val)
            mov_targs_ramps.append(mov_targs_sersic)

        #moving target using an extended object
        if self.runStep['movingTargetsExtended']:
            #print("Extended moving targets!!!")
            mov_targs_ext = self.movingTargetInputs(self.params['simSignals']['movingTargetExtended'],'extended',MT_tracking=tracking,ra_vel=tracking_ra_vel,dec_vel=tracking_dec_val)
            mov_targs_ramps.append(mov_targs_ext)

        mov_targs_integration = None
        if self.runStep['movingTargets'] or self.runStep['movingTargetsSersic'] or self.runStep['movingTargetsExtended']:
            #Combine the ramps of the moving targets if there is more than one type
            mov_targs_integration = mov_targs_ramps[0]
            if len(mov_targs_ramps) > 1:
                for i in xrange(1,len(mov_targs_ramps)):
                    mov_targs_integration += mov_targs_ramps[0]
        return mov_targs_integration

            
    def calcFrameTime(self):
        #calculate the exposure time of a single frame of the proposed output ramp
        #based on the size of the croped dark current integration
        #numint,numgrp,yd,xd = self.dark.data.shape
        yd,xd = self.nominal_dims
        self.frametime = (xd/self.params['Readout']['namp'] + 12.) * (yd+1) * 10.00 * 1.e-6
        

    def calcCoordAdjust(self):
        # Calculate the factors by which to expand the output array size, as well as the coordinate
        # offsets between the nominal output array and the input lists if the observation being
        # modeled is TSO or TSO+grism output modes

        dtor = math.radians(1.)

        # Normal imaging with grism image requested
        if self.params['Output']['grism_source_image']:
            self.coord_adjust['x'] = self.grism_direct_factor
            self.coord_adjust['y'] = self.grism_direct_factor
            self.coord_adjust['xoffset'] = np.int((self.grism_direct_factor - 1.) * (self.subarray_bounds[2]-self.subarray_bounds[0]+1) / 2.)
            self.coord_adjust['yoffset'] = np.int((self.grism_direct_factor - 1.) * (self.subarray_bounds[3]-self.subarray_bounds[1]+1) / 2.)


            
    def non_sidereal_seed(self):
        # Create a seed RAMP in the case where NIRCam is tracking
        # a non-sidereal target
        
        # Create a count rate image containing only the non-sidereal target(s)
        # These will be stationary in the fov 
        nonsidereal_countrate,self.ra_vel,self.dec_vel,vel_flag = self.nonsidereal_CRImage(self.params['simSignals']['movingTargetToTrack'])

        # Expand into a RAPID ramp and convert from signal rate to signals
        ns_yd,ns_xd = nonsidereal_countrate.shape
        totframes = self.params['Readout']['ngroup'] * (self.params['Readout']['nframe']+self.params['Readout']['nskip'])
        tmptimes = self.frametime * np.arange(1,totframes+1)
        non_sidereal_ramp = np.zeros((totframes,ns_yd,ns_xd))
        for i in xrange(totframes):
            non_sidereal_ramp[i,:,:] = nonsidereal_countrate * tmptimes[i]
        #non_sidereal_zero = non_sidereal_ramp[0,:,:]

        #Now we need to collect all the other sources (point sources, galaxies, extended)
        #in the other input files, and treat them as targets which will move across
        #the field of view during the exposure.
        mtt_data_list = []
        #mtt_zero_list = [] 

        if self.runStep['pointsource']:
            # Now ptsrc is a list, which we need to get into movingTargetInputs...
            mtt_ptsrc = self.movingTargetInputs(self.params['simSignals']['pointsource'],'pointSource',MT_tracking=True,tracking_ra_vel=self.ra_vel,tracking_dec_vel=self.dec_vel,trackingPixVelFlag=vel_flag)
            #mtt_ptsrc_zero = mtt_ptsrc[0,:,:]
            mtt_data_list.append(mtt_ptsrc)
            #mtt_zero_list.append(mtt_ptsrc_zero)
            print("Done with creating moving targets from {}".format(self.params['simSignals']['pointsource']))

        if self.runStep['galaxies']:
            mtt_galaxies = self.movingTargetInputs(self.params['simSignals']['galaxyListFile'],'galaxies',MT_tracking=True,tracking_ra_vel=self.ra_vel,tracking_dec_vel=self.dec_vel,trackingPixVelFlag=vel_flag)

            mtt_data_list.append(mtt_galaxies)
            print("Done with creating moving targets from {}".format(self.params['simSignals']['galaxyListFile']))

        if self.runStep['extendedsource']:
            mtt_ext = self.movingTargetInputs(self.params['simSignals']['extended'],'extended',MT_tracking=True,tracking_ra_vel=self.ra_vel,tracking_dec_vel=self.dec_vel,trackingPixVelFlag=vel_flag)
                
            mtt_data_list.append(mtt_ext)

        # Add in the other objects which are not being tracked on 
        # (i.e. the sidereal targets)
        if len(mtt_data_list) > 0:
            for i in xrange(len(mtt_data_list)):
                non_sidereal_ramp += mtt_data_list[i]
                #non_sidereal_zero += mtt_zero_list[i]

        return non_sidereal_ramp#, non_sidereal_zero

    
    def readMTFile(self,file):
        #read in moving target list file
        mtlist = ascii.read(file,comment='#')

        #convert all columns to floats
        for col in mtlist.colnames:
            if mtlist[col].dtype in ['int64','int']:
                mtlist[col] = mtlist[col].data * 1.
                
        #check to see whether the position is in x,y or ra,dec
        pixelflag = False
        try:
            if 'position_pixels' in mtlist.meta['comments'][0]:
                pixelflag = True
        except:
            pass

        #if present, check whether the velocity entries are pix/sec
        #or arcsec/sec.
        pixelvelflag = False
        try:
            if 'velocity_pixels' in mtlist.meta['comments'][1]:
                pixelvelflag = True
        except:
            pass

        return mtlist,pixelflag,pixelvelflag

        
    def movingTargetInputs(self,file,input_type,MT_tracking=False,tracking_ra_vel=None,tracking_dec_vel=None,trackingPixVelFlag=False):

        #read in listfile of moving targets and perform needed calculations to get inputs
        #for moving_targets.py

        #input_type can be 'pointSource','galaxies', or 'extended'

        #get countrate for mag 15 source, for scaling later
        try:
            if self.params['Readout']['pupil'][0].upper() == 'F':
                usephot = 'pupil'
            else:
                usephot = 'filter'
            mag15rate = self.countvalues[self.params['Readout'][usephot]]
        except:
            print("Unable to find mag 15 countrate for {} filter in {}.".format(self.params['Readout'][usephot],self.params['Reffiles']['phot']))
            print("Fix me!")
            sys.exit()
            
        #read input file - should be able to use for all modes
        mtlist,pixelFlag,pixvelflag = self.readMTFile(file)
        
        if MT_tracking == True:
            # Here, we are tracking a non-sidereal target.
            try:
                # If there are moving targets on top of the non-
                # sidereal tracking (e.g. tracking Io but Europa
                # comes into the fov), then the velocity vector
                # of the moving target needs to be adjusted.
                # If the input catalog already contains
                # 'x_or_RA_velocity' then we know we have a moving
                # target. If it doesn't, then we have sidereal
                # targets, and we can simply set their velocity
                # as the inverse of that being tracked.
                mtlist['x_or_RA_velocity'] -= tracking_ra_vel * (1./365.25/24.)
                mtlist['y_or_Dec_velocity'] -= tracking_dec_vel * (1./365.25/24.)
                pixvelflag = trackingPixVelFlag
            except:
                print('Setting velocity of targets equal to the non-sidereal tracking velocity')
                mtlist['x_or_RA_velocity'] = 0. - tracking_ra_vel * (1./365.25/24.)
                mtlist['y_or_Dec_velocity'] = 0. - tracking_dec_vel * (1./365.25/24.)
                pixvelflag = trackingPixVelFlag

        #get necessary information for coordinate transformations
        coord_transform = None
        if self.runStep['astrometric']:

            #Read in the CRDS-format distortion reference file
            with AsdfFile.open(self.params['Reffiles']['astrometric']) as dist_file:
                coord_transform = dist_file.tree['model']

        #Using the requested RA,Dec of the reference pixel, along with the 
        #V2,V3 of the reference pixel, and the requested roll angle of the telescope
        #create a matrix that can be used to translate between V2,V3 and RA,Dec
        #for any pixel.
        #v2,v3 need to be in arcsec, and RA, Dec, and roll all need to be in degrees
        #attitude_matrix = rotations.attitude(self.refpix_pos['v2'],self.refpix_pos['v3'],self.ra,self.dec,self.params['Telescope']["rotation"])
        attitude_matrix = self.getAttitudeMatrix()
        
        #exposure times for all frames
        frameexptimes = self.frametime * np.arange(-1,self.params['Readout']['ngroup'] * (self.params['Readout']['nframe'] + self.params['Readout']['nskip']))
            
        #output image dimensions
        #dims = np.array(self.dark.data[0,0,:,:].shape)
        dims = self.nominal_dims
        newdimsx = np.int(dims[1] * self.coord_adjust['x'])
        newdimsy = np.int(dims[0] * self.coord_adjust['y'])

        mt_integration = np.zeros((len(frameexptimes)-1,newdimsy,newdimsx))

        for entry in mtlist:

            #for each object, calculate x,y or ra,dec of initial position
            pixelx,pixely,ra,dec,ra_str,dec_str = self.getPositions(entry['x_or_RA'],entry['y_or_Dec'],attitude_matrix,coord_transform,pixelFlag)

            #now generate a list of x,y position in each frame
            if pixvelflag == False:
                #calculate the RA,Dec in each frame
                #input velocities are arcsec/hour. ra/dec are in units of degrees,
                #so divide velocities by 3600^2.
                ra_frames = ra + (entry['x_or_RA_velocity']/3600./3600.) * frameexptimes
                dec_frames = dec + (entry['y_or_Dec_velocity']/3600./3600.) * frameexptimes

                #print('RA,Dec velocities per sec are: {},{}'.format(entry['x_or_RA_velocity']/3600.,entry['y_or_Dec_velocity']/3600.))
                #print('RA frames {}'.format(ra_frames))
                #print('Dec frames {}'.format(dec_frames))

                x_frames = []
                y_frames = []
                for in_ra,in_dec in zip(ra_frames,dec_frames):
                    #calculate the x,y position at each frame
                    px,py,pra,pdec,pra_str,pdec_str = self.getPositions(in_ra,in_dec,attitude_matrix,coord_transform,False)
                    x_frames.append(px)
                    y_frames.append(py)
                x_frames = np.array(x_frames)
                y_frames = np.array(y_frames)
                    
            else:
                #if input velocities are pixels/hour, then generate the list of
                #x,y in each frame directly
                x_frames = pixelx + (entry['x_or_RA_velocity']/3600.) * frameexptimes
                y_frames = pixely + (entry['y_or_Dec_velocity']/3600.) * frameexptimes
                      
                      
            #If the target never falls on the detector, then move on to the next target
            xfdiffs = np.fabs(x_frames - (newdimsx/2))
            yfdiffs = np.fabs(y_frames - (newdimsy/2))
            if np.min(xfdiffs) > (newdimsx/2) or np.min(yfdiffs) > (newdimsy/2):
                continue

            #So now we have xinit,yinit, a list of x,y positions for each frame, and the frametime.
            #subsample factor can be hardwired for now. outx and outy are also known. So all we need is the stamp
            #image, then we can call moving_targets.py and feed it these things, which contain all the info needed

            if input_type == 'pointSource':
                stamp = self.centerpsf

            elif input_type == 'extended':
                stamp,header = self.basicGetImage(entry['filename'])
                if entry['pos_angle'] != 0.:
                    stamp = self.basicRotateImage(stamp,entry['pos_angle'])

                #convolve with NIRCam PSF if requested
                if self.params['simSignals']['PSFConvolveExtended']:
                    stamp = s1.fftconvolve(stamp,self.centerpsf,mode='same')

            elif input_type == 'galaxies':
                stamp = self.create_galaxy(entry['radius'],entry['ellipticity'],entry['sersic_index'],entry['pos_angle'],1.) #entry['counts_per_frame_e'])
                #convolve the galaxy with the NIRCam PSF
                stamp = s1.fftconvolve(stamp,self.centerpsf,mode='same')


            #normalize the PSF to a total signal of 1.0
            totalsignal = np.sum(stamp)
            stamp /= totalsignal

            #Scale the stamp image to the requested magnitude
            scale = 10.**(0.4*(15.0-entry['magnitude']))
            rate = scale * mag15rate
            stamp *= rate

            #each entry will have stamp image as array, ra_init,dec_init,ra_velocity,dec_velocity,frametime,numframes,subsample_factor,outputarrayxsize,outputarrayysize (maybe without the values that will be the same to each entry.

            #entryList = (stamp,ra,dec,entry[3]/3600.,entry[4]/3600.,self.frametime,numframes,subsample_factor,outx,outy)
            #entryList = (stamp,x_frames,y_frames,self.frametime,subsample_factor,outx,outy)
            mt = moving_targets.MovingTarget()
            mt.subsampx = 3
            mt.subsampy = 3
            mt_integration += mt.create(stamp,x_frames,y_frames,self.frametime,newdimsx,newdimsy)

        #save the moving target input integration
        #if self.params['Output']['save_intermediates']:
        #    h0 = fits.PrimaryHDU(mt_integration)
        #    hl = fits.HDUList([h0])
        #    mtoutname = self.params['Output']['file'][0:-5] + '_movingTargetIntegration.fits'
        #    hl.writeto(mtoutname,overwrite=True)
        #    print("Integration showing only moving targets saved to {}".format(mtoutname))

        return mt_integration

    
    def getPositions(self,inx,iny,matrix,transform,pixelflag):
        #input a row containing x,y or ra,dec values, and figure out
        #x,y, RA, Dec, and RA string and Dec string 
        try:
            entry0 = float(inx)
            entry1 = float(iny)
            if not pixelflag:
                ra_str,dec_str = self.makePos(entry0,entry1)
                ra = entry0
                dec = entry1
        except:
            #if inputs can't be converted to floats, then 
            #assume we have RA/Dec strings. Convert to floats.
            ra_str = inx
            dec_str = iny
            ra,dec = self.parseRADec(ra_str,dec_str)

        #Case where point source list entries are given with RA and Dec
        if not pixelflag:

            #If distortion is to be included - either with or without the full set of coordinate
            #translation coefficients
            if self.runStep['astrometric']:
                pixelx,pixely = self.RADecToXY_astrometric(ra,dec,matrix,transform)
            else:
                #No distortion at all - "manual mode"
                pixelx,pixely = self.RADecToXY_manual(ra,dec)

        else:
            #Case where the point source list entry locations are given in units of pixels
            #In this case we have the source position, and RA/Dec are calculated only so 
            #they can be written out into the output source list file.

            pixelx = entry0
            pixely = entry1

            ra,dec,ra_str,dec_str = self.XYToRADec(pixelx,pixely,matrix,transform)
        return pixelx,pixely,ra,dec,ra_str,dec_str


    
    def nonsidereal_CRImage(self,file):
        #create countrate image of non-sidereal sources that are being tracked.

        #get countrate for mag 15 source, for scaling later
        try:
            if self.params['Readout']['pupil'][0].upper() == 'F':
                usephot = 'pupil'
            else:
                usephot = 'filter'
            mag15rate = self.countvalues[self.params['Readout'][usephot]]
        except:
            print("Unable to find mag 15 countrate for {} filter in {}.".format(self.params['Readout'][usephot],self.params['Reffiles']['phot']))
            sys.exit()


        totalCRList = []

        #read in file containing targets
        targs,pixFlag,velFlag = self.readMTFile(self.params['simSignals']['movingTargetToTrack'])
        
        #We need to keep track of the proper motion of the target being tracked, because
        #all other objects in the field of view will be moving at the same rate in the 
        #opposite direction
        track_ra_vel = targs[0]['x_or_RA_velocity']
        track_dec_vel = targs[0]['y_or_Dec_velocity']

        #sort the targets by whether they are point sources, galaxies, extended
        ptsrc_rows = []
        galaxy_rows = []
        extended_rows = []
        for i,line in enumerate(targs):
            if 'point' in line['object'].lower():
                ptsrc_rows.append(i)
            elif 'sersic' in line['object'].lower():
                galaxy_rows.append(i)
            else:
                extended_rows.append(i)

        #then re-use functions for the sidereal tracking situation
        if len(ptsrc_rows) > 0:
            ptsrc = targs[ptsrc_rows]
            if pixFlag:
                meta0 = 'position_pixel'
            else:
                meta0 = ''
            if velFlag:
                meta1 = 'velocity_pixel'
            else:
                meta1 = ''
            meta2 = 'Point sources with non-sidereal tracking. File produced by ramp_simulator.py'
            meta3 = 'from run using non-sidereal moving target list {}.'.format(self.params['simSignals']['movingTargetToTrack'])
            ptsrc.meta['comments'] = [meta0,meta1,meta2,meta3]
            ptsrc.write(os.path.join(self.params['Output']['directory'],'temp_non_sidereal_point_sources.list'),format='ascii',overwrite=True)

            ptsrc = self.getPointSourceList('temp_non_sidereal_point_sources.list')
            ptsrcCRImage = self.makePointSourceImage(ptsrc)
            totalCRList.append(ptsrcCRImage)

        if len(galaxy_rows) > 0:
            galaxies = targs[galaxy_rows]
            if pixFlag:
                meta0 = 'position_pixel'
            else:
                meta0 = ''
            if velFlag:
                meta1 = 'velocity_pixel'
            else:
                meta1 = ''
            meta2 = 'Galaxies (2d sersic profiles) with non-sidereal tracking. File produced by ramp_simulator.py'
            meta3 = 'from run using non-sidereal moving target list {}.'.format(self.params['simSignals']['movingTargetToTrack'])
            galaxies.meta['comments'] = [meta0,meta1,meta2,meta3]
            galaxies.write(os.path.join(self.params['Output']['directory'],'temp_non_sidereal_sersic_sources.list'),format='ascii',overwrite=True)
            
            #read in the PSF file with centered source
            wfe = self.params['simSignals']['psfwfe']
            wfegroup = self.params['simSignals']['psfwfegroup']
            basename = self.params['simSignals']['psfbasename'] + '_'
            if wfe == 0:
                psfname=basename+self.params['simSignals'][usefilt].lower()+'_zero'
            else:
                psfname=basename+self.params['Readout'][usefilt].lower()+"_"+str(wfe)+"_"+str(wfegroup)

            galaxyCRImage = self.makeGalaxyImage('temp_non_sidereal_sersic_sources.list',self.centerpsf)
            totalCRList.append(galaxyCRImage)

        if len(extended_rows) > 0:
            extended = targs[extended_rows]

            if pixFlag:
                meta0 = 'position_pixel'
            else:
                meta0 = ''
            if velFlag:
                meta1 = 'velocity_pixel'
            else:
                meta1 = ''
            meta2 = 'Extended sources with non-sidereal tracking. File produced by ramp_simulator.py'
            meta3 = 'from run using non-sidereal moving target list {}.'.format(self.params['simSignals']['movingTargetToTrack'])
            extended.meta['comments'] = [meta0,meta1,meta2,meta3]
            extended.write(os.path.join(self.params['Output']['directory'],'temp_non_sidereal_extended_sources.list'),format='ascii',overwrite=True)
            
            #read in the PSF file with centered source
            wfe = self.params['simSignals']['psfwfe']
            wfegroup = self.params['simSignals']['psfwfegroup']
            basename = self.params['simSignals']['psfbasename'] + '_'
            if wfe == 0:
                psfname=basename+self.params['simSignals'][usefilt].lower()+'_zero'
            else:
                psfname=basename+self.params['Readout'][usefilt].lower()+"_"+str(wfe)+"_"+str(wfegroup)

            extlist,extstamps = self.getExtendedSourceList('temp_non_sidereal_extended_sources.list')

            #translate the extended source list into an image
            extCRImage = self.makeExtendedSourceImage(extlist,extstamps)

            #if requested, convolve the stamp images with the NIRCam PSF
            if self.params['simSignals']['PSFConvolveExtended']:
                extCRImage = s1.fftconvolve(extCRImage,self.centerpsf,mode='same')

            totalCRList.append(extCRImage)

        #Now combine into a final countrate image of non-sidereal sources (that are being tracked)
        if len(totalCRList) > 0:
            totalCRImage = totalCRList[0]
            if len(totalCRList) > 1:
                for i in xrange(1,len(totalCRList)):
                    totalCRImage += totalCRList[i]
        else:
            print("No non-sidereal countrate targets produced.")
            print("You shouldn't be here.")
            sys.exit()

        return totalCRImage,track_ra_vel,track_dec_vel,velFlag
        

            
    def addedSignals(self):
        # Generate a signal rate image from input sources
        if self.params['Output']['grism_source_image'] == False:
            signalimage = np.zeros(self.nominal_dims)
            segmentation_map = np.zeros(self.nominal_dims)
        else:
            xd = np.int(self.nominal_dims[1] * self.coord_adjust['x'])
            yd = np.int(self.nominal_dims[0] * self.coord_adjust['y'])
            signalimage = np.zeros((yd,xd),dtype=np.float)
            segmentation_map = np.zeros((yd,xd))

        #yd,xd = signalimage.shape
        arrayshape = signalimage.shape

        #MASK IMAGE
        #Create a mask so that we don't add signal to masked pixels
        #Initially this includes only the reference pixels
        #Keep the mask image equal to the true subarray size, since this 
        #won't be used to make a requested grism source image
        maskimage = np.zeros((self.ffsize,self.ffsize),dtype=np.int)
        maskimage[4:self.ffsize-4,4:self.ffsize-4] = 1.
        
        #crop the mask to match the requested output array
        if self.params['Readout']['array_name'] != "FULL":
            maskimage = maskimage[self.subarray_bounds[1]:self.subarray_bounds[3]+1,self.subarray_bounds[0]:self.subarray_bounds[2]+1]

        #POINT SOURCES
        #Read in the list of point sources to add
        #Adjust point source locations using astrometric distortion
        #Translate magnitudes to counts in a single frame
        if self.runStep['pointsource'] == True:
            pslist = self.getPointSourceList(self.params['simSignals']['pointsource'])

            #translate the point source list into an image
            psfimage, ptsrc_segmap = self.makePointSourceImage(pslist)
            
            #save the point source image for examination by user
            if self.params['Output']['save_intermediates'] == True:
                psfImageName = self.basename + '_pointSourceRateImage_elec_per_sec.fits'
                h0 = fits.PrimaryHDU(psfimage)
                h1 = fits.ImageHDU(ptsrc_segmap)
                hlist = fits.HDUList([h0,h1])
                hlist.writeto(psfImageName,overwrite=True)
                #ptsrc_seg_file = psfImageName[0:-5] + '_segmap.fits'
                #self.saveSingleFits(psfimage,psfImageName)
                print("Point source image and segmap saved as {}".format(psfImageName))
                #self.saveSingleFits(ptsrc_segmap,ptsrc_seg_file)
                #print("Point source segmentation map saved as {}".format(ptsrc_seg_file))
                
            #Add the point source image to the overall image
            signalimage = signalimage + psfimage
            segmentation_map += ptsrc_segmap
            
        #Simulated galaxies
        #Read in the list of galaxy positions/magnitudes to simulate
        #and create a countrate image of those galaxies.
        if self.runStep['galaxies'] == True:
            galaxyCRImage, galaxy_segmap = self.makeGalaxyImage(self.params['simSignals']['galaxyListFile'],self.centerpsf)
            # Check segmentation map values. If the index numbers overlap
            # between point sources and galaxies, then bump up the values
            # in the galaxy list, and tell the user
            mingalseg = np.min(galaxy_segmap > 0)
            maxseg = np.max(segmentation_map)
            if mingalseg < maxseg:
                diff = maxseg - mingalseg
                print("WARNING: index numbers in {} overlap".format(self.params['simSignals']['galaxyListFile']))
                print("with those from {}".format(self.params['simSignals']['pointsource']))
                print("Adding {} to the galaxy object indexes.".format(diff+1))
                galaxy_segmap[galaxy_segmap > 0] += (diff+1)
            # Add galaxy segmentation map to the master copy
            segmentation_map += galaxy_segmap
                
            #save the galaxy image for examination by the user
            if self.params['Output']['save_intermediates'] == True:
                galImageName = self.basename + '_galaxyRateImage_elec_per_sec.fits'
                h0 = fits.PrimaryHDU(galaxyCRImage)
                h1 = fits.ImageHDU(galaxy_segmap)
                hlist = fits.HDUList([h0,h1])
                hlist.writeto(galImageName,overwrite=True)
                #self.saveSingleFits(galaxyCRImage,galImageName)
                print("Simulated galaxy image and segmap saved as {}".format(galImageName))

            #add the galaxy image to the signalimage
            signalimage = signalimage + galaxyCRImage

        #read in extended signal image and add the image to the overall image
        if self.runStep['extendedsource'] == True:
            #print("Reading in extended image from {} and adding to the simulated data rate image.".format(self.params['simSignals']['extended']))

            extlist,extstamps = self.getExtendedSourceList(self.params['simSignals']['extended'])

            #translate the extended source list into an image
            extimage, ext_segmap = self.makeExtendedSourceImage(extlist,extstamps)

            # Check segmentation map values. If the index numbers overlap
            # between master copy and extended sources, then bump up the values
            # in the extended source list, and tell the user
            minextseg = np.min(ext_segmap > 0)
            maxseg = np.max(segmentation_map)
            if minextseg < maxseg:
                diff = maxseg - minextseg
                print("WARNING: index numbers in {} overlap".format(self.params['simSignals']['extended']))
                print("with those from previously added sources.")
                print("Adding {} to the extended object indexes.".format(diff+1))
                ext_segmap[ext_segmap > 0] += (diff+1)
            # Add galaxy segmentation map to the master copy
            segmentation_map += ext_segmap

            # Save extended source image and segmap 
            if self.params['Output']['save_intermediates'] == True:
                extImageName = self.basename + '_extendedObject_elec_per_sec.fits'
                h0 = fits.PrimaryHDU(extimage)
                h1 = fits.ImageHDU(ext_segmap)
                hlist = fits.HDUList([h0,h1])
                hlist.writeto(extImageName,overwrite=True)
                print("Extended object image and segmap saved as {}".format(extImageName))

            #convolution now done inside makeextendedsourceimage
            #if requested, convolve the stamp images with the NIRCam PSF
            #if self.params['simSignals']['PSFConvolveExtended']:
            #    extimage = s1.fftconvolve(extimage,self.centerpsf,mode='same')

            #add the extended image to the synthetic signal rate image
            signalimage = signalimage + extimage

        #ZODIACAL LIGHT
        if self.runStep['zodiacal'] == True:
            zodiangle = self.eclipticangle() - self.params['Telescope']['rotation']
            zodiacalimage,zodiacalheader = self.getImage(self.params['simSignals']['zodiacal'],arrayshape,True,zodiangle,arrayshape/2)
            signalimage = signalimage + zodiacalimage*self.params['simSignals']['zodiscale']

        #SCATTERED LIGHT - no rotation here. 
        if self.runStep['scattered']:
            scatteredimage,scatteredheader = self.getImage(self.params['simSignals']['scattered'],arrayshape,False,0.0,arrayshape/2)
            signalimage = signalimage + scatteredimage*self.params['simSignals']['scatteredscale']

        #CONSTANT BACKGROUND
        signalimage = signalimage + self.params['simSignals']['bkgdrate']


        #Save the image containing all of the added sources from the 'sky'
        #if self.params['Output']['save_intermediates'] == True:
        #    sourcesImageName = self.params['Output']['file'][0:-5] + '_AddedSourcesRateImage_elec_per_sec.fits'
        #    self.saveSingleFits(signalimage,sourcesImageName)
        #    print("Image of added sources from the 'sky' saved as {}".format(sourcesImageName))


        #Save the final rate image of added signals
        if self.params['Output']['save_intermediates'] == True:
            #rateImageName = self.params['Output']['file'][0:-5] + '_AddedSourcesPlusDetectorEffectsRateImage_elec_per_sec.fits'
            rateImageName = self.basename + '_AddedSources_elec_per_sec.fits'
            self.saveSingleFits(signalimage,rateImageName)
            print("Signal rate image of all added sources (plus flats and IPC applied if requested) saved as {}".format(rateImageName))

        return signalimage, segmentation_map


    def getDistortionCoefficients(self,table,from_sys,to_sys,aperture):
        '''from the table of distortion coefficients, get the coeffs that correspond
        to the requested transformation and return as a list for x and another for y
        '''
        match = table['AperName'] == aperture
        if np.any(match) == False:
            print("Aperture name {} not found in input CSV file.".format(aperture))
            sys.exit()

        row = table[match]

        if ((from_sys == 'science') & (to_sys == 'ideal')):
            label = 'Sci2Idl'
        elif ((from_sys == 'ideal') & (to_sys == 'science')):
            label = 'Idl2Sci'
        else:
            print("WARNING: from_sys of {} and to_sys of {} not a valid transformation.".format(from_sys,to_sys))
            sys.exit()
        
        #get the coefficients, return as list
        X_cols = [c for c in row.colnames if label+'X' in c]
        Y_cols = [c for c in row.colnames if label+'Y' in c]
        x_coeffs = [row[c].data[0] for c in X_cols]
        y_coeffs = [row[c].data[0] for c in Y_cols]

        #Also get the V2,V3 values of the reference pixel
        v2ref = row['V2Ref'].data[0]
        v3ref = row['V3Ref'].data[0]

        #Get parity and V3 Y angle info 
        parity = row['VIdlParity'].data[0]
        yang = row['V3IdlYAngle'].data[0]
        v3scixang = row['V3SciXAngle'].data[0]
        
        #Get pixel scale info - not used but needs to be in output
        xsciscale = row['XSciScale'].data[0]
        ysciscale = row['YSciScale'].data[0]
        
        return x_coeffs,y_coeffs,v2ref,v3ref,parity,yang,xsciscale,ysciscale,v3scixang

    
    def getPointSourceList(self,filename):
        #read in the list of point sources to add, and adjust the
        #provided positions for astrometric distortion

        #find the array sizes of the PSF files in the library. Assume they are all the same.
        #We want the distance from the PSF peak to the edge, assuming the peak is centered
        if self.params['simSignals']['psfwfe'] != 0:
            numstr = str(self.params['simSignals']['psfwfe'])
        else:
            numstr = 'zero'
        psflibfiles = glob.glob(self.params['simSignals']['psfpath'] +'*')


        #If a PSF library is specified, then just get the dimensions from one of the files
        if self.params['simSignals']['psfpath'] != None:
            h = fits.open(psflibfiles[0])
            edgex = h[0].header['NAXIS1'] / 2 - 1
            edgey = h[0].header['NAXIS2'] / 2 - 1
            self.psfhalfwidth = np.array([edgex,edgey])
            h.close()
        else:
            #if no PSF library is specified, then webbpsf will be creating the PSF on the 
            #fly. In this case, we assume webbpsf's default output size of 301x301 pixels?
            edgex = int(301 / 2)
            edgey = int(301 / 2)
            print("INFO: no PSF library specified, but point sources are to be added to")
            print("the output. PSFs will be generated by WebbPSF on the fly")
            print("Not yet implemented.")
            sys.exit()


        pointSourceList = Table(names=('index','pixelx','pixely','RA','Dec','RA_degrees','Dec_degrees','magnitude','countrate_e/s','counts_per_frame_e'),dtype=('i','f','f','S14','S14','f','f','f','f','f'))
        
        try:
            lines,pixelflag = self.readPointSourceFile(filename)
            if pixelflag:
                print("Point source list input positions assumed to be in units of pixels.")
            else:
                print("Point list input positions assumed to be in units of RA and Dec.") 
        except:
            print("WARNING: Unable to open the point source list file {}".format(filename))
            sys.exit()

        #File to save adjusted point source locations
        psfile = self.params['Output']['file'][0:-5] + '_pointsources.list'
        pslist = open(psfile,'w')
    
        # If the input catalog has an index column
        # use that, otherwise add one
        if 'index' in lines.colnames:
            print('Using point source catalog index numbers')
            indexes = lines['index']
        else:
            print('No point source catalog index numbers. Adding to output: {}.'.format(psfile))
            indexes = np.arange(len(lines['x_or_RA']))
            
        dtor = math.radians(1.)
        nx = (self.subarray_bounds[2]-self.subarray_bounds[0])+1
        ny = (self.subarray_bounds[3]-self.subarray_bounds[1])+1
        xc = (self.subarray_bounds[2]+self.subarray_bounds[0])/2.
        yc = (self.subarray_bounds[3]+self.subarray_bounds[1])/2.

        #Location of the subarray's reference pixel. 
        xrefpix = self.refpix_pos['x']
        yrefpix = self.refpix_pos['y']

        # center positions, sub-array sizes in pixels
        # now offset the field center to array center for astrometric distortion corrections
        coord_transform = None
        if self.runStep['astrometric']:

            #Read in the CRDS-format distortion reference file
            with AsdfFile.open(self.params['Reffiles']['astrometric']) as dist_file:
                coord_transform = dist_file.tree['model']

        #Using the requested RA,Dec of the reference pixel, along with the 
        #V2,V3 of the reference pixel, and the requested roll angle of the telescope
        #create a matrix that can be used to translate between V2,V3 and RA,Dec
        #for any pixel.
        #v2,v3 need to be in arcsec, and RA, Dec, and roll all need to be in degrees
        attitude_matrix = self.getAttitudeMatrix()
             
        #Define the min and max source locations (in pixels) that fall onto the subarray
        #Inlude the effects of a requested grism_direct image, and also keep sources that
        #will only partially fall on the subarray
        #pixel coords here can still be negative and kept if the grism image is being made

        #First, coord limits for just the subarray
        miny = 0
        maxy = self.subarray_bounds[3] - self.subarray_bounds[1] 
        minx = 0
        maxx = self.subarray_bounds[2] - self.subarray_bounds[0] 
        #print('before adjusting for grism, min/max x {} {}, min/max y {} {}'.format(minx,maxx,miny,maxy))
        
        #Expand the limits if a grism direct image is being made
        if self.params['Output']['grism_source_image'] == True:
            extrapixy = np.int((maxy+1)/2 * (self.coord_adjust['y'] - 1.))
            miny -= extrapixy
            maxy += extrapixy
            extrapixx = np.int((maxx+1)/2 * (self.coord_adjust['x'] - 1.))
            minx -= extrapixx
            maxx += extrapixx

        #Now, expand the dimensions again to include point sources that fall only partially on the 
        #subarray
        miny -= edgey
        maxy += edgey
        minx -= edgex
        maxx += edgex

        #Write out the RA and Dec of the field center to the output file
        #Also write out column headers to prepare for source list
        pslist.write("# Field center (degrees): %13.8f %14.8f y axis rotation angle (degrees): %f  image size: %4.4d %4.4d\n" % (self.ra,self.dec,self.params['Telescope']['rotation'],nx,ny))
        pslist.write('#\n')
        pslist.write("#    Index   RA_(hh:mm:ss)   DEC_(dd:mm:ss)   RA_degrees      DEC_degrees     pixel_x   pixel_y    magnitude   counts/sec    counts/frame\n")

        #Loop over input lines in the source list 
        for index,values in zip(indexes,lines):
            #try:
            #line below (if 1>0) used to keep the block of code below at correct indent for the try: above
            #the try: is commented out for code testing.
            if 1>0:
                try:
                    entry0 = float(values['x_or_RA'])
                    entry1 = float(values['y_or_Dec'])
                    if not pixelflag:
                        ra_str,dec_str = self.makePos(entry0,entry1)
                        ra = entry0
                        dec = entry1
                except:
                    #if inputs can't be converted to floats, then 
                    #assume we have RA/Dec strings. Convert to floats.
                    ra_str = values['x_or_RA']
                    dec_str = values['y_or_Dec']
                    ra,dec = self.parseRADec(ra_str,dec_str)

                #Case where point source list entries are given with RA and Dec
                if not pixelflag:

                    #If distortion is to be included - either with or without the full set of coordinate
                    #translation coefficients
                    if self.runStep['astrometric']:
                        pixelx,pixely = self.RADecToXY_astrometric(ra,dec,attitude_matrix,coord_transform)
                    else:
                        #No distortion at all - "manual mode"
                        pixelx,pixely = self.RADecToXY_manual(ra,dec)

                else:
                    #Case where the point source list entry locations are given in units of pixels
                    #In this case we have the source position, and RA/Dec are calculated only so 
                    #they can be written out into the output source list file.

                    #Assume that the input x and y values are coordinate values
                    #WITHIN THE SPECIFIED SUBARRAY. So for example, a source in the file
                    #at 0,0 when you are making a SUB160 ramp will fall on the lower left
                    #corner of the SUB160 subarray, NOT the lower left corner of the full
                    #frame. 

                    pixelx = entry0
                    pixely = entry1

                    ra,dec,ra_str,dec_str = self.XYToRADec(pixelx,pixely,attitude_matrix,coord_transform)
                            

                #Get the input magnitude of the point source
                mag=float(values['magnitude'])

                if pixely > miny and pixely < maxy and pixelx > minx and pixelx < maxx:
                            
                    #set up an entry for the output table
                    entry = [index,pixelx,pixely,ra_str,dec_str,ra,dec,mag]

                    #translate magnitudes to countrate
                    scale = 10.**(0.4*(15.0-mag))

                    #get the countrate that corresponds to a 15th magnitude star for this filter
                    if self.params['Readout']['pupil'][0].upper() == 'F':
                        usefilt = 'pupil'
                    else:
                        usefilt = 'filter'
                    cval = self.countvalues[self.params['Readout'][usefilt]]

                    #DEAL WITH THIS LATER, ONCE PYSYNPHOT IS INCLUDED WITH PIPELINE DIST?
                    if cval == 0:
                        print("Countrate value for {} is zero in {}.".format(self.params['Readout'][usefilt],self.parameters['phot_file']))
                        print("Eventually attempting to calculate value using pysynphot.")
                        print("but pysynphot is not present in jwst build 6, so pushing off to later...")
                        sys.exit()
                        cval = self.findCountrate(self.params['Readout'][usefilt])

                    #translate to counts in single frame at requested array size
                    framecounts = scale*cval*self.frametime
                    countrate = scale*cval

                    #add the countrate and the counts per frame to pointSourceList
                    #since they will be used in future calculations
                    entry.append(countrate)
                    entry.append(framecounts)

                    #add the good point source, including location and counts, to the pointSourceList
                    pointSourceList.add_row(entry)

                    #write out positions, distances, and counts to the output file
                    pslist.write("%i %s %s %14.8f %14.8f %9.3f %9.3f  %9.3f  %13.6e   %13.6e\n" % (index,ra_str,dec_str,ra,dec,pixelx,pixely,mag,countrate,framecounts))

        print("Number of point sources found within the requested aperture: {}".format(len(pointSourceList)))
        #close the output file
        pslist.close()
        
        #If no good point sources were found in the requested array, alert the user
        if len(pointSourceList) < 1:
            print("INFO: no point sources within the requested array.")
            #print("The point source image option is being turned off")
            #self.runStep['pointsource']=False
            #if self.runStep['extendedsource'] == False and self.runStep['cosmicray'] == False:
            #    print("Error: no input point sources, extended image, nor cosmic rays specified")
            #    print("Exiting...")
            #    sys.exit()
        
        return pointSourceList

    
    def makePointSourceImage(self,pointSources):
        dims = np.array(self.nominal_dims)
        #dims = np.array(self.dark.data[0,0,:,:].shape)

        #offset that needs to be applied to the x,y positions of the
        #source list to account for case where we make a point
        #source image that is extra-large, to be used as a grism 
        #direct image
        deltax = 0
        deltay = 0

        newdimsx = np.int(dims[1] * self.coord_adjust['x'])
        newdimsy = np.int(dims[0] * self.coord_adjust['y'])
        deltax = self.coord_adjust['xoffset']
        deltay = self.coord_adjust['yoffset']
        dims = np.array([newdimsy,newdimsx])

        #create the empty image
        psfimage = np.zeros((dims[0],dims[1]))

        # create empty segmentation map
        seg = segmap.SegMap()
        seg.xdim = newdimsx
        seg.ydim = newdimsy
        seg.initialize_map()
        
        #Loop over the entries in the point source list
        for entry in pointSources:
            #adjust x,y position if the grism output image is requested
            xpos = entry['pixelx'] + deltax
            ypos = entry['pixely'] + deltay

            #desired counts per second in the point source
            counts = entry['countrate_e/s'] #/ self.frametime

            #find sub-pixel offsets in position from the center of the pixel
            xoff = math.floor(xpos)
            yoff = math.floor(ypos)
            xfract = abs(xpos-xoff)
            yfract = abs(ypos-yoff)

            #Now we need to determine the proper PSF file to read in from the library
            #This depends on the sub-pixel offsets above
            interval = self.params['simSignals']['psfpixfrac']
            numperpix = int(1./interval)
            a = interval * int(numperpix*xfract + 0.5) - 0.5
            b = interval * int(numperpix*yfract + 0.5) - 0.5

            if a < 0:
                astr = str(a)[0:4]
            else:
                astr = str(a)[0:3]
            if b < 0:
                bstr = str(b)[0:4]
            else:
                bstr = str(b)[0:3]

            #generate the psf file name based on the center of the point source
            #in units of fraction of a pixel
            frag = astr + '_' + bstr
            frag = frag.replace('-','m')
            frag = frag.replace('.','p')
            
            #now create the PSF image. If no PSF library is supplied
            #then webbpsf will be called to create a PSF. In that case, return
            #zeros right now for the PSF
            if self.params['simSignals']['psfpath'] is None:
                webbpsfimage = self.psfimage
            else:
                #case where PSF library location is specified. 
                #Read in the appropriate PSF file
                try:
                    psffn = self.psfname+'_'+frag+'.fits'
                    webbpsfimage = fits.getdata(psffn)
                except:
                    print("ERROR: Could not load PSF file {} from library".format(psffn))
                    sys.exit()

            #Extract the appropriate subarray from the PSF image if necessary
            #Assume that the brightest pixel corresponds to the peak of the psf
            nyshift,nxshift = np.where(webbpsfimage == np.max(webbpsfimage))
            nyshift = nyshift[0]
            nxshift = nxshift[0]

            psfdims = webbpsfimage.shape
            nx = int(xoff)
            ny = int(yoff)
            i1 = max(nx-nxshift,0)
            i2 = min(nx+1+nxshift,dims[1])
            j1 = max(ny-nyshift,0)
            j2 = min(ny+1+nyshift,dims[0])
            k1 = nxshift-(nx-i1)
            k2 = nxshift+(i2-nx)
            l1 = nyshift-(ny-j1)
            l2 = nyshift+(j2-ny)

            #if the cutout for the psf is larger than
            #the psf array, truncate it, along with the array
            #in the source image where it will be placed
            if l2 > psfdims[0]:
                l2 = psfdims[0]
                j2 = j1 + (l2-l1)

            if k2 > psfdims[1]:
                k2 = psfdims[1]
                i2 = i1 + (k2-k1)

            #At this point coordinates are in the final output array coordinate system, so there
            #should be no negative values, nor values larger than the output array size
            if j1 < 0 or i1<0 or l1<0 or k1<0:
                print(j1,i1,l1,k1)
                print('bad low')
            if j2>(dims[0]+1) or i2>(dims[1]+1) or l2>(psfdims[1]+1) or k2>(psfdims[1]+1):
                print(j2,i2,l2,k2)
                print('bad high')

            try:
                psfimage[j1:j2,i1:i2] = psfimage[j1:j2,i1:i2] + webbpsfimage[l1:l2,k1:k2]*counts
                # Divide readnoise by 100 sec, which is a 10 group RAPID ramp?
                noiseval = self.single_ron / 100. + self.params['simSignals']['bkgdrate']
                if self.params['Inst']['mode'].lower() == 'wfss':
                    noiseval += self.grism_background
                seg.add_object_noise(webbpsfimage[l1:l2,k1:k2]*counts,j1,i1,entry['index'],noiseval)
            except:
                #In here we catch sources that are off the edge
                #of the detector. These may not necessarily be caught in
                #getpointsourcelist because if the PSF is not centered
                #in the webbpsf stamp, then the area to be pulled from
                #the stamp may shift off of the detector.
                #print(psffn,dims,psfdims,xoff,yoff,nx,ny,nxshift,nyshift)
                #print(j1,j2,i1,i2,l1,l2,k1,k2)
                #print(entry)
                #sys.exit()
                pass
                
        return psfimage,seg.segmap
    

    def cropPSF(self,psf):
        '''take an array containing a psf and crop it such that the brightest
        pixel is in the center of the array'''
        nyshift,nxshift = np.where(psf == np.max(psf))
        nyshift = nyshift[0]
        nxshift = nxshift[0]
        py,px = psf.shape

        xl = nxshift - 0
        xr = px - nxshift - 1
        if xl <= xr:
            xdist = xl
        if xr < xl:
            xdist = xr

        yl = nyshift - 0
        yr = py - nyshift - 1
        if yl <= yr:
            ydist = yl
        if yr < yl:
            ydist = yr

        return psf[nyshift-ydist:nyshift+ydist+1,nxshift-xdist:nxshift+xdist+1]

    
    def readPointSourceFile(self,filename):
        # Read in the point source list
        try:
            #read table
            gtab = ascii.read(filename,comment='#')
            #print(gtab)
            #Look at the header lines to see if inputs
            #are in units of pixels or RA,Dec
            pflag = False
            try:
                if 'position_pixel' in gtab.meta['comments'][0]:
                    pflag = True
            except:
                pass

        except:
            print("WARNING: Unable to open the source list file {}".format(filename))
            sys.exit()

        return gtab,pflag

    
    def getAttitudeMatrix(self):
        #create an attitude_matrix from the distortion reference file model and other info
        #calculate a local roll angle for the aperture
        self.local_roll = set_telescope_pointing.compute_local_roll(self.params['Telescope']['rotation'],self.ra,self.dec,self.v2_ref,self.v3_ref)
        #create attitude_matrix    
        attitude_matrix = rotations.attitude(self.refpix_pos['v2'],self.refpix_pos['v3'],self.ra,self.dec,self.local_roll)
        return attitude_matrix


    def makePos(self,alpha1,delta1):
        #given a numerical RA/Dec pair, convert to string
        #values hh:mm:ss
        if alpha1 < 0.: 
            alpha1=alpha1+360.
        if delta1 < 0.: 
            sign="-"
            d1=abs(delta1)
        else:
            sign="+"
            d1=delta1
        decd=int(d1)
        value=60.*(d1-float(decd))
        decm=int(value)
        decs=60.*(value-decm)
        a1=alpha1/15.0
        radeg=int(a1)
        value=60.*(a1-radeg)
        ramin=int(value)
        rasec=60.*(value-ramin)
        alpha2="%2.2d:%2.2d:%7.4f" % (radeg,ramin,rasec)
        delta2="%1s%2.2d:%2.2d:%7.4f" % (sign,decd,decm,decs)
        alpha2=alpha2.replace(" ","0")
        delta2=delta2.replace(" ","0")
        return alpha2,delta2


    def parseRADec(self,rastr,decstr):
        #convert the input RA and Dec strings to floats
        try:
            rastr=rastr.lower()
            rastr=rastr.replace("h",":")
            rastr=rastr.replace("m",":")
            rastr=rastr.replace("s","")
            decstr=decstr.lower()
            decstr=decstr.replace("d",":")
            decstr=decstr.replace("m",":")
            decstr=decstr.replace("s","")

            values=rastr.split(":")
            ra0=15.*(int(values[0])+int(values[1])/60.+float(values[2])/3600.)

            values=decstr.split(":")
            if "-" in values[0]:
                sign=-1
                values[0]=values[0].replace("-"," ")
            else:
                sign=+1
            dec0=sign*(int(values[0])+int(values[1])/60.+float(values[2])/3600.)
            return ra0,dec0
        except:
            print("Error parsing RA,Dec strings: {} {}".format(rastr,decstr))
            sys.exit()


    def RADecToXY_astrometric(self,ra,dec,attitude_matrix,coord_transform):
        #Translate backwards, RA,Dec to V2,V3
        pixelv2,pixelv3 = rotations.getv2v3(attitude_matrix,ra,dec)

        if self.runStep['distortion_coeffs']:
            #If the full set of distortion coefficients are provided, then
            #use those to make the exact transformation from the 'ideal'
            #to 'science' coordinate systems

            #Now V2,V3 to undistorted angular distance from the reference pixel
            xidl = self.v2v32idlx(pixelv2-self.v2_ref,pixelv3-self.v3_ref)
            yidl = self.v2v32idly(pixelv2-self.v2_ref,pixelv3-self.v3_ref)
                                
            #Finally, undistorted distances to distorted pixel values
            deltapixelx, deltapixely, err, iter = polynomial.invert(self.x_sci2idl,self.y_sci2idl,xidl,yidl,5)

            pixelx = deltapixelx + self.refpix_pos['x']
            pixely = deltapixely + self.refpix_pos['y']

        else:
            #If the full set of distortion coefficients are not provided,
            #then we fall back to the coordinate transform provided by the
            #distortion reference file. These results are not exact, and
            #become less accurate the farther the source is from the center
            #of the detector. Results can be incorrect by ~20 pixels in the
            #corners of the detector.
    
            #Now go backwards from V2,V3 to distorted pixels
            #deltapixelx,deltapixely = coord_transform.inverse(pixelv2-self.refpix_pos['v2'],pixelv3-self.refpix_pos['v3'])
            pixelx,pixely = coord_transform.inverse(pixelv2,pixelv3)
            

        return pixelx,pixely


    def RADecToXY_manual(self,ra,dec):
        #In this case, the sources are provided as an RA,Dec list, 
        #but no astrometry information is provided. So assume an average
        #pixel scale and calculate the pixel position of the source from that.
        #This obviously does not include distortion, and is kind of a last
        #resort.
        ra_source = ra * 3600.
        dec_source = dec * 3600.
        
        dist_between,deltaang = self.dist([self.ra,self.dec],[ra_source,dec_source])

        #Now translate to deltax and deltay if the 
        #position angle is non-zero
        tot_ang = deltaang + (0. - self.params['Telescope']['rotation'] * np.pi / 180.)

        deltax = dist_between * np.sin(tot_ang) / self.pixscale[0]
        deltay = dist_between * np.cos(tot_ang) / self.pixscale[0]
                            
        pixelx = self.refpix_pos['x'] + deltax
        pixely = self.refpix_pos['y'] + deltay

        return pixelx,pixely


    def XYToRADec(self,pixelx,pixely,attitude_matrix,coord_transform):
        #Translate a given x,y location on the detector
        #to RA,Dec

        #If distortion is to be included
        #if self.runStep['astrometric']:
        if coord_transform is not None:
            #Transform distorted pixels to V2,V3
            #deltav2,deltav3 = coord_transform(pixelx-self.refpix_pos['x'],pixely-self.refpix_pos['y'])
            #pixelv2 = deltav2 + self.refpix_pos['v2']
            #pixelv3 = deltav3 + self.refpix_pos['v3']
            pixelv2,pixelv3 = coord_transform(pixelx,pixely)
            
            #Now translate V2,V3 to RA,Dec
            ra,dec = rotations.pointing(attitude_matrix,pixelv2,pixelv3)

        else:
            #Without including distortion. 
            #Fall back to "manual" calculations
            dist_between = np.sqrt((pixelx-self.refpix_pos['x'])**2 + (pixely-self.refpix_pos['y'])**2)
            deltaang = np.arctan2(pixely,pixelx)

            tot_ang = deltaang + (self.parms['Telescope']['rotation'] * np.pi / 180.)

            deltara = dist_between * np.sin(tot_ang) / self.pixoscale[0]
            deltadec = dist_between * np.cos(tot_ang) / self.pixscale[0]
                            
            ra = self.ra + deltara
            dec = self.dec + deltadec
                            
        #Translate the RA/Dec floats to strings
        ra_str,dec_str = self.makePos(ra,dec)

        return ra,dec,ra_str,dec_str


    def readGalaxyFile(self,filename):
        # Read in the galaxy source list
        try:
            #read table
            gtab = ascii.read(filename,comment='#')

            #Look at the header lines to see if inputs
            #are in units of pixels or RA,Dec
            pflag = False
            rpflag = False
            try:
                if 'position_pixel' in gtab.meta['comments'][0]:
                    pflag = True
            except:
                pass
            try:
                if 'radius_pixel' in gtab.meta['comments'][1]:
                    rpflag = True
            except:
                pass

        except:
            print("WARNING: Unable to open the galaxy source list file {}".format(filename))
            sys.exit()

        return gtab,pflag,rpflag


    def filterGalaxyList(self,galaxylist,pixelflag,radiusflag):
        #given a list of galaxies (location, size, orientation, magnitude)
        #keep only those which will fall fully or partially on the output array
        
        filteredList = Table(names=('index','pixelx','pixely','RA','Dec','RA_degrees','Dec_degrees','V2','V3','radius','ellipticity','pos_angle','sersic_index','magnitude','countrate_e/s','counts_per_frame_e'),dtype=('i','f','f','S14','S14','f','f','f','f','f','f','f','f','f','f','f'))

        #each entry in galaxylist is:
        #x_or_RA  y_or_Dec  radius  ellipticity  pos_angle  sersic_index  magnitude
        #remember that x/y are interpreted as coordinates in the output subarray
        #NOT full frame coordinates. This is the same as the point source list coords

        #First, begin to define the pixel limits beyond which a galaxy will be completely
        #outside of the field of view
        #First, coord limits for just the subarray
        miny = 0
        maxy = self.subarray_bounds[3] - self.subarray_bounds[1] 
        minx = 0
        maxx = self.subarray_bounds[2] - self.subarray_bounds[0] 
        ny = self.subarray_bounds[3] - self.subarray_bounds[1]
        nx = self.subarray_bounds[2] - self.subarray_bounds[0]
        #print('before adjusting for grism, min/max x {} {}, min/max y {} {}'.format(minx,maxx,miny,maxy))
        
        #Expand the limits if a grism direct image is being made
        if self.params['Output']['grism_source_image'] == True:
            extrapixy = np.int((maxy+1)/2 * (self.grism_direct_factor - 1.))
            miny -= extrapixy
            maxy += extrapixy
            extrapixx = np.int((maxx+1)/2 * (self.grism_direct_factor - 1.))
            minx -= extrapixx
            maxx += extrapixx

            nx = np.int(nx * self.grism_direct_factor)
            ny = np.int(ny * self.grism_direct_factor)

        #Create transform matrix for galaxy sources
        #Read in the CRDS-format distortion reference file
        coord_transform = None
        if self.runStep['astrometric']:
            with AsdfFile.open(self.params['Reffiles']['astrometric']) as dist_file:
                coord_transform = dist_file.tree['model']

        #Using the requested RA,Dec of the reference pixel, along with the 
        #V2,V3 of the reference pixel, and the requested roll angle of the telescope
        #create a matrix that can be used to translate between V2,V3 and RA,Dec
        #for any pixel
        #v2,v3 need to be in arcsec, and RA, Dec, and roll all need to be in degrees
        #attitude_matrix = rotations.attitude(self.refpix_pos['v2'],self.refpix_pos['v3'],self.ra,self.dec,self.params['Telescope']["rotation"])
        attitude_matrix = self.getAttitudeMatrix()

        # If an index column is present use that, otherwise
        # create one
        if 'index' in galaxylist.colnames:
            indexes = galaxylist['indexes']
        else:
            indexes = np.arange(len(galaxylist['radius']))

        
        #Loop over galaxy sources
        for index,source in zip(indexes,galaxylist):

            #If galaxy radii are given in units of arcseconds, translate to pixels 
            if radiusflag == False:
                source['radius'] /= self.pixscale[0]

            #how many pixels beyond the nominal subarray edges can a source be located and
            #still have it fall partially on the subarray? Galaxy stamps are nominally set to
            #have a length and width equal to 100 times the requested radius.
            edgex = source['radius'] * 100 / 2 - 1
            edgey = source['radius'] * 100 / 2 - 1

            #reset the field of view limits for the size of the current stamp image
            outminy = miny - edgey
            outmaxy = maxy + edgey
            outminx = minx - edgex
            outmaxx = maxx + edgex
        
            try:
                entry0 = float(source['x_or_RA'])
                entry1 = float(source['y_or_Dec'])
                if not pixelflag:
                    ra_str,dec_str = self.makePos(entry0,entry1)
                    #dec_str = self.makePos(entry0,entry1)
                    ra = entry0
                    dec = entry1
            except:
                #if inputs can't be converted to floats, then 
                #assume we have RA/Dec strings. Convert to floats.
                ra_str = source['x_or_RA']
                dec_str = source['y_or_Dec']
                ra,dec = self.parseRADec(ra_str,dec_str)

            #case where point source list entries are given with RA and Dec
            if not pixelflag:

                #if distortion is to be included
                if self.runStep['astrometric']:

                    pixelx,pixely = self.RADecToXY_astrometric(ra,dec,attitude_matrix,coord_transform)

                else:
                    #No distortion. Fall back to "manual" calculations
                    pixelx,pixely = self.RADecToXY_manual(ra,dec)
                    
            else:
                #case where the point source list entry locations are given in units of pixels
                #In this case we have the source position, and RA/Dec are calculated only so 
                #they can be written out into the output source list file.

                #Assume that the input x and y values are coordinate values
                #WITHIN THE SPECIFIED SUBARRAY. So for example, a source in the file
                #at 0,0 when you are making a SUB160 ramp will fall on the lower left
                #corner of the SUB160 subarray, NOT the lower left corner of the full
                #frame. 

                pixelx = entry0
                pixely = entry1

                ra,dec,ra_str,dec_str = self.XYToRADec(pixelx,pixely,attitude_matrix,coord_transform)

            #only keep the source if the peak will fall within the subarray
            if pixely > outminy and pixely < outmaxy and pixelx > outminx and pixelx < outmaxx:

                pixelv2,pixelv3 = rotations.getv2v3(attitude_matrix,ra,dec)
                entry = [index,pixelx,pixely,ra_str,dec_str,ra,dec,pixelv2,pixelv3,source['radius'],source['ellipticity'],source['pos_angle'],source['sersic_index']]

                #Now look at the input magnitude of the point source
                #append the mag and pixel position to the list of ra,dec
                mag = float(source['magnitude'])
                entry.append(mag)
    
                #translate magnitudes to countrate
                scale = 10.**(0.4*(15.0-mag))

                #get the countrate that corresponds to a 15th magnitude star for this filter
                if self.params['Readout']['pupil'][0].upper() == 'F':
                    usefilt = 'pupil'
                else:
                    usefilt = 'filter'
                cval = self.countvalues[self.params['Readout'][usefilt]]

                #DEAL WITH THIS LATER, ONCE PYSYNPHOT IS INCLUDED WITH PIPELINE DIST?
                if cval == 0:
                    print("Countrate value for {} is zero in {}.".format(self.params['Readout'][usefilt],self.parameters['phot_file']))
                    print("Eventually attempting to calculate value using pysynphot.")
                    print("but pysynphot is not present in jwst build 6, so pushing off to later...")
                    sys.exit()
                    cval = self.findCountrate(self.params['Readout'][usefilt])

                #translate to counts in single frame at requested array size
                framecounts = scale*cval*self.frametime
                rate = scale*cval

                #add the countrate and the counts per frame to pointSourceList
                #since they will be used in future calculations
                entry.append(rate)
                entry.append(framecounts)

                #add the good point source, including location and counts, to the pointSourceList
                filteredList.add_row(entry)

        #Write the results to a file
        filteredList.meta['comments'] = ["Field center (degrees): %13.8f %14.8f y axis rotation angle (degrees): %f  image size: %4.4d %4.4d\n" % (self.ra,self.dec,self.params['Telescope']['rotation'],nx,ny)]
        filteredOut = self.basename + '_galaxySources.list'
        filteredList.write(filteredOut,format='ascii',overwrite=True)
        return filteredList


    def create_galaxy(self,radius,ellipticity,sersic,posang,totalcounts):
        #given relevent parameters, create a model sersic image with a given radius, eccentricity, 
        #position angle, and total counts.

        #create the grid of pixels
        meshmax = np.min([np.int(self.ffsize*self.coord_adjust['y']),radius*100.])
        x,y = np.meshgrid(np.arange(meshmax), np.arange(meshmax))

        #center the galaxy in the array
        xc = meshmax/2
        yc = meshmax/2
        
        #create model
        #print('posang is {}'.format(posang))
        mod = Sersic2D(amplitude = 1,r_eff = radius, n=sersic, x_0=xc, y_0=yc, ellip=ellipticity, theta=posang)

        #create instance of model
        img = mod(x, y)

        #check to see if you've cropped too small and there is still significant signal
        #at the edges
        
        mxedge = np.max(np.array([np.max(img[:,-1]),np.max(img[:,0]),np.max(img[0,:]),np.max(img[-1,:])]))
        if mxedge > 0.001:
            print('Too small!')

        #scale such that the total number of counts in the galaxy matches the input
        summedcounts = np.sum(img)
        factor = totalcounts / summedcounts
        img = img * factor
        return img
    
    

    def makeGalaxyImage(self,file,psf):
        #Using the entries in the 'simSignals' 'galaxyList' file, create a countrate image
        #of model galaxies (sersic profile)

        #Read in the list of galaxies (positions and magnitides)
        glist, pixflag, radflag = self.readGalaxyFile(file)
        if pixflag:
            print("Galaxy list input positions assumed to be in units of pixels.")
        else:
            print("Galaxy list input positions assumed to be in units of RA and Dec.") 

        if radflag:
            print("Galaxy list input radii assumed to be in units of pixels.")
        else:
            print("Galaxy list input radii assumed to be in units of arcsec.") 


        #Extract and save only the entries which will land (fully or partially) on the
        #aperture of the output
        galaxylist = self.filterGalaxyList(glist,pixflag,radflag)

        #galaxylist is a table with columns:
        #'pixelx','pixely','RA','Dec','RA_degrees','Dec_degrees','radius','ellipticity','pos_angle','sersic_index','magnitude','countrate_e/s','counts_per_frame_e'

        #final output image
        origyd,origxd = self.nominal_dims
        #origyd,origxd = self.dark.data[0,0,:,:].shape
        yd = origyd
        xd = origxd
        
        #expand if a grism source image is being made
        xfact = 1
        yfact = 1
        if self.params['Output']['grism_source_image']:
            #xfact = self.grism_direct_factor
            #yfact = self.grism_direct_factor
            #elif 
            yd = np.int(origyd * self.coord_adjust['y'])
            xd = np.int(origxd * self.coord_adjust['x'])

        #create the final galaxy countrate image
        galimage = np.zeros((yd,xd))
        dims = galimage.shape

        # Create corresponding segmentation map
        segmentation = segmap.SegMap()
        segmentation.xdim = xd
        segmentation.ydim = yd
        segmentation.initialize_map()

        #Adjust the coordinate system of the galaxy list if working with a grism direct image output
        deltax = 0
        deltay = 0
        if self.params['Output']['grism_source_image']:
            deltax = np.int((dims[1] - origxd) / 2)
            deltay = np.int((dims[0] - origyd) / 2)

        #create attitude matrix so we can calculate the North->V3 angle for
        #each galaxy
        attitude_matrix = self.getAttitudeMatrix()
            
        #For each entry, create an image, and place it onto the final output image
        for entry in galaxylist:

            #Get position angle in the correct units. Inputs for each
            #source are degrees east of north. So we need to find the
            #angle between north and V3, and then the angle between
            #V3 and the y-axis on the detector. The former can be found
            #using rotations.posang(attitude_matrix,v2,v3). The latter
            #is just V3SciYAngle in the SIAF (I think???)
            #v3SciYAng is measured in degrees, from V3 towards the Y axis,
            #measured from V3 towards V2.
            north_to_east_V3ang = rotations.posangle(attitude_matrix,entry['V2'],entry['V3'])
            #xposang = (0-self.v3scixang) - (north_to_east_V3ang - entry['pos_angle'])            
            xposang = 0. - (self.v3scixang - north_to_east_V3ang + self.local_roll - entry['pos_angle'] + 90. + self.params['Telescope']['rotation'])
            
            #first create the galaxy image
            stamp = self.create_galaxy(entry['radius'],entry['ellipticity'],entry['sersic_index'],xposang*np.pi/180.,entry['counts_per_frame_e'])

            #convolve the galaxy with the NIRCam PSF
            stamp = s1.fftconvolve(stamp,psf,mode='same')
            
            #Now add the stamp to the main image
            #Extract the appropriate subarray from the galaxy image if necessary
            galdims = stamp.shape

            #print('requested radius: {}  stamp size: {}'.format(entry['radius'],galdims))

            nyshift = galdims[0] / 2
            nxshift = galdims[1] / 2

            nx = int(entry['pixelx']+deltax)
            ny = int(entry['pixely']+deltay)
            i1 = max(nx-nxshift,0)
            i2 = min(nx+1+nxshift,dims[1])
            j1 = max(ny-nyshift,0)
            j2 = min(ny+1+nyshift,dims[0])
            k1 = nxshift-(nx-i1)
            k2 = nxshift+(i2-nx)
            l1 = nyshift-(ny-j1)
            l2 = nyshift+(j2-ny)

            #if the cutout for the psf is larger than
            #the psf array, truncate it, along with the array
            #in the source image where it will be placed
            if l2 > galdims[0]:
                l2 = galdims[0]
                j2 = j1 + (l2-l1)

            if k2 > galdims[1]:
                k2 = galdims[1]
                i2 = i1 + (k2-k1)

            #At this point coordinates are in the final output array coordinate system, so there
            #should be no negative values, nor values larger than the output array size
            if j1 < 0 or i1<0 or l1<0 or k1<0:
                print(j1,i1,l1,k1)
                print('bad low')
                #sys.exit()
            if j2>(dims[0]+1) or i2>(dims[1]+1) or l2>(galdims[1]+1) or k2>(galdims[1]+1):
                print(j2,i2,l2,k2)
                print('bad high')
                #sys.exit()

            if ((j2 > j1) and (i2 > i1) and (l2 > l1) and (k2 > k1) and (j1 < dims[0]) and (i1 < dims[0])):
                galimage[j1:j2,i1:i2] = galimage[j1:j2,i1:i2] + stamp[l1:l2,k1:k2]
                # Divide readnoise by 100 sec, which is a 10 group RAPID ramp?
                noiseval = self.single_ron / 100. + self.params['simSignals']['bkgdrate']
                if self.params['Inst']['mode'].lower() == 'wfss':
                    noiseval += self.grism_background
                segmentation.add_object_noise(stamp[l1:l2,k1:k2],j1,i1,entry['index'],noiseval)

            else:
                print("Source located entirely outside the field of view. Skipping.")

        return galimage,segmentation.segmap


    def getExtendedSourceList(self,filename):
        #read in the list of point sources to add, and adjust the
        #provided positions for astrometric distortion

        extSourceList = Table(names=('pixelx','pixely','RA','Dec','RA_degrees','Dec_degrees','magnitude','countrate_e/s','counts_per_frame_e'),dtype=('f','f','S14','S14','f','f','f','f','f'))

        try:
            lines,pixelflag = self.readPointSourceFile(filename)
            if pixelflag:
                print("Extended source list input positions assumed to be in units of pixels.")
            else:
                print("Extended list input positions assumed to be in units of RA and Dec.") 
        except:
            print("WARNING: Unable to open the extended source list file {}".format(filename))
            sys.exit()

        #File to save adjusted point source locations
        eslist = open(self.params['Output']['file'][0:-5] + '_extendedsources.list','w')

        dtor = math.radians(1.)
        nx = (self.subarray_bounds[2]-self.subarray_bounds[0])+1
        ny = (self.subarray_bounds[3]-self.subarray_bounds[1])+1
        xc = (self.subarray_bounds[2]+self.subarray_bounds[0])/2.
        yc = (self.subarray_bounds[3]+self.subarray_bounds[1])/2.

        #Location of the subarray's reference pixel. 
        xrefpix = self.refpix_pos['x']
        yrefpix = self.refpix_pos['y']

        # center positions, sub-array sizes in pixels
        # now offset the field center to array center for astrometric distortion corrections
        coord_transform = None
        if self.runStep['astrometric']:

            #Read in the CRDS-format distortion reference file
            with AsdfFile.open(self.params['Reffiles']['astrometric']) as dist_file:
                coord_transform = dist_file.tree['model']

        #Using the requested RA,Dec of the reference pixel, along with the 
        #V2,V3 of the reference pixel, and the requested roll angle of the telescope
        #create a matrix that can be used to translate between V2,V3 and RA,Dec
        #for any pixel.
        #v2,v3 need to be in arcsec, and RA, Dec, and roll all need to be in degrees
        #attitude_matrix = rotations.attitude(self.refpix_pos['v2'],self.refpix_pos['v3'],self.ra,self.dec,self.params['Telescope']["rotation"])
        attitude_matrix = self.getAttitudeMatrix()
        
        #Write out the RA and Dec of the field center to the output file
        #Also write out column headers to prepare for source list
        eslist.write("# Field center (degrees): %13.8f %14.8f y axis rotation angle (degrees): %f  image size: %4.4d %4.4d\n" % (self.ra,self.dec,self.params['Telescope']['rotation'],nx,ny))
        eslist.write('#\n')
        eslist.write("#    RA_(hh:mm:ss)   DEC_(dd:mm:ss)   RA_degrees      DEC_degrees     pixel_x   pixel_y    magnitude   counts/sec    counts/frame\n")

        #Loop over input lines in the source list 
        all_stamps = []
        for values in lines:
            #try:
            #line below (if 1>0) used to keep the block of code below at correct indent for the try: above
            #the try: is commented out for code testing.
            if 1>0:
                try:
                    entry0 = float(values['x_or_RA'])
                    entry1 = float(values['y_or_Dec'])
                    print(entry0,entry1,pixelflag)
                    if not pixelflag:
                        ra_str,dec_str = self.makePos(entry0,entry1)
                        ra = entry0
                        dec = entry1
                except:
                    #if inputs can't be converted to floats, then 
                    #assume we have RA/Dec strings. Convert to floats.
                    ra_str = values['x_or_RA']
                    dec_str = values['y_or_Dec']
                    ra,dec = self.parseRADec(ra_str,dec_str)

                #Case where point source list entries are given with RA and Dec
                if not pixelflag:

                    #If distortion is to be included - either with or without the full set of coordinate
                    #translation coefficients
                    if self.runStep['astrometric']:
                        pixelx,pixely = self.RADecToXY_astrometric(ra,dec,attitude_matrix,coord_transform)
                    else:
                        #No distortion at all - "manual mode"
                        pixelx,pixely = self.RADecToXY_manual(ra,dec)

                else:
                    #Case where the point source list entry locations are given in units of pixels
                    #In this case we have the source position, and RA/Dec are calculated only so 
                    #they can be written out into the output source list file.

                    #Assume that the input x and y values are coordinate values
                    #WITHIN THE SPECIFIED SUBARRAY. So for example, a source in the file
                    #at 0,0 when you are making a SUB160 ramp will fall on the lower left
                    #corner of the SUB160 subarray, NOT the lower left corner of the full
                    #frame. 

                    pixelx = entry0
                    pixely = entry1

                    ra,dec,ra_str,dec_str = self.XYToRADec(pixelx,pixely,attitude_matrix,coord_transform)
                            

                #Get the input magnitude of the point source
                try:
                    mag=float(values['magnitude'])
                except:
                    mag = None

                #Now find out how large the extended source image is, so we
                #know if all, part, or none of it will fall in the field of view
                ext_stamp = fits.getdata(values['filename'])
                if len(ext_stamp.shape) != 2:
                    ext_stamp = fits.getdata(values['filename'],1)

                eshape = np.array(ext_stamp.shape)
                if len(eshape) == 2:
                    edgey,edgex = eshape / 2
                else:
                    print("WARNING, extended source image {} is not 2D! Not sure how to proceed. Quitting.".format(values['filename']))
                    sys.exit()
      

                #Define the min and max source locations (in pixels) that fall onto the subarray
                #Inlude the effects of a requested grism_direct image, and also keep sources that
                #will only partially fall on the subarray
                #pixel coords here can still be negative and kept if the grism image is being made
              
                #First, coord limits for just the subarray
                miny = 0
                maxy = self.subarray_bounds[3] - self.subarray_bounds[1] 
                minx = 0
                maxx = self.subarray_bounds[2] - self.subarray_bounds[0] 
        
                #Expand the limits if a grism direct image is being made
                if self.params['Output']['grism_source_image'] == True:
                    extrapixy = np.int((maxy+1)/2 * (self.coord_adjust['y'] - 1.))
                    miny -= extrapixy
                    maxy += extrapixy
                    extrapixx = np.int((maxx+1)/2 * (self.coord_adjust['x'] - 1.))
                    minx -= extrapixx
                    maxx += extrapixx

                #Now, expand the dimensions again to include point sources that fall only partially on the 
                #subarray
                miny -= edgey
                maxy += edgey
                minx -= edgex
                maxx += edgex

                #Keep only sources within the appropriate bounds
                if pixely > miny and pixely < maxy and pixelx > minx and pixelx < maxx:
                            
                    #set up an entry for the output table
                    entry = [pixelx,pixely,ra_str,dec_str,ra,dec,mag]

                    #save the stamp image after normalizing to a total signal of 1.
                    # and convolving with PSF if requested

                    if self.params['simSignals']['PSFConvolveExtended']:
                        ext_stamp = s1.fftconvolve(ext_stamp,self.centerpsf,mode='same')
                    
                    norm_factor = np.sum(ext_stamp)
                    ext_stamp /= norm_factor
                    all_stamps.append(ext_stamp)

                    #If a magnitude is given then adjust the countrate to match it
                    if mag is not None:
                        #translate magnitudes to countrate
                        scale = 10.**(0.4*(15.0-mag))

                        #get the countrate that corresponds to a 15th magnitude star for this filter
                        if self.params['Readout']['pupil'][0].upper() == 'F':
                            usefilt = 'pupil'
                        else:
                            usefilt = 'filter'
                        cval = self.countvalues[self.params['Readout'][usefilt]]

                        #DEAL WITH THIS LATER, ONCE PYSYNPHOT IS INCLUDED WITH PIPELINE DIST?
                        if cval == 0:
                            print("Countrate value for {} is zero in {}.".format(self.params['Readout'][usefilt],self.parameters['phot_file']))
                            print("Eventually attempting to calculate value using pysynphot.")
                            print("but pysynphot is not present in jwst build 6, so pushing off to later...")
                            sys.exit()
                            cval = self.findCountrate(self.params['Readout'][usefilt])

                        #translate to counts in single frame at requested array size
                        framecounts = scale*cval*self.frametime
                        countrate = scale*cval
                        magwrite = mag

                    else:
                        #In this case, no magnitude is given in the extended input list
                        #Assume the input stamp image is in units of e/sec then.
                        print("No magnitude given for extended source in {}.".format(values['filename']))
                        print("Assuming the original file is in units of counts per sec.")
                        print("Multiplying original file values by 'extendedscale'.")
                        countrate = norm_factor * self.params['simSignals']['extendedscale']
                        framecounts = countrate*self.frametime
                        magwrite = 99.99999

                    #add the countrate and the counts per frame to pointSourceList
                    #since they will be used in future calculations
                    #entry.append(scale)
                    entry.append(countrate)
                    entry.append(framecounts)

                    #add the good point source, including location and counts, to the pointSourceList
                    #self.pointSourceList.append(entry)
                    extSourceList.add_row(entry)

                    #write out positions, distances, and counts to the output file
                    eslist.write("%s %s %14.8f %14.8f %9.3f %9.3f  %9.3f  %13.6e   %13.6e\n" % (ra_str,dec_str,ra,dec,pixelx,pixely,magwrite,countrate,framecounts))
                #except:
                #    print("ERROR: bad point source line %s. Skipping." % (line))
        print("Number of extended sources found within the requested aperture: {}".format(len(extSourceList)))
        #close the output file
        eslist.close()
        
        #If no good point sources were found in the requested array, alert the user
        if len(extSourceList) < 1:
            print("Warning: no non-sidereal extended sources within the requested array.")
            print("The extended source image option is being turned off")
        
        return extSourceList,all_stamps

    
    def makeExtendedSourceImage(self,extSources,extStamps):
        dims = np.array(self.nominal_dims)
        #dims = np.array(self.dark.data[0,0,:,:].shape)

        #offset that needs to be applied to the x,y positions of the
        #source list to account for case where we make a point
        #source image that is extra-large, to be used as a grism 
        #direct image
        deltax = 0
        deltay = 0

        newdimsx = np.int(dims[1] * self.coord_adjust['x'])
        newdimsy = np.int(dims[0] * self.coord_adjust['y'])
        deltax = self.coord_adjust['xoffset']
        deltay = self.coord_adjust['yoffset']
        dims = np.array([newdimsy,newdimsx])

        #create the empty image
        extimage = np.zeros((dims[0],dims[1]))

        # Create corresponding segmentation map
        segmentation = segmap.SegMap()
        segmentation.xdim = newdimsx
        segmentation.ydim = newdimsy
        segmentation.initialize_map()

        #Loop over the entries in the point source list
        for entry,stamp in zip(extSources,extStamps):
            #adjust x,y position if the grism output image is requested
            xpos = entry['pixelx'] + deltax
            ypos = entry['pixely'] + deltay
            xoff = math.floor(xpos)
            yoff = math.floor(ypos)

            #desired counts per second in the source
            counts = entry['countrate_e/s'] #/ self.frametime

            #Extract the appropriate subarray from the PSF image if necessary
            #Assume that the brightest pixel corresponds to the peak of the source
            psfdims = stamp.shape
            nyshift,nxshift = np.array(psfdims) / 2
            nx = int(xoff)
            ny = int(yoff)
            i1 = max(nx-nxshift,0)
            i2 = min(nx+1+nxshift,dims[1])
            j1 = max(ny-nyshift,0)
            j2 = min(ny+1+nyshift,dims[0])
            k1 = nxshift-(nx-i1)
            k2 = nxshift+(i2-nx)
            l1 = nyshift-(ny-j1)
            l2 = nyshift+(j2-ny)

            #if the cutout for the psf is larger than
            #the psf array, truncate it, along with the array
            #in the source image where it will be placed
            if l2 > psfdims[0]:
                l2 = psfdims[0]
                j2 = j1 + (l2-l1)

            if k2 > psfdims[1]:
                k2 = psfdims[1]
                i2 = i1 + (k2-k1)

            #At this point coordinates are in the final output array coordinate system, so there
            #should be no negative values, nor values larger than the output array size
            if j1 < 0 or i1<0 or l1<0 or k1<0:
                print(j1,i1,l1,k1)
                print('bad low')
                sys.exit()
            if j2>(dims[0]+1) or i2>(dims[1]+1) or l2>(psfdims[1]+1) or k2>(psfdims[1]+1):
                print(j2,i2,l2,k2)
                print('bad high')
                sys.exit()

            #Add stamp image to the extended source countrate image
            extimage[j1:j2,i1:i2] = extimage[j1:j2,i1:i2] + stamp[l1:l2,k1:k2]*counts
            # Divide readnoise by 100 sec, which is a 10 group RAPID ramp?
            noiseval = self.single_ron / 100. + self.params['simSignals']['bkgdrate']
            if self.params['Inst']['mode'].lower() == 'wfss':
                noiseval += self.grism_background
            segmentation.add_object_noise(stamp[l1:l2,k1:k2]*counts,j1,i1,entry['index'],noiseval)

            
        return extimage, segmentation.segmap


    def makeFilterTable(self):
        #Create the table that contains the possible filter list, quantum yields, and countrates for a
        #star with vega magnitude of 15 in each filter. Do this by reading in phot_file
        #listed in the parameter file. 

        #FUTURE WORK: If the countrates are left as 0, then pysynphot will
        #be used to calculate them later
        try:
            cvals_tab = ascii.read(self.params['Reffiles']['phot'])
            instrumentfilternames = cvals_tab['filter'].data
            stringcountrates = cvals_tab['countrate_for_vegamag15'].data
            instrumentmag15countrates = [float(s) for s in stringcountrates]
            strinstrumentqy = cvals_tab['quantum_yield'].data
            qy = [float(s) for s in strinstrumentqy]
            self.countvalues=dict(zip(instrumentfilternames,instrumentmag15countrates))
            self.qydict=dict(zip(instrumentfilternames,qy))

        except:
            print("WARNING: Unable to read in {}.".format(self.params['Reffiles']['phot']))
            sys.exit()


    def readParameterFile(self):
        #read in the parameter file
        try:
            with open(self.paramfile,'r') as infile:
                self.params = yaml.load(infile)
        except:
            print("WARNING: unable to open {}".format(self.paramfile))
            sys.exit()

            
    def checkParams(self):
        #check instrument name
        if self.params['Inst']['instrument'].lower() not in inst_list:
            print("WARNING: {} instrument not implemented within ramp simulator")
            sys.exit()

        #check entred mode: 
        possibleModes = modes[self.params['Inst']['instrument'].lower()]
        self.params['Inst']['mode'] = self.params['Inst']['mode'].lower()
        if self.params['Inst']['mode'] in possibleModes:
            pass
        else:
            print("WARNING: unrecognized mode {} for {}. Must be one of: {}".format(self.params['Inst']['mode'],self.params['Inst']['instrument'],possibleModes))
            sys.exit()

        #check the number of amps. Full frame data will always be collected using 4 amps. 
        #Subarray data will always be 1 amp, except for the grism subarrays which span the
        #entire width of the detector. Those can be read out using 1 or 4 amps.
        

        #Make sure that the requested number of groups is less than or equal to the maximum
        #allowed. If you're continuing on with an unknown readout pattern (not recommended)
        #then assume a max of 10 groups.
        #For science operations, ngroup is going to be limited to 10 for all readout patterns
        #except for the DEEP patterns, which can go to 20.
        #match = self.readpatterns['name'] == self.params['Readout']['readpatt'].upper()
        #if sum(match) == 1:
        #    maxgroups = self.readpatterns['maxgroups'].data[match][0]
        #if sum(match) == 0:
        #    print("Unrecognized readout pattern {}. Assuming a maximum allowed number of groups of 10.".format(self.params['Readout']['readpatt']))
        #    maxgroups = 10

        #if (self.params['Readout']['ngroup'] > maxgroups):
        #    print("WARNING: {} is limited to a maximum of {} groups. Proceeding with ngroup = {}.".format(self.params['Readout']['readpatt'],maxgroups,maxgroups))
        #    self.params['Readout']['readpatt'] = maxgroups

        #check for entries in the parameter file that are None or blank,
        #indicating the step should be skipped. Create a dictionary of steps 
        #and populate with True or False
        self.runStep = {}
        #self.runStep['superbias'] = self.checkRunStep(self.params['Reffiles']['superbias'])
        #self.runStep['nonlin'] = self.checkRunStep(self.params['Reffiles']['linearity'])
        #self.runStep['gain'] = self.checkRunStep(self.params['Reffiles']['gain'])
        self.runStep['phot'] = self.checkRunStep(self.params['Reffiles']['phot'])
        self.runStep['pixelflat'] = self.checkRunStep(self.params['Reffiles']['pixelflat'])
        self.runStep['illuminationflat'] = self.checkRunStep(self.params['Reffiles']['illumflat'])
        self.runStep['astrometric'] = self.checkRunStep(self.params['Reffiles']['astrometric'])
        self.runStep['distortion_coeffs'] = self.checkRunStep(self.params['Reffiles']['distortion_coeffs'])
        self.runStep['ipc'] = self.checkRunStep(self.params['Reffiles']['ipc'])
        self.runStep['crosstalk'] = self.checkRunStep(self.params['Reffiles']['crosstalk'])
        self.runStep['occult'] = self.checkRunStep(self.params['Reffiles']['occult'])
        self.runStep['pointsource'] = self.checkRunStep(self.params['simSignals']['pointsource'])
        self.runStep['galaxies'] = self.checkRunStep(self.params['simSignals']['galaxyListFile'])
        self.runStep['extendedsource'] = self.checkRunStep(self.params['simSignals']['extended'])
        self.runStep['movingTargets'] = self.checkRunStep(self.params['simSignals']['movingTargetList'])
        self.runStep['movingTargetsSersic'] = self.checkRunStep(self.params['simSignals']['movingTargetSersic'])
        self.runStep['movingTargetsExtended'] = self.checkRunStep(self.params['simSignals']['movingTargetExtended'])
        self.runStep['MT_tracking'] = self.checkRunStep(self.params['simSignals']['movingTargetToTrack'])
        self.runStep['zodiacal'] = self.checkRunStep(self.params['simSignals']['zodiacal'])
        self.runStep['scattered'] = self.checkRunStep(self.params['simSignals']['scattered'])
        self.runStep['fwpw'] = self.checkRunStep(self.params['Reffiles']['filtpupilcombo'])
        self.runStep['pixelAreaMap'] = self.checkRunStep(self.params['Reffiles']['pixelAreaMap'])

        #create table that will contain filters/quantum yield/and vegamag=15 countrates
        self.makeFilterTable()

        #make sure the requested filter is allowed. For imaging, all filters are allowed.
        #In the future, other modes will be more restrictive
        if self.params['Readout']['pupil'][0].upper() == 'F':
            usefilt = 'pupil'
        else:
            usefilt = 'filter'
        if self.params['Readout'][usefilt] not in self.qydict:
            print("WARNING: requested filter {} is not in the list of possible filters.".format(self.params['Readout'][usefilt]))
            sys.exit()

        #PSF: generate the name of the PSF file to use
        #if the psf path has been left blank or set to 'None'
        #then assume the user does not want to add point sources
        if self.params['simSignals']['psfpath'] is not None:
            if self.params['simSignals']['psfpath'][-1] != '/':
                self.params['simSignals']['psfpath']=self.params['simSignals']['psfpath']+'/'

            wfe = self.params['simSignals']['psfwfe']
            wfegroup = self.params['simSignals']['psfwfegroup']
            basename = self.params['simSignals']['psfbasename'] + '_'
            if wfe == 0:
                psfname=basename+self.params['simSignals'][usefilt].lower()+'_zero'
                self.params['simSignals']['psfpath']=self.params['simSignals']['psfpath']+self.params['simSignals'][usefilt].lower()+'/zero/'
            else:
                if wfe in [115,123,136,155] and wfegroup > -1 and wfegroup < 10:
                    psfname=basename+self.params['Readout'][usefilt].lower()+"_"+str(wfe)+"_"+str(wfegroup)
                    self.params['simSignals']['psfpath']=self.params['simSignals']['psfpath']+self.params['Readout'][usefilt].lower()+'/'+str(wfe)+'/'
                else:
                    print("WARNING: wfe value of {} not allowed. Library does not exist.".format(wfe))
                    sys.exit()
                self.psfname = self.params['simSignals']['psfpath'] + psfname

        else:
            #case where psfPath is None. In this case, create a PSF on the fly to use
            #for adding sources
            print("update this to include a call to WebbPSF????????")
            self.psfimage=np.zeros((5,5),dtype=np.float32)
            sum1=0
            for i in range(5):
                for j in range(5):
                    self.psfimage[i,j]=(0.02**(abs((i-2)))*(0.02**abs(j-2)))
                    sum1=sum1+self.psfimage[i,j]
            self.psfimage=self.psfimage/sum1
            self.psfname = None
                
        # Read in the 'central' PSF file. This is the file that
        # has the PSF centered on the pixel. This will be used
        # if there are sersic or extended sources that need to
        # be convolved with the NIRCam PSF before adding
        centerpsffile = os.path.join(self.params['simSignals']['psfpath'],psfname+'_0p0_0p0.fits')
        self.centerpsf = fits.getdata(centerpsffile)
        self.centerpsf = self.cropPSF(self.centerpsf)

        #normalize the PSF to a total signal of 1.0
        totalsignal = np.sum(self.centerpsf)
        self.centerpsf /= totalsignal
   
        #ASTROMETRY
        #Read in the distortion coefficients file if present. These will provide a more exact
        #transform from RA,Dec to x,y than the astrometric distortion reference file above.
        #The file above can be off by ~20 pixels in the corners of the array. This file will give
        #exact answers
        if self.runStep['distortion_coeffs'] == True:
            if os.path.isfile(self.params['Reffiles']['distortion_coeffs']):
                distortionTable = ascii.read(self.params['Reffiles']['distortion_coeffs'],header_start=1)
            else:
                print("WARNING: Input distortion coefficients file {} does not exist.".format(self.params['Reffile']['distortion_coeffs']))
                sys.exit()

            #read in coefficients for the forward 'science' to 'ideal' coordinate transformation.
            #'science' is in units of distorted pixels, while 'ideal' is the undistorted
            #angular distance from the reference pixel
            ap_name = self.params['Readout']['array_name']

            self.x_sci2idl,self.y_sci2idl,self.v2_ref,self.v3_ref,self.parity,self.v3yang,self.xsciscale,self.ysciscale,self.v3scixang = self.getDistortionCoefficients(distortionTable,'science','ideal',ap_name)
            
            #Generate the coordinate transform for V2,V3 to 'ideal'
            siaf = ascii.read(self.params['Reffiles']['distortion_coeffs'],header_start=1)
            match = siaf['AperName'] == ap_name
            if np.any(match) == False:
                print("Aperture name {} not found in input CSV file.".
                      format(aperture))
                sys.exit()

            siaf_row = siaf[match]


            #self.v2v32idlx, self.v2v32idly = read_siaf_table.get_siaf_v2v3_transform(self.params['Reffiles']['distortion_coeffs'],ap_name,to_system='ideal')
            self.v2v32idlx, self.v2v32idly = read_siaf_table.get_siaf_v2v3_transform(siaf_row,ap_name,to_system='ideal')

            
        #convert the input RA and Dec of the pointing position into floats
        #check to see if the inputs are in decimal units or hh:mm:ss strings
        try:
            self.ra = float(self.params['Telescope']['ra'])
            self.dec = float(self.params['Telescope']['dec'])
        except:
            self.ra,self.dec=self.parseRADec(self.params['Telescope']['ra'],self.params['Telescope']['dec'])

        if abs(self.dec) > 90. or self.ra < 0. or self.ra > 360. or self.ra is None or self.dec is None:
            print("WARNING: bad requested RA and Dec {} {}".format(self.ra,self.dec))
            sys.exit()

        #make sure the rotation angle is a float
        try:
            self.params['Telescope']["rotation"]=float(self.params['Telescope']["rotation"])
        except:
            print("ERROR: bad rotation value {}, setting to zero.".format(self.params['Telescope']["rotation"]))
            self.params['Telescope']["rotation"]=0.


        #check that the various scaling factors are floats and within a reasonable range
        #self.params['cosmicRay']['scale'] = self.checkParamVal(self.params['cosmicRay']['scale'],'cosmicRay',0,100,1)
        self.params['simSignals']['extendedscale'] = self.checkParamVal(self.params['simSignals']['extendedscale'],'extendedEmission',0,10000,1)
        self.params['simSignals']['zodiscale'] = self.checkParamVal(self.params['simSignals']['zodiscale'],'zodi',0,10000,1)
        self.params['simSignals']['scatteredscale'] = self.checkParamVal(self.params['simSignals']['scatteredscale'],'scatteredLight',0,10000,1)

        #make sure the requested output format is an allowed value
        if self.params['Output']['format'] not in allowedOutputFormats:
            print("WARNING: unsupported output format {} requested. Possible options are {}.".format(self.params['Output']['format'],allowedOutputFormats))
            sys.exit()

        #Entries for creating the grims input image
        if not isinstance(self.params['Output']['grism_source_image'],bool):
            if self.params['Output']['grism_source_image'].lower() == 'none':
                self.params['Output']['grism_source_image'] = False
            else:
                print("WARNING: grism_source_image needs to be True or False")
                sys.exit()

        #if not isinstance(self.params['Output']['grism_input_only'],bool):
        #    if self.params['Output']['grism_source_image'].lower() == 'none':
        #        self.params['Output']['grism_input_only'] = False
        #    else:
        #        print("WARNING: grism_input_only needs to be True or False")
        #        sys.exit()

        #Location of extended image on output array, pixel x,y values.
        try:
            self.params['simSignals']['extendedCenter'] = np.fromstring(self.params['simSignals']['extendedCenter'], dtype=int, sep=",")
        except:
            print("WARNING: not able to parse the extendedCenter list {}. It should be a comma-separated list of x and y pixel positions.".format(self.params['simSignals']['extendedCenter']))
            sys.exit()

        #Time series settings
        #if self.params['Inst']['mode'] == 'tso':
        #    
        #    #make sure slew rate and angle are floats
        #    try:
        #        self.params['Telescope']['slewRate'] = np.float(self.params['Telescope']['slewRate'])
        #    except:
        #        print("WARNING: input slew rate {} is not an integer or float.".format(self.params['Telescope']['slewRate']))
        #        sys.exit()
        #
        #    try:
        #        self.params['Telescope']['slewAngle'] = np.float(self.params['Telescope']['slewAngle'])
        #    except:
        #        print("WARNING: input slew angle {} is not an integer or float.".format(self.params['Telescope']['slewAngle']))
        #        sys.exit()


        #check the output metadata, including visit and observation numbers, obs_id, etc
        #
        #kwchecks = ['program_number','visit_number','visit_group',
        #            'sequence_id','activity_id','exposure_number','observation_number','obs_id','visit_id']
        #for quality in kwchecks:
        #    try:
        #        self.params['Output'][quality] = str(self.params['Output'][quality])
        #    except:
        #        print("WARNING: unable to convert {} to string. This is required.".format(self.params['Output'][quality]))
        #        sys.exit()

    def checkRunStep(self,filename):
        #check to see if a filename exists in the parameter file.
        if ((len(filename) == 0) or (filename.lower() == 'none')):
            return False
        else:
            return True
        

    def checkParamVal(self,value,typ,vmin,vmax,default):
        #make sure the input value is a float and between min and max
        try:
            value = float(value)
        except:
            print("WARNING: {} for {} is not a float.".format(value,typ))
            sys.exit()
            
        if ((value >= vmin) & (value <= vmax)):
            return value
        else:
            print("ERROR: {} for {} is not within reasonable bounds. Setting to {}".format(value,typ,default))
            return default


    def readSubarrayDefinitionFile(self):
        #read in the file that contains a list of subarray names and positions on the detector
        try:
            self.subdict = ascii.read(self.params['Reffiles']['subarray_defs'],data_start=1,header_start=0)
        except:
            print("Error: could not read in subarray definitions file.")
            sys.exit()

        
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


    def instrument_specific_dicts(self,instrument):
        #get instrument-specific values for things that
        #don't need to be in the parameter file

        #array size of a full frame image
        self.ffsize = full_array_size[instrument]

        #pixel scale - return as a 2-element list, with pixscale for x and y.
        if instrument.lower() == 'nircam':
            filt = self.params['Readout']['filter']
            fnum = int(filt[1:4])
            if fnum < 230:
                channel = 'sw'
            else:
                channel = 'lw'
            self.pixscale = [pixelScale[instrument][channel],pixelScale[instrument][channel]]
        else:
            self.pixscale = [pixelScale[instrument],pixelScale[instrument]]

            
    def saveSingleFits(self,image,name,key_dict=None,image2=None,image2type=None):
        #save an array into the first extension of a fits file
        h0 = fits.PrimaryHDU()
        h1 = fits.ImageHDU(image)
        if image2 is not None:
            h2 = fits.ImageHDU(image2)
            if image2type is not None:
                h2.header['EXTNAME'] = image2type

        #if a keyword dictionary is provided, put the 
        #keywords into the 0th and 1st extension headers
        if key_dict is not None:
            for key in key_dict:
                h0.header[key] = key_dict[key]
                h1.header[key] = key_dict[key]

        if image2 is None:
            hdulist = fits.HDUList([h0,h1])
        else:
            hdulist = fits.HDUList([h0,h1,h2])
        hdulist.writeto(name,overwrite=True)


    def add_options(self,parser=None,usage=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,description='Create seed image via catalogs')
        parser.add_argument("paramfile",help='File describing the input parameters and instrument settings to use. (YAML format).')
        parser.add_argument("--param_example",help='If used, an example parameter file is output.')
        return parser


if __name__ == '__main__':

    usagestring = 'USAGE: catalog_seed_image.py inputs.yaml'

    seed = Catalog_seed()
    parser = seed.add_options(usage = usagestring)
    args = parser.parse_args(namespace=seed)
    seed.run()
