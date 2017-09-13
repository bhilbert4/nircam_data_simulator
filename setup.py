#! /usr/bin/env python

from distutils.core import setup

setup(name='nircam_data_simulator',
      version = '1.0',
      description='JWST Multiaccum Ramp Simulator',
      long_description='A tool to create simulated raw NIRCAm ramps based on dark current ramps from ground testing in combination with source catalogs.',
      url='https://github.com/bhilbert4/nircam_data_simulator',
      author='Bryan Hilbert',
      author_email='hilbert@stsci.edu',
      py_modules = ['catalog_seed_image','polynomial','rotations',
                    'dark_prep','read_fits','moving_targets',
                    'read_siaf_table','segmentation_map',
                    'set_telescope_pointing_separated'],
    )
