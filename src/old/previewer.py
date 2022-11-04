#!/usr/bin/env python

import os
import glob
from astropy.visualization import ZScaleInterval
from astropy.io import fits
from matplotlib import pyplot as plt

# MODIFY THESE FIELDS AS NEEDED!
# input path *with* ending forward slash
input_path = './'
# output path *with* ending forward slash
output_path = './preview/'

# does output directory exist? if not, create it
try:
    os.mkdir(output_path)
except:
    pass

# get a list of all FITS files in the input directory
fits_files = glob.glob(input_path+'*.fits') + glob.glob(input_path+'*.fit')
# loop through all FITS files and generate PNG images
for fits_file in sorted(fits_files):
    data = fits.getdata(fits_file) # get FITS data
    zscale = ZScaleInterval() # scale using ZScale
    fig, ax = plt.subplots(figsize=plt.figaspect(data)) # remove borders
    fig.subplots_adjust(0,0,1,1)    
    plt.imshow(zscale(data), cmap="gray") # create image
    plt.savefig(output_path + os.path.basename(fits_file) + '.jpg') # save image
