import sys
import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import pylab                   # Needed to plot histograms.
from darepype.drp.datafits import DataFits     # Gets the function that makes DataFits io objects.
from astropy.io import fits                    # Need this if you want to use astropy.io io objects.
from ipywidgets import interact                # Need this for interactive plots.
from matplotlib.colors import LogNorm          # Machinery for LogNorm scaling of intensities.
from matplotlib.colors import SymLogNorm       # Machinery for SymLogNorm scaling of intensities.
from matplotlib.colors import PowerNorm        # Machinery for LogNorm (e.g., square root) scaling of intensities.
from astropy.stats import mad_std              # The median absolute deviation, a more robust estimator than std.
import scipy.ndimage as nd                        # Various algorithms for image transformations.
from astropy.time import Time
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
from astropy import wcs
from astroscrappy import detect_cosmics
from drizzle import drizzle as drz

# Specify data paths
newkeys = True    # True if the input files have new keywords introduced on 2021-09-02.
reduced_by = 'imagemath_jonah'   # Version of notebook used in this reduction.
datapath = ''
outpath  = ''
darkpath = ''
biaspath = ''
flatpath = ''
# Specify camera
camera = 'SBIG'
# Set config file (avoids error messages that occur if config isn't
# specified in Datafits() statement).
config = os.getcwd() + '/content/pipeline/config/pipeconf_SEO.txt'