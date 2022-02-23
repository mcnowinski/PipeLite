from setup import *
# Load data
datain = []
newname = ""
files = []

"""
load_data (internal)
Loads data from the outfolder with the fiven filter
"""
def load_data(inpath = outpath, filter = "WCS"):
    somefiles = [f for f in os.listdir(inpath) if '.fits' in f and filter in f]
    print(somefiles)
    global files
    files = sorted(somefiles)   ## Put files in alphabatical order, if necessary.
    for i in range(len(somefiles)):
        print( i, somefiles[i])
    global newname
    newname = f"{files[0][:-11]}_DRZ.fits"
    global datain
    datain = [DataFits(f"{inpath}/{f}") for f in files]
    print('Data loaded')

"""
display_images:
    Display all the images in the outpath folder.
    This is useful for checking that the data is loaded correctly.
"""
def display_images(folder = outpath, filter = "WCS"):
    load_data(folder, filter)
    show_all = True
    show_grid = False
    figx, figy = 8,8   # Use these to get a quick look and save space.
    meds, mads = np.zeros((len(datain))), np.zeros((len(datain)))
    airmass = np.zeros((len(datain)))
    for i in range(len(datain)):
        img = datain[i].image
        med, mad = np.nanmedian(img), mad_std(img,ignore_nan=True)
        vmin, vmax = med - mad * 3.0, med + mad * 20.0
        meds[i], mads[i] = med, mad
        airmass[i] = datain[i].header['airmass']
        # print(datain[i].header['RA'])
        # print(datain[i].header['TRAKHA'], datain[i].header['TRAKDEC'])
        if show_all == True:
            plt.figure(figsize = (figx,figy))
            if show_grid == True: plt.grid()
            plt.title(f"{files[i]} Median = {str(med)}  mad_std = {str(mad)}\n RA = {str(datain[i].header['RA'])}  DEC = {str(datain[i].header['DEC'])}\n  TRAKHA = {str(datain[i].header['TRAKHA'])} TRAKDEC = {str(datain[i].header['TRAKDEC'])} \n")
            plt.imshow(datain[i].image,'gray_r',interpolation = 'nearest',vmax=vmax,vmin=vmin)

"""
drizzle_helper:
    - Takes a single image and drizzles it.
    - Saves the drizzled image to the outpath.
"""
def drizzle_helper(filter, infolder, outfolder):
    load_data(inpath = infolder, filter=filter)
    kernel = 'square'
    pixfrac = 1.0
    resolution = 1.0
    pad = 0.0
    fillval = np.nan
    drizzleweights = 'exptime'
    outangle = 0.0
    try:
        det = np.linalg.det(wcs.WCS(datain[0].header).wcs.cd)
        pscale = np.sqrt(np.abs(det))*3600.
    except:
        try:
            det = np.linalg.det(wcs.WCS(datain[0].header).wcs.pc)
            pscale = np.sqrt(np.abs(det))*3600.
        except:
            pscale = self.datain[0].header['PIXSCAL']

    #calculations necessary for updating wcs information
    px = []
    py = []

    #in order to avoid NaN interactions, creating weight map
    weights=[]
    for f in datain:
        weights.append((np.where(np.isnan(f.image) == True, 0, 1)))

    for f in datain:
        px.extend(wcs.WCS(f.header).calc_footprint()[:,0])
        py.extend(wcs.WCS(f.header).calc_footprint()[:,1])
    x0 = (max(px)+min(px))/2.
    y0 = (max(py)+min(py))/2.
    sx = (max(px)-min(px))*np.cos(y0/180*np.pi) # arcsec
    sy = (max(py)-min(py)) # arcsec
    size = (sx*3600+(pad)*2, sy*3600+(pad)*2)
    xpix = size[0]//pscale
    ypix = size[1]//pscale
    cdelt = [pscale/3600.]*2
    dataout = DataFits()
    dataout.header = datain[0].header.copy()
    dataout.header['CRPIX1'] = xpix/2
    dataout.header['CRPIX2'] = ypix/2
    dataout.header['CRVAL1'] = x0
    dataout.header['CRVAL2'] = y0
    dataout.header['CD1_1'] = -cdelt[0]
    dataout.header['CD1_2'] = dataout.header['CD2_1'] = 0.
    dataout.header['CD2_2'] = cdelt[1]
    dataout.header['NAXIS1'] = int(xpix)
    dataout.header['NAXIS2'] = int(ypix)
    dataout.header['CTYPE1'] = 'RA---TAN-SIP'
    dataout.header['CTYPE2'] = 'DEC--TAN-SIP'
    dataout.header['RADESYS'] = 'ICRS'
    dataout.header['EQUINOX'] = 2000
    dataout.header['LATPOLE'] = datain[0].header['CRVAL2']
    dataout.header['LONPOLE'] = 180
    dataout.header['PIXASEC'] = pscale

    theta_rad = np.deg2rad(outangle)
    rot_matrix = np.array([[np.cos(theta_rad), -np.sin(theta_rad)], 
                    [np.sin(theta_rad),  np.cos(theta_rad)]])
    rot_cd = np.dot(rot_matrix, np.array([[dataout.header['CD1_1'], 0.],[0., dataout.header['CD2_2']]]))
    for i in [0,1]:
        for j in [0,1]:
            dataout.header['CD{0:d}_{1:d}'.format(i+1, j+1)] = rot_cd[i,j]

    ##check drizzle arguments
    if kernel == 'smoothing':
        kernel = 'lanczos3'
    elif kernel in ['square', 'point', 'gaussian', 'tophat']:
        kernel = kernel
    else:
        print('Kernel name not recognized, using default = square')
        kernel = 'square'
    if drizzleweights == 'uniform':
        driz_wt = ''
    elif drizzleweights in ['exptime', 'expsq']:
        driz_wt = drizzleweights
    else:
        print('Drizzle weighting not recognized, using default = null string')
        driz_wt = ''
    
    ##create drizzle object and add input images
    fullwcs = wcs.WCS(dataout.header)
    print('Starting drizzle')
    driz = drz.Drizzle(outwcs = fullwcs, pixfrac=pixfrac, \
                    kernel=kernel, fillval='10000', wt_scl=driz_wt)
    for i,f in enumerate(datain):
        print(i, 'Adding %s to drizzle stack' % f.filename)
        driz.add_image(f.imgdata[0], wcs.WCS(f.header), inwht=weights[i])

    try:
        fillval=float(fillval)
    except:
        fillval=np.nan
        print('Fillvalue not recognized or missing, using default = np.nan')

    dataout.imageset(np.where(driz.outsci == 10000, fillval, driz.outsci))
    dataout.imageset(driz.outwht,'OutWeight', dataout.header)
    dataout.filename = datain[0].filename

    #add history
    dataout.setheadval('HISTORY','Coadd: %d files combined with %s kernel, pixfrac %f at %f times resolution' \
                            % (len(datain), kernel, pixfrac, resolution))

    outf = os.path.join(outfolder, newname)
    dataout.filename = outf
    print(dataout.filename)
    outd = dataout
    outd.save(outf)

"""
drizzle:
    Drizzles all WCS images in the output folder with given filters
    inputs: 
        - filter (lst str): list of filter strings to be drizzled. Will combine all images containing each filter into a single image.
        - infolder (str): path to folder containing images to stack (default: outpath set in setup.py)
        Note: Infolder is set to outpath by default so it can take the output of the previous steps
        - outfolder (str): path to output folder (default: outpath set in setup.py)
    outputs:
        - Saves drizzled images to outpath
    example:
        drizzle(filter=['g-band','i-band'], outpath) drizzles all images containing the string "g-band into their respective files
"""
def drizzle(filter = [""], infolder = outpath, outfolder = outpath):
    for i in filter:
        drizzle_helper(i, infolder, outfolder)
