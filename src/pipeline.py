from setup import *
from darepype.drp import PipeLine
import subprocess

'''File and path definitions.'''
# Location of the pipeline source code.
codefolder = os.getcwd() + '/content'
logfilename = codefolder + '/pipelog.txt'
# Location the data to be processed.
datafolder = outpath
# Location of the baseconfig file.
baseconfig = os.path.join(codefolder,'pipeline/config','pipeconf_SEO.txt')
# Location of delta config file.
dconfig = os.path.join(os.getcwd(), 'dconf_dah.txt')
# Set the path to the pipeline source code.
sys.path.append(os.path.join(codefolder,'pipeline/source'))
print('Pipeline setup complete')

## For ALL the files in the directory.
"""
run_pipeline:
    Runs the pipeline with the hotpix and astrometry steps
    on all the files in the directory.

    Inputs:
    folder (str): The folder containing the files to be processed. (Default: outpath)

    Outputs:
    Saves HPX and WPS files in the outpath.
"""
def run_pipeline(folder = outpath):
    files = sorted([f for f in os.listdir(folder) if '.fits' in f and 'HDR' in f])
    print("Running pipeline on the following files:")
    for i in range(len(files)):
        print(files[i])

    for file in files:
        infile = os.path.join(folder,file)
        pipe = PipeLine(config =[baseconfig, dconfig]) # Make a pipe object
        result = pipe(infile, pipemode='postbdf', force=True) # Run the pipeline
        result.save()

def run_local_astrometry(folder = outpath):
    files = sorted([f for f in os.listdir(folder) if '.fits' in f and 'HDR' in f and "WCS" not in f])
    print("Running local astrometry on the following files:")
    for i in range(len(files)):
        print(files[i])
    nospace = folder.replace(' ','\ ')
    out = subprocess.Popen(f"cd {nospace} && solve-field --crpix-center --no-verify --scale-units arcsecperpix --scale-low 0.8 --scale-high 1.0 --downsample 4 --overwrite *.fits", shell=True, stdout=subprocess.PIPE)
    subprocess_return = out.stdout.read()
    
    # rename all files with .new extension to .fits extension
    for file in os.listdir(folder):
        if file.endswith('.new'):
            newfile = file.replace('HDR.new','WCS.fits')
            os.rename(os.path.join(folder,file), os.path.join(folder,newfile))
        elif not "fits" in file:
            # remove all files that are not .fits
            os.remove(os.path.join(folder,file))
    # subprocess run with output