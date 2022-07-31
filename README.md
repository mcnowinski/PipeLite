<h1>Welcome to PipeLite for ASTR 21200</h1>

The goal of this repo is to abstract as much behind-the-scenes functionality as possible to make the image processing pipeline simple and easy to run.

It combines image_math_BDF, the HotPix and Astrometry steps of the pipeline, and image_stack.

<h2> Running Locally </h2>

1. Fetch the repository. If you have Git on your computer, you can run ```git fetch https://github.com/jonahdf/pipeLite ```
2. Create a new Python environment. Note: You must have Python 3.3+ do use venv. You can follow the instructions for that [here](md/Python-Virtual-Environments.md).
3. Once inside your virtual environment, run ```conda install requirements.txt``` (assuming you are running Anaconda).
4. Everything else you need to know should be in src/main.ipynb!
**The only configuration you must set is your astrometry.net API key. Set this in the /src/dconf_dah.txt file**

<h2> Running in Colab</h2>

1. In [Colab](https://colab.research.google.com), open File -> Open Notebook -> Github, and paste in the [link to this repo](https://github.com/jonahdf/pipeLite), and open pipeLite-colab.ipynb

2. Make a copy with File -> Save a copy in Drive
3. When you run the first cell, it will fetch this Github repo and add it to your working directory. From there, you can add your astrometry.net API key.

<h2> Functions Implemented </h2>

```
"""
batch_process(datapaths, outfolder, darkfolder, flatfolder, biasfolder)
    Batch process HDR images.
    -Information for which darks,flats,biases is available here: https://wiki.uchicago.edu/display/2HA/220118
    -The datapaths must contain the 2 dynamic range files for every image, i.e bin1H and bin1L in the filenames
    Inputs:
        -datapaths (optional list str): list of paths to data files (default: datapath from input.py)
        -outfolder (optional str): path to output folder (default: outpath from input.py)
        -darkfolder (optional str): path to dark folder (default: darkpath from input.py)
        -flatfolder (optional str): path to flat folder (default: flatpath from input.py)
        -biasfolder (optional str): path to bias folder (default: biaspath from input.py)
    Outputs:
        Saves HDR images to output path.
""" 
```
```
"""
run_pipeline(folder)
    Runs the pipeline with the hotpix and astrometry steps
    on all the files in the directory.

    Inputs:
        -folder (str): The folder containing the files to be processed. (Default: outpath from setup.py)
        Note: Only HDR images will be processed by the pipeline

    Outputs:
        Saves HPX and WCS files in the outpath.
"""
```
```
"""
drizzle(filter, infolder, outfolder):
    Drizzles all WCS images in the output folder with given filters
    inputs: 
        - filter (optional lst str): list of filter strings to be drizzled. Will combine all images containing each filter into a single image.
        If not specified, will combine all images in the input folder
        - infolder (str): path to folder containing images to stack (default: outpath set in setup.py)
        Note: Infolder is set to outpath by default so it can take the output of the previous steps
        - outfolder (str): path to output folder (default: outpath set in setup.py)
    outputs:
        - Saves drizzled images to outpath
    example:
        drizzle(filter=['g-band','i-band'], infolder = out, outfolder = out) drizzles all images containing the strings "g-band" and "i-band" into their respective files
"""
```
