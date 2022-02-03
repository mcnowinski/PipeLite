<h1>Welcome to PipeLite, formally known as the Pretty Picture Processing Pipeline for ASTR 21200</h1>

The goal of this repo is to abstract as much behind-the-scenes functionality as possible to make the image processing pipeline simple and easy to run

<h2> Running in Colab</h2>

1. In [Colab](https://colab.research.google.com), open File -> Open Notebook -> Github, and paste in the [link to this repo](https://github.com/jonahdf/PPPP), and open PPPP-colab.ipynb

2. Make a copy with File -> Save a copy in Drive
3. When you run the first cell, it will fetch this Github repo and add it to your working directory. From there, you can add your astrometry.net API key.
To get started, clone the repo and open the main.py notebook. It explains how to run the pipeline.

<h2> Running Locally </h2>

1. Fetch the repository. If you have Git on your computer, you can run ```git fetch https://github.com/jonahdf/PPPP ```
2. Create a new Python environment. Note: You must have Python 3.3+ do use venv. You can follow the instructions for that [here](https://wiki.uchicago.edu/display/2HA/Python+Virtual+Environments)
3. Once inside your virtual environment, run ```pip install -r requirements.txt``` or ```conda install requirements.txt``` if using Conda.
4. Everything else you need to know should be in src/main.ipynb!
**The only configuration you must set is your astrometry.net API key. Set this in the /src/dconf_dah.txt file**
