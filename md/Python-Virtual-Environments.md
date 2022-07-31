<h1>Python Virtual Environments</h1>

When working on a large python project, you often need to install third party packages. The larger a project, the more packages you need, often with sub-dependencies. It can quickly get quite complicated. 

Third-party software gets updated periodically, and it is not uncommon for projects to break due to updated packages. This is why it can be helpful to take a snapshot of the exact versions of each third-party package that a project uses. You will often find something that looks like this, called requirements.txt in Python projects:

![requirements.txt](images/requirements.txt.jog.jpg)

When older or newer versions of packages are already installed onto machines, code can often break. That is why it is helpful to use a tool called **virtual environments**. When you create a virtual environment, you are installing a python environment that is a completely blank slate, with no packages installed. It is common to create separate virtual environments for each software project you run locally on your computer. With that blank slate, you can then install the exact version of each package that the developer has on their computer, ensuring that the code will run correctly.

Luckily, Python 3.3+ comes preinstalled with a virtual environment manager called venv. This tutorial will show you how to get started with it. Note that I am using a Mac, so running these steps on Windows might be a bit complicated.

[1] Make sure you have the correct version of python installed. **If you are using Anaconda, skip to after step 3.**

  In your terminal of choice (the default one works well on Mac, and powershell/cmd can work on Windows, run the following command:

```
python3 --version
```

And make sure the output is 3.3 or above. If not, you will have to install Python 3, which you can learn to do [here](https://www.python.org/downloads/). 

[2] Navigate to the folder you want to create your virtual environment in. 

[3] Create and run the virtual environment:

```
python3 -m venv myenv
```

Where the last argument is the name of your environment, e.g "astro." This will create a folder containing your environment.

Next, you can start running your new environment in the terminal with the following command:

```
source myenv /bin/activate
```

Now, you can tell you are in your new environment if the name appears in parentheses in your terminal:

[python command line](images/python_environment.png)

**Note: If you are using anaconda, the same thing can be achieved in these steps from any directory:**

* Create conda env: `conda create --name myenv`

* Activate your environment: `conda activate myenv`

[4] Next, if you want your new python environment to be recognized by Jupyter Notebook, you must install the ipykernel package to your new environment:

```
pip install ipykernel
```

[5] Once this is done, add the environment to your Jupyter notebook path with the following:

```
python -m ipykernel install --user --name=myenv
```

Check to make sure your Jupyter recognizes your new environment, and you should be good to install packages!

You can install from a requirements.txt by running this command in your virtual environment.

```
pip install -r path /to/requirements.txt
```

or, if using Conda:

```
conda install --file path/to/requirements.txt
```
