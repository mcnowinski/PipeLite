<h1>Python Virtual Environments</h1>

When working on a large python project, you often need to install third party packages. The larger a project, the more packages you need, often with sub-dependencies. It can quickly get quite complicated. 

Third-party software gets updated periodically, and it is not uncommon for projects to break due to updated packages. This is why it can be helpful to take a snapshot of the exact versions of each third-party package that a project uses. You will often find something that looks like this, called requirements.txt in Python projects:

![requirements.txt](images/requirements.txt.jog.jpg)

When older or newer versions of packages are already installed onto machines, code can often break. That is why it is helpful to use a tool called **virtual environments**. When you create a virtual environment, you are installing a python environment that is a completely blank slate, with no packages installed. It is common to create separate virtual environments for each software project you run locally on your computer. With that blank slate, you can then install the exact version of each package that the developer has on their computer, ensuring that the code will run correctly.

Luckily, Python 3.3+ comes preinstalled with a virtual environment manager called venv. This tutorial will show you how to get started with it. Note that I am using a Windows machine, so running these steps on a Mac might be a bit complicated.

[1] Install [Anaconda](https://repo.anaconda.com/archive/Anaconda3-2022.05-Windows-x86_64.exe) (Python v3.3 or later is required). 

[2] Open the Anaconda Prompt and navigate to the folder in which you want to create your virtual environment. 

[3] Create the virtual environment:

```
conda create --name pipelite
```

where the last argument is the name of your environment, e.g "pipelite." This will create a folder containing your environment.

[4] Next, you can start running your new environment in the terminal with the following command:

```
conda activate pipelite
```

Now, you can tell you are in your new environment if the name appears in parentheses in your terminal:

![python command line](images/python_environment.jpg)

[5] Next, if you want your new python environment to be recognized by Jupyter Notebook, you must install the ipykernel package to your new environment:

```
conda install ipykernel
```

[6] Once this is done, add the environment to your Jupyter notebook path with the following:

```
python -m ipykernel install --user --name=pipelite
```

[7] Install the required packages for your project, running this command in your virtual environment.

```
conda install --file path/to/requirements.txt
```
