<h1>Python Virtual Environments</h1>

When working on a large python project, you often need to install third party packages. The larger a project, the more packages you need, often with sub-dependencies. It can quickly get quite complicated. 

Third-party software gets updated periodically, and it is not uncommon for projects to break due to updated packages. This is why it can be helpful to take a snapshot of the exact versions of each third-party package that a project uses. You will often find something that looks like this, called requirements.txt in Python projects:

![requirements.txt](images/requirements.txt.jog.jpg)

When older or newer versions of packages are already installed onto machines, code can often break. That is why it is helpful to use a tool called **virtual environments**. When you create a virtual environment, you are installing a python environment that is a completely blank slate, with no packages installed. It is common to create separate virtual environments for each software project you run locally on your computer. With that blank slate, you can then install the exact version of each package that the developer has on their computer, ensuring that the code will run correctly.

Luckily, Python 3.3+ comes preinstalled with a virtual environment manager called venv. This tutorial will show you how to get started with it. Note that I am using a Windows machine, so running these steps on a Mac might be a bit complicated.

[1] Install [Python](https://www.python.org/ftp/python/3.10.5/python-3.10.5-amd64.exe) (Python v3.3 or later is required).

[1b] For Windows, you will need to install the [Microsoft Build Tools for Visual Studio](https://wiki.python.org/moin/WindowsCompilers#Microsoft_Visual_C.2B-.2B-_14.2_standalone:_Build_Tools_for_Visual_Studio_2019_.28x86.2C_x64.2C_ARM.2C_ARM64.29)
* In Build tools, ensure the latest versions of MSVCvXXX - VS YYYY C++ x64/x86 build tools and Windows 10 SDK are checked.

[2] Ensure that Python v3.3+ is installed:

```
python --version
```

[3] Open the a command prompt and navigate to the folder in which you want to create your virtual environment. 

[4] Create the virtual environment:

```
python -m venv [name]
```

where the [name] is the name of your environment, e.g "astro" This will create a sub-folder (e.g., pipelite) containing your environment.

[4] Next, you can start running your new environment in the terminal with the following command:

```
[name]\Scripts\activate.bat
```

Now, you can tell you are in your new environment if the name appears in parentheses in your terminal:

![python command line](images/python_environment.jpg)

[5] Next, if you want your new python environment to be recognized by Jupyter Notebook, you must install the ipykernel package to your new environment:

```
pip install ipykernel
```

[6] Once this is done, add the environment to your Jupyter notebook path with the following:

```
python -m ipykernel install --user --name=[name]
```

[7] Install the required packages for your project, running this command in your virtual environment. NOTE: For Windows, you made need to install [Microsoft C++ Build Tools Visual Studio 2022]([https://visualstudio.microsoft.com/downloads/](https://visualstudio.microsoft.com/thank-you-downloading-visual-studio/?sku=BuildTools&rel=16).

```
pip install -r path/to/requirements.txt
```
