Getting Started (Python)
========================

In order to use Python version of ScientiMate library, first, Python programming language, and then, the ScientiMate library should be installed.


Installing
----------

To use Python version of ScientiMate library:

* Install Python
* Install ScientiMate

**1) Install Python**

First, you need to install Python programming language.

* Method 1:
    Install Python from https://www.python.org and then use the **pip** command to install required packages
* Method 2 (Recommended):
    Install Anaconda Python distribution from https://www.anaconda.com and then use the **conda** command to install required packages

**2) Install ScientiMate**

After Python is installed, you need to install ScientiMate library.

If you installed Python, then you need to install ScientiMate via pip (https://pypi.org/project/scientimate). To do that, open the Command Prompt (or Terminal) and run:

.. code:: python

    pip install scientimate

If you installed Anaconda Python distribution, then you need to install ScientiMate via Anaconda cloud (https://anaconda.org/akarimp/scientimate). To do that, open the Anaconda Prompt and run:

.. code:: python

     conda install -c akarimp scientimate


Operating System
----------------

This code can be run on Microsoft Windows, Mac, and Linux. However, make sure any given path is compatible with a running operating system. In particular, "\\" is used in Windows path, while "/" is used in Mac or Linux path. For example, if a path is "C:\\" on Windows machine, it would be "C:/" on Mac or Linux.


Required Programing Language
----------------------------

The Python version of this library can be run by using Python 3 or later (https://www.python.org or https://www.anaconda.com).


Required Package for Python
---------------------------

Following packages are required:

* NumPy (https://numpy.org)
* SciPy (https://www.scipy.org)
* pandas (https://pandas.pydata.org)
* Matplotlib (https://matplotlib.org)


Quick Start
-----------

.. code:: python

    import scientimate as sm
    import numpy as np

    print(sm.__version__)

    x=np.linspace(1,10,10)
    y=np.zeros((10,2))
    y[:,0]=1+np.random.rand(10)
    y[:,1]=2+np.random.rand(10)
    sm.plot2d(x,y,'line_confid','blue_red','large')
