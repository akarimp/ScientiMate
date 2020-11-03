Getting Started
===============

In order to use ScientiMate package, first, Python programming language, and then, the ScientiMate package should be installed.

Installation
------------

To use ScientiMate package:

* Install Python
* Install ScientiMate

**1) Install Python**

First, we need to install Python programming language.

* Method 1:
    Install pure Python from https://www.python.org and then use the **pip** command to install required packages
* Method 2 (Recommended):
    Install Anaconda Python distribution from https://www.anaconda.com and then use the **conda** command to install required packages

**2) Install ScientiMate**

After Python is installed, we need to install ScientiMate package.

To install ScientiMate via pip (https://pypi.org/project/scientimate):

.. code:: python

    pip install scientimate

To install ScientiMate via Anaconda cloud (https://anaconda.org/akarimp/scientimate):

.. code:: python

     conda install -c akarimp scientimate


Operating System
----------------

This code can be run on Windows, Mac, and Linux.


Required Programing Language
----------------------------

This package can be run by using Python 3 or later (https://www.python.org or https://www.anaconda.com).


Required Package for Python
---------------------------

Following packages are required:

* NumPy (https://numpy.org)
* SciPy (https://www.scipy.org)
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

