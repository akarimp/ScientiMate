.. YA LATIF

ScientiMate
===========

ScientiMate is a **coastal** and **ocean** data analysis library.
ScientiMate is a collection of functions and tools (in **MATLAB** and **Python**) developed for **metocean**, oceanography, coastal engineering, and earth science data analysis.
The examples of ScientiMate's applications are (but not limited to) ocean wave analysis, time series data analysis, signal processing, wind engineering, etc.

:Name: ScientiMate
:Description: Coastal and Ocean Data Analysis Library
:Version: 2.0
:Requirements: MATLAB or GNU Octave | Python (3 or later)
:Developer: Arash Karimpour | https://www.arashkarimpour.com
:Documentation: https://scientimate.readthedocs.io
:Tutorial Video: `YouTube Playlist <https://www.youtube.com/playlist?list=PLcrFHi9M_GZRTCshcgujlK7y5ZPim6afM>`_
:Source Code: https://github.com/akarimp/ScientiMate
:Report Issues: https://github.com/akarimp/ScientiMate/issues

MATLAB Version
==============

Installing (MATLAB Version)
---------------------------

To use MATLAB version of ScientiMate library:

* Install MATLAB or GNU Octave

    * MATLAB: https://www.mathworks.com
    * GNU Octave: https://octave.org

* Download ScientiMate:

    * Version 1.1 (GitHub): https://github.com/akarimp/ScientiMate/releases/download/2.0/scientimate.zip

* Unzip ScientiMate in any location you choose such as "C:\\"
* Open MATLAB or GNU Octave
* Add ScientiMate folder to MATLAB or GNU Octave path

Add ScientiMate folder to MATLAB or GNU Octave path
---------------------------------------------------

* Open MATLAB or GNU Octave
* Change a current folder (working directory) to a folder that contains ScientiMate files, for example "C:\\scientimate", in MATLAB or GNU Octave.
* Run a file named ``add_scientimate_to_path.m`` in MATLAB or GNU Octave to add ScientiMate folder to MATLAB or GNU Octave path.

Required Package for MATLAB
---------------------------

MATLAB users may need to install additional MATLAB Toolboxes such as Signal Processing Toolbox for some functions.

Required Package for GNU Octave
-------------------------------

GNU Octave users may need to install/load additional packages such as GNU Octave Signal package for some functions.
To find the list of the GNU Octave's pre-installed packages, run the following command in the Command Window:

.. code:: octave
    
    >> pkg list

For example, GNU Octave comes with Signal package but it needs to loaded every time GNU Octave starts. The Signal package can be loaded inside GNU Octave by running the following command in the Command Window (This should be done every time GNU Octave is opened):

.. code:: octave
    
    >> pkg load signal

If GNU Octave Signal Package is not already installed, it should be first installed from https://packages.octave.org, and then get loaded by running the following commands in the Command Window:

.. code:: octave

    >> pkg install "https://downloads.sourceforge.net/project/octave/Octave%20Forge%20Packages/Individual%20Package%20Releases/signal-1.4.5.tar.gz"
    >> pkg load signal

Quick Start (MATLAB Version)
----------------------------

.. code:: matlab

    x(:,1)=linspace(1,10,10);
    y(:,1)=1+rand(10,1);
    y(:,2)=2+rand(10,1);
    plot2d(x,y,'line_confid','blue_red','large')

Python Version
==============

Installing (Python Version)
---------------------------

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

Required Package for Python
---------------------------

Following packages are required:

* NumPy (https://numpy.org)
* SciPy (https://www.scipy.org)
* pandas (https://pandas.pydata.org)
* Matplotlib (https://matplotlib.org)

Quick Start (Python Version)
----------------------------

.. code:: python

    import scientimate as sm
    import numpy as np

    print(sm.__version__)

    x=np.linspace(1,10,10)
    y=np.zeros((10,2))
    y[:,0]=1+np.random.rand(10)
    y[:,1]=2+np.random.rand(10)
    sm.plot2d(x,y,'line_confid','blue_red','large')

About
=====

Operating System
----------------

ScientiMate code can be run on Microsoft Windows, Mac, and Linux. However, make sure any given path is compatible with a running operating system. In particular, "\\" is used in Windows path, while "/" is used in Mac or Linux path. For example, if a path is "C:\\" on Windows machine, it would be "C:/" on Mac or Linux.

Required Programing Language
----------------------------

This library can be run by using MATLAB (https://www.mathworks.com), GNU Octave (https://octave.org), or Python (https://www.python.org). 

Citation
--------

Cite ScientiMate as:

Karimpour, A. (2023). ScientiMate, Coastal and Ocean Data Analysis Library (Version 2.0) [Computer software]. https://github.com/akarimp/ScientiMate

Recommended Books
-----------------

.. list-table::
   :header-rows: 1
   :align: center

   * - .. figure:: figures/Figure_Book_Coastal.jpg
     - .. figure:: figures/Figure_Book_Matlab.jpg
     - .. figure:: figures/Figure_Book_Python.jpg

   * - | **Ocean Wave Data Analysis**
       | Introduction to Time Series Analysis, Signal Processing, and Wave Prediction.
       |
       | Order at Amazon: https://www.amazon.com/dp/0692109978
       |
       | Read Online: https://github.com/akarimp/Ocean-Wave-Data-Analysis
     - | **Fundamentals of Data Science with MATLAB**
       | Introduction to Scientific Computing, Data Analysis, and Data Visualization.
       |
       | Order at Amazon: https://www.amazon.com/dp/1735241016
       |
       | Read Online: https://github.com/akarimp/Fundamentals-of-Data-Science-with-MATLAB
     - | **Principles of Data Science with Python**
       | Introduction to Scientific Computing, Data Analysis, and Data Visualization.
       |
       | Order at Amazon: https://www.amazon.com/dp/1735241008
       |
       | Read Online: https://github.com/akarimp/Principles-of-Data-Science-with-Python

Recommended Applications
------------------------

.. list-table::
   :header-rows: 1
   :align: center

   * - .. figure:: figures/Figure_Oceanlyz_Logo.png
     - .. figure:: figures/Figure_ScientiMate_Logo.png
     - .. figure:: figures/Figure_AsanPlot_Screenshot.jpg

   * - | **OCEANLYZ**
       | Ocean Wave Analyzing Toolbox
       |
       | Download: https://github.com/akarimp/Oceanlyz
     - | **ScientiMate**
       | Coastal and Ocean Data Analysis Library
       |
       | Download: https://github.com/akarimp/ScientiMate
     - | **AsanPlot**
       | Data cleaning and plotting software
       |
       | Download: https://github.com/akarimp/AsanPlot

License Agreement and Disclaimer
--------------------------------

ScientiMate: Coastal and Ocean Data Analysis Library

Copyright (c) 2023 Arash Karimpour

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
