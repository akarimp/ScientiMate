.. YA LATIF

ScientiMate
===========

ScientiMate is a library for earth-science data analysis. This library can be used for wide range of data analysis including a time series analysis, signal processing, and geo-data calculation.

:Name: ScientiMate
:Description: Earth-Science Data Analysis Library
:Version: 1.1.2
:Requirements: MATLAB, or GNU Octave, or Python (3 or later)
:Developer: Arash Karimpour (http://www.arashkarimpour.com)
:Documentation: https://scientimate.readthedocs.io
:Tutorial Video: `YouTube Playlist <https://www.youtube.com/playlist?list=PLcrFHi9M_GZRTCshcgujlK7y5ZPim6afM>`_
:Source Code: https://github.com/akarimp/scientimate
:Report Issues: https://github.com/akarimp/ScientiMate/issues


Installation (MATLAB Version)
-----------------------------

To use MATLAB version of ScientiMate library:

* Install MATLAB or GNU Octave
* Download ScientiMate:

    * Version 1.1 (GitHub): https://github.com/akarimp/ScientiMate/releases/download/1.1/scientimate.zip

* Unzip ScientiMate in any location you choose such as "C:\\"
* Open MATLAB or GNU Octave
* Add ScientiMate folder to MATLAB or GNU Octave path

**Add ScientiMate folder to MATLAB or GNU Octave path using add_scientimate_to_path.m**

* Open MATLAB or GNU Octave
* Change a current folder (working directory) to a folder that contains ScientiMate files, for example "C:\\scientimate", in MATLAB or GNU Octave.
* Run a file named add_scientimate_to_path.m in MATLAB or GNU Octave to add ScientiMate folder to MATLAB or GNU Octave path.


Installation (Python Version)
-----------------------------

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


Required Programing Language
----------------------------

This library can be run by using MATLAB (https://www.mathworks.com), GNU Octave (https://www.gnu.org/software/octave), or Python (https://www.python.org). 


Required Package for MATLAB
---------------------------

MATLAB users may need to install additional MATLAB Toolboxes such as Signal Processing Toolbox for some functions.


Required Package for GNU Octave
-------------------------------

GNU Octave users may need to install/load additional packages such as GNU Octave Signal Package for some functions.

For example, GNU Octave Signal Package can be loaded inside GNU Octave by using a following command in a command window (This should be done every time GNU Octave is opened):


.. code:: octave
    
    >> pkg load signal


If GNU Octave Signal Package is not already installed, it should be first installed from Octave Forge (https://octave.sourceforge.io), and then get loaded by using the following commands in a command window:

.. code:: octave

    >> pkg install -forge signal
    >> pkg load signal


Required Package for Python
---------------------------

Following packages are required:

* NumPy (https://numpy.org)
* SciPy (https://www.scipy.org)
* Matplotlib (https://matplotlib.org)


Quick Start (MATLAB Version)
----------------------------

.. code:: matlab

    x(:,1)=linspace(1,10,10);
    y(:,1)=1+rand(10,1);
    y(:,2)=2+rand(10,1);
    plot2d(x,y,'line_confid','blue_red','large')


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


Recommended Books
-----------------

* | **Ocean Wave Data Analysis**
  | Introduction to Time Series Analysis, Signal Processing, and Wave Prediction.
  | Order at Amazon: https://www.amazon.com/dp/0692109978
  |
* | **Principles of Data Science with Python**
  | Introduction to Scientific Computing, Data Analysis, and Data Visualization.
  | Order at Amazon: https://www.amazon.com/dp/1735241008
  |
* | **Fundamentals of Data Science with MATLAB**
  | Introduction to Scientific Computing, Data Analysis, and Data Visualization.
  | Order at Amazon: https://www.amazon.com/dp/1735241016

Citation
--------

Cite this package as:

Karimpour, A. (2022). ScientiMate, Earth-Science Data Analysis Library.


License Agreement and Disclaimer
--------------------------------

ScientiMate: Earth-Science Data Analysis Library

Copyright (c) 2022 Arash Karimpour

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
