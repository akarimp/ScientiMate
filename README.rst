.. YA LATIF

ScientiMate
===========

ScientiMate is a library for earth-science data analysis. This library can be used for wide range of data analysis including a time series analysis, signal processing, and geo-data calculation.

:Name: ScientiMate
:Description: Earth-Science Data Analysis Library
:Version: 1.0
:Requirements: Python (3 or later), NumPy, SciPy, Matplotlib
:Developer: Arash Karimpour (http://www.arashkarimpour.com)
:Documentation: https://scientimate.readthedocs.io
:Tutorial Video: `YouTube Playlist <https://www.youtube.com/playlist?list=PLcrFHi9M_GZRTCshcgujlK7y5ZPim6afM>`_
:Source Code: https://github.com/akarimp/scientimate
:Report Issues: https://github.com/akarimp/scientimate/issues


Installation
------------

Available soon

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

Karimpour, A. (2020). ScientiMate, Earth-Science Data Analysis Library.


License Agreement and Disclaimer
--------------------------------

ScientiMate: Earth-Science Data Analysis Library

Copyright (c) 2020 Arash Karimpour

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
