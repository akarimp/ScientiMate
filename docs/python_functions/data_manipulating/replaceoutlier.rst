.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2020-05-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

scientimate.replaceoutlier
==========================

.. code:: python

    xReplaced, outlier_Indx = scientimate.replaceoutlier(x, WindowSize=15, zscore_threshold=2, interpMethod='linear', dispout='no')

Description
-----------

Remove outliers in the time series using moving z-score window

Inputs
------

x
    Input data
WindowSize=15
    Window size (number of adjacent elements) that is used for moving window, should be equal or larger than 3
zscore_threshold=2
    | z-score threshold to define outliers
    | data in range of x < (xmean-zscore_threshold*std) or x > (xmean+zscore_threshold*std) considered outliers
interpMethod='linear'
    | Interpolation method for replacing spike points:
    | Matlab/Octave: 'linear', 'nearest', 'next', 'previous', 'pchip', 'cubic', 'spline'
    | Python/Scipy : 'linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic'
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

xReplaced
    Replaced data
outlier_Indx
    Logical index of replaced points

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np
    import scipy as sp
    from scipy import signal

    fs=128
    t=np.linspace(0,9.5,10*fs)
    x=np.sin(2*np.pi*0.3*t)+0.1*np.sin(2*np.pi*4*t)
    spikeloc=np.arange(10,len(t),100)
    x[np.int64(spikeloc+np.round(2*np.random.randn(len(spikeloc))))]=np.sign(np.random.randn(len(spikeloc)))
    x[220:225]=1.5
    x=x+5
    xReplaced,outlier_Indx=sm.replaceoutlier(x,37,2,'linear','yes')

    fs=2
    t=np.linspace(0,1023.5,1024*fs)
    x=sp.signal.detrend(0.5*np.cos(2*np.pi*0.2*t)+(-0.1+(0.1-(-0.1)))*np.random.rand(1024*fs))
    spikeloc=np.arange(10,len(t),100)
    x[np.int64(spikeloc+np.round(2*np.random.randn(len(spikeloc))))]=np.sign(np.random.randn(len(spikeloc)))
    xReplaced,outlier_Indx=sm.replaceoutlier(x,21,2,'linear','yes')

References
----------


.. License & Disclaimer
.. --------------------
..
.. Copyright (c) 2020 Arash Karimpour
..
.. http://www.arashkarimpour.com
..
.. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
.. IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
.. FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
.. AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
.. LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
.. OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
.. SOFTWARE.
