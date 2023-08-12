.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-02-01/2020-02-01                  +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

scientimate.replacemissing1d
============================

.. code:: python

    xReplaced, NaN_Indx = scientimate.replacemissing1d(x, what2replace='both', interpMethod='linear', dispout='no')

Description
-----------

Replace missing data points in 1d data such as time series

Inputs
------

x
    Input data
what2replace='both'
    | What needs to be replaced
    | 'NaN': replacing NaN data points
    | 'Inf': replacing Inf data points
    | 'both': replacing NaN and Inf data points
    | Number: replacing data points equal to Number
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
NaN_Indx
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
    x[np.int64(spikeloc+np.round(2*np.random.randn(len(spikeloc))))]=np.nan
    x=x+5
    x[220:225]=np.nan
    xReplaced,NaN_Indx=sm.replacemissing1d(x,'NaN','linear','yes')

    fs=2
    t=np.linspace(0,1023.5,1024*fs)
    x=sp.signal.detrend(0.5*np.cos(2*np.pi*0.2*t)+(-0.1+(0.1-(-0.1)))*np.random.rand(1024*fs))
    spikeloc=np.arange(10,len(t),100)
    x=x+1
    x[np.int64(spikeloc+np.round(2*np.random.randn(len(spikeloc))))]=0
    xReplaced,NaN_Indx=sm.replacemissing1d(x,0,'linear','yes')

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
