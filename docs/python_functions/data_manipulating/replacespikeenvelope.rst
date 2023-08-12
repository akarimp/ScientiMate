.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-02-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

scientimate.replacespikeenvelope
================================

.. code:: python

    xDespiked, Indx = scientimate.replacespikeenvelope(x, lowbound, upbound, nrpeat=1, interpMethod='linear', dispout='no')

Description
-----------

Remove spikes in the time series that are outside a defined envelope

Inputs
------

x
    Input data
lowbound
    Lower boundary of the data, all data points should be larger than that
upbound
    Upper boundary of the data, all data points should be smaller than that
nrpeat=1
    Number of time despiking procedure is repeating
interpMethod='linear'
    | Interpolation method for replacing spike points:
    | Matlab/Octave: 'linear', 'nearest', 'next', 'previous', 'pchip', 'cubic', 'spline'
    | Python/Scipy : 'linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic'
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

xDespiked
    Dispiked data
Indx
    Index of despiked points

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
    x[220:225]=1.5
    x=x+5
    xDespiked,Indx=sm.replacespikeenvelope(x,3.9,6.1,1,'linear','yes')

    fs=2
    t=np.linspace(0,1023.5,1024*fs)
    x=sp.signal.detrend(0.5*np.cos(2*np.pi*0.2*t)+(-0.1+(0.1-(-0.1)))*np.random.rand(1024*fs))
    spikeloc=np.arange(10,len(t),100)
    x[np.int64(spikeloc+np.round(2*np.random.randn(len(spikeloc))))]=np.sign(np.random.randn(len(spikeloc)))
    xDespiked,Indx=sm.replacespikeenvelope(x,-0.6,0.6,1,'linear','yes')

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
