.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-04-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

scientimate.wavefromsurfaceelevzcross
=====================================

.. code:: python

    Hs, Ts, Hz, Tz, Hrms, H, T = scientimate.wavefromsurfaceelevzcross(Eta, fs, dispout='no')

Description
-----------

Calculate wave properties from water surface elevation by using an upward zero crossing method

Inputs
------

Eta
    Water surface elevation time series data in (m)
fs
    Sampling frequency that data collected at in (Hz)
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

Hs
    Significant Wave Height (m)
Ts
    Significant Wave Period (second)
Hz
    Zero Crossing Mean Wave Height (m)
Tz
    Zero Crossing Mean Wave Period (second)
Hrms
    Root Mean Square Wave Height (m)
H
    Wave Height Data Series array (m)
T
    Wave Period Data Series array (second)

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np
    import scipy as sp
    from scipy import signal

    fs=2 #Sampling frequency
    duration=1024 #Duration of the data
    N=fs*duration #Total number of points
    df=fs/N #Frequency difference 
    dt=1/fs #Time difference, dt=1/fs
    t=np.linspace(0,duration-dt,N) #Time
    Eta=sp.signal.detrend(0.5*np.cos(2*np.pi*0.2*t)+(-0.1+(0.1-(-0.1)))*np.random.rand(N))
    Hs,Ts,Hz,Tz,Hrms,H,T=sm.wavefromsurfaceelevzcross(Eta,fs,'yes')

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
