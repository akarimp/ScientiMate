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

scientimate.wavefromvelocityzcross
==================================

.. code:: python

    Hs, Ts, Hz, Tz, Hrms, H, T, Eta, t = scientimate.wavefromvelocityzcross(Ux, Uy, fs, h, heightfrombed=0, Kuvmin=0.15, kCalcMethod='beji', dispout='no')

Description
-----------

Calculate wave properties from wave orbital velocity by using an upward zero crossing method

Inputs
------

Ux
    Wave horizontal orbital velocity data in x direction in (m/s)
Uy
    Wave horizontal orbital velocity data in y direction in (m/s)
fs
    Sampling frequency that data collected at in (Hz)
h
    Water depth in (m)
heightfrombed=0
    Height from bed that data collected at in (m)
Kuvmin=0.15
    | Minimum acceptable value for an orbital velocity converstion factor
    | If Kuvmin=0.15, it avoid wave amplification larger than 6 times (1/0.15)
kCalcMethod='beji'
    | Wave number calculation method 
    | 'hunt': Hunt (1979), 'beji': Beji (2013), 'vatankhah': Vatankhah and Aghashariatmadari (2013) 
    | 'goda': Goda (2010), 'exact': calculate exact value 
Rho=1000
    Water density (kg/m^3)
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
Eta
    Water surface elevation time series in (m)
t
    Time (s)

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
    hfrombed=4
    h=5
    k=0.2
    Ux=(np.pi/5)*(2*Eta)*(np.cosh(k*hfrombed)/np.sinh(k*h)) 
    Uy=0.2*Ux
    Hs,Ts,Hz,Tz,Hrms,H,T,Eta,t=sm.wavefromvelocityzcross(Ux,Uy,fs,5,4,0.15,'beji','yes')

References
----------

Beji, S. (2013). 
Improved explicit approximation of linear dispersion relationship for gravity waves. 
Coastal Engineering, 73, 11-12.

Goda, Y. (2010). 
Random seas and design of maritime structures. 
World scientific.

Hunt, J. N. (1979). 
Direct solution of wave dispersion equation. 
Journal of the Waterway Port Coastal and Ocean Division, 105(4), 457-459.

Vatankhah, A. R., & Aghashariatmadari, Z. (2013). 
Improved explicit approximation of linear dispersion relationship for gravity waves: A discussion. 
Coastal engineering, 78, 21-22.

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
