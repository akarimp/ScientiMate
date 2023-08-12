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

scientimate.pressure2surfaceelevzcross
======================================

.. code:: python

    Eta, t = scientimate.pressure2surfaceelevzcross(P, fs, h, heightfrombed=0, Kpmin=0.15, kCalcMethod='beji', Rho=1000, dispout='no')

Description
-----------

Calculate water surface elevation time series from water pressure time series by using an upward zero crossing method

Inputs
------

P
    Water pressure time series data in (N/m^2)
fs
    Sampling frequency that data collected at in (Hz)
h
    Water depth in (m)
heightfrombed=0
    Height from bed that data collected at in (m)
Kpmin=0.15
    | Minimum acceptable value for a pressure response factor
    | If Kpmin=0.15, it avoid wave amplification larger than 6 times (1/0.15)
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
    P=Eta*9.81*1000*(np.cosh(k*hfrombed)/np.cosh(k*h))
    Eta,t=sm.pressure2surfaceelevzcross(P,fs,5,4,0.15,'beji',1025,'yes')

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
