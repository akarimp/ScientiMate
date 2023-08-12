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

scientimate.velocity2surfaceelevfft
===================================

.. code:: python

    Eta, t = scientimate.velocity2surfaceelevfft(U, fs, duration, h, heightfrombed=0, fKuvmin=None, fcL=0, fcH=None, KuvafterfKuvmin='constant', kCalcMethod='beji', dispout='no')

Description
-----------

Calculate water surface elevation time series from wave orbital velocity time series by using Fast Fourier Transform

Inputs
------

U
    Wave horizontal orbital velocity data in (m/s)
fs
    Sampling frequency that data collected at in (Hz)
duration
    Duration time that data are collected (second)
h
    Water depth in (m)
heightfrombed=0
    Height from bed that data collected at in (m)
fKuvmin=fs/2
    Frequency that a velocity conversion factor (Kuv) at that frequency is considered as a minimum limit for Kuv
fcL=0
    Low cut-off frequency, between 0*fs to 0.5*fs (Hz)
fcH=fs/2
    High cut-off frequency, between 0*fs to 0.5*fs (Hz)
KuvafterfKuvmin='constant'
    | Define conversion factor, Kuv, value for frequency larger than fKuvmin
    | 'nochange': Kuv is not changed for frequency larger than fKuvmin 
    | 'one': Kuv=1 for frequency larger than fKuvmin 
    | 'constant': Kuv for f larger than fKuvmin stays equal to Kuv at fKuvmin (constant)
kCalcMethod='beji'
    | Wave number calculation method 
    | 'hunt': Hunt (1979), 'beji': Beji (2013), 'vatankhah': Vatankhah and Aghashariatmadari (2013) 
    | 'goda': Goda (2010), 'exact': calculate exact value 
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
    U=(np.pi/5)*(2*Eta)*(np.cosh(k*hfrombed)/np.sinh(k*h)) 
    Eta1,t=sm.velocity2surfaceelevfft(U,fs,duration,5,4,0.6,0,fs/2,'constant','beji','yes')

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

Wiberg, P. L., & Sherwood, C. R. (2008). 
Calculating wave-generated bottom orbital velocities from surface-wave parameters. 
Computers & Geosciences, 34(10), 1243-1262.

Welch, P. (1967). 
The use of fast Fourier transform for the estimation of power spectra: a method based on time averaging over short, modified periodograms. 
IEEE Transactions on audio and electroacoustics, 15(2), 70-73.

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
