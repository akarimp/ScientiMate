.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2018-05-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

scientimate.wavespectrum2timeseries
===================================

.. code:: python

    Eta, t, Hm0, fp, fEta, SxxEta, a, w, Phi = scientimate.wavespectrum2timeseries(f, Sxx, fs=2, dispout='no')

Description
-----------

| Generate random water wave data from a given water wave spectrum using wave superposition
| For more options use psd2timeseries

Inputs
------

f
    Frequency (Hz)
Sxx
    Wave power spectral density (m^2s)
fs=2
    Sampling frequency that data are collected at in (Hz)
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

Eta
    Water surface level time series in (m)
t
    Time in (s)
Hm0
    Zero moment wave height (m)
fp
    Peak wave frequency (Hz), fp=1/Tp (Tp: Peak wave period (s))
fEta
    Frequency from generated time series(Hz)
SxxEta
    Power spectral density from generated time series (m^2s)
a
    Wave amplitude for for one-sided spectrum (0<fEta<fs/2) from generated time series (m)
w
    Wave angular frequency for for one-sided spectrum (0<fEta<fs/2) from generated time series (rad/s)
Phi
    Wave random phase for for one-sided spectrum (0<fEta<fs/2) from generated time series (rad)

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np

    N=2**11
    fs=8
    df=fs/N #Frequency difference 
    f=np.arange(0,fs/2+df,df) #Frequency vector 
    f[0]=f[1]/2 #Assign temporarily non-zero value to fisrt element of f to prevent division by zero
    Sxx=0.016*9.81**2/((2*np.pi)**4*(f**5))*np.exp(-1.25*(0.33/f)**4) #Calculating Spectrum 
    f[0]=0
    Sxx[0]=0
    Eta,t,Hm0,fp,fEta,SxxEta,a,w,Phi=sm.wavespectrum2timeseries(f,Sxx,fs,'yes')

References
----------

Branlard, E. (2010).
Generation of time series from a spectrum.
Technical University Denmark. National Laboratory for Sustainable Energy.

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
