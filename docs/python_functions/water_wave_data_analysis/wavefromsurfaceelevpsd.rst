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

scientimate.wavefromsurfaceelevpsd
==================================

.. code:: python

    Hm0, fp, Tp, Tm01, Tm02, f, Syy, m0 = scientimate.wavefromsurfaceelevpsd(Eta, fs, fcL=0, fcH=None, nfft=None, SegmentSize=256, OverlapSize=128, dispout='no')

Description
-----------

Calculate wave properties from water surface elevation power spectral density

Inputs
------

Eta
    Water surface elevation time series data in (m)
fs
    Sampling frequency that data collected at in (Hz)
fcL=0
    Low cut-off frequency, between 0*fs to 0.5*fs (Hz)
fcH=fs/2
    High cut-off frequency, between 0*fs to 0.5*fs (Hz)
nfft=length(Eta)
    Total number of points between 0 and fs that spectrum reports at is (nfft+1)
SegmentSize=256
    Segment size, data are divided into the segments each has a total element equal to SegmentSize
OverlapSize=128
    Number of data points that are overlaped with data in previous segments 
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

Hm0
    Zero-Moment Wave Height (m)
fp
    Peak wave frequency (Hz)
Tp
    Peak wave period (second)
Tm01
    Wave Period from m01 (second), Mean Wave Period
Tm02
    Wave Period from m02 (second), Mean Zero Crossing Period
f
    Frequency (Hz)
Syy
    Power spectral density (m^2/Hz)
m0
    Zero-Moment of the power spectral density (m^2)

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
    Hm0,fp,Tp,Tm01,Tm02,f,Syy,m0=sm.wavefromsurfaceelevpsd(Eta,fs,0,fs/2,N,256,128,'yes')

References
----------

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
