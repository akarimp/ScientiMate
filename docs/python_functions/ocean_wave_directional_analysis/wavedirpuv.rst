.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-05-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

scientimate.wavedirpuv
======================

.. code:: python

    Wavedir, theta1, theta2, f = scientimate.wavedirpuv(P, Ux, Uy, fs, h, Pheightfrombed=0, UVheightfrombed=0, dirCalcMethod='puv1', coordinatesys='xyz', fmaxpcorr=None, fminpcorr=0, fKuvmin=None, fcL=0, fcH=None, fmaxpcorrCalcMethod='auto', Kpafterfmaxpcorr='constant', KuvafterfKuvmin='constant', kCalcMethod='beji', Rho=1000, nfft=None, SegmentSize=256, OverlapSize=128, dispout='no')

Description
-----------

Calculate wave direction using pressure and horizontal orbital velocity

Inputs
------

P
    Water pressure time series data in (N/m^2)
Ux
    Wave horizontal orbital velocity data in x direction in (m/s)
Uy
    Wave horizontal orbital velocity data in y direction in (m/s)
fs
    Sampling frequency that data collected at in (Hz)
h
    Water depth in (m)
Pheightfrombed=0
    Pressure sensor height from bed that data collected at in (m)
UVheightfrombed=0
    Velocity sensor height from bed that data collected at in (m)
dirCalcMethod='puv1'
    | Wave number calculation method 
    | 'puv1': PUV Method 1, 'puv2': PUV Method 2, 'puv3': PUV Method 3 
coordinatesys='xyz'
    | Define the coordinate system 
    | 'xyz': XYZ coordinate system, 'enu': ENU (East North Up) coordinate system 
    | If coordinatesys='enu', then x is East and y is North  
    | If coordinatesys='enu', results are reported with respect to true north  
    | In true north coordinate system, wave comes from as:
    | 0 degree: from north, 90 degree: from east, 180 degree: from south, 270 degree: from west  
fmaxpcorr=fs/2
    | Maximum frequency that a pressure attenuation factor applies up on that (Hz)
    | If fmaxpcorrCalcMethod='user', then the smaller of calculated and user defined fmaxpcorr will be chosen
fminpcorr=0
    | Minimum frequency that is used for defining fmaxpcorr if fmaxpcorrCalcMethod='auto' (Hz)
    | fminpcorr should be smaller than fp 
    | If swell energy exists, fminpcorr should be smaller than fp of wind sea (fpsea) and larger than fp of swell (fpswell) if there swell 
fKuvmin=fs/2
    Frequency that a velocity conversion factor (Kuv) at that frequency is considered as a minimum limit for Kuv
fcL=0
    Low cut-off frequency, between 0*fs to 0.5*fs (Hz)
fcH=fs/2
    High cut-off frequency, between 0*fs to 0.5*fs (Hz)
fmaxpcorrCalcMethod='auto'
    | Define if to calculate fmaxpcorr and ftail or to use user defined
    | 'user': use user defined value for fmaxpcorr
    | 'auto': automatically define value for fmaxpcorr
Kpafterfmaxpcorr='constant'
    | Define a apressure response factor, Kp, value for frequency larger than fmaxpcorr
    | 'nochange': Kp is not changed for frequency larger than fKuvmin 
    | 'one': Kp=1 for frequency larger than fmaxpcorr 
    | 'constant': Kp for f larger than fmaxpcorr stays equal to Kp at fmaxpcorr (constant)
KuvafterfKuvmin='constant'
    | Define conversion factor, Kuv, value for frequency larger than fKuvmin
    | 'nochange': Kuv is not changed for frequency larger than fKuvmin 
    | 'one': Kuv=1 for frequency larger than fKuvmin 
    | 'constant': Kuv for f larger than fKuvmin stays equal to Kuv at fKuvmin (constant)
kCalcMethod='beji'
    | Wave number calculation method 
    | 'hunt': Hunt (1979), 'beji': Beji (2013), 'vatankhah': Vatankhah and Aghashariatmadari (2013) 
    | 'goda': Goda (2010), 'exact': calculate exact value 
Rho=1000
    Water density (kg/m^3)
nfft=length(P)
    Total number of points between 0 and fs that spectrum reports at is (nfft+1)
SegmentSize=256
    Segment size, data are divided into the segments each has a total element equal to SegmentSize
OverlapSize=128
    Number of data points that are overlaped with data in previous segments 
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

Wavedir
    Mean wave direction (Degree)
theta1
    Mean wave direction as a function of frequency (Degree)
theta2
    Principal wave direction as a function of frequency (Degree)
f
    Frequency (Hz)

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
    Ux=(np.pi/5)*(2*Eta)*(np.cosh(k*hfrombed)/np.sinh(k*h)) 
    Uy=0.2*Ux
    Wavedir,theta1,theta2,f=sm.wavedirpuv(P,Ux,Uy,fs,h,4,4,'puv1','xyz',0.7,0,0.7,0,fs/2,'auto','constant','constant','beji',1025,N,256,128,'yes')

References
----------

Beji, S. (2013). 
Improved explicit approximation of linear dispersion relationship for gravity waves. 
Coastal Engineering, 73, 11-12.

Deo, M. C., Gondane, D. S., & Sanil Kumar, V. (2002). 
Analysis of wave directional spreading using neural networks. 
Journal of waterway, port, coastal, and ocean engineering, 128(1), 30-37.

Earle, M. D., McGehee, D., & Tubman, M. (1995). 
Field Wave Gaging Program, Wave Data Analysis Standard (No. WES/IR/CERC-95-2). 
ARMY ENGINEER WATERWAYS EXPERIMENT STATION VICKSBURG MS.

Ewans, K. C. (1998). 
Observations of the directional spectrum of fetch-limited waves. 
Journal of Physical Oceanography, 28(3), 495-512.

Goda, Y. (2010). 
Random seas and design of maritime structures. 
World scientific.

Grosskopf, W., Aubrey, D., Mattie, M., & Mathiesen, M. (1983). 
Field intercomparison of nearshore directional wave sensors. 
IEEE Journal of Oceanic Engineering, 8(4), 254-271.

Herbers, T. H. C., Elgar, S., & Guza, R. T. (1999). 
Directional spreading of waves in the nearshore. 
Journal of Geophysical Research: Oceans, 104(C4), 7683-7693.

Hunt, J. N. (1979). 
Direct solution of wave dispersion equation. 
Journal of the Waterway Port Coastal and Ocean Division, 105(4), 457-459.

Vatankhah, A. R., & Aghashariatmadari, Z. (2013). 
Improved explicit approximation of linear dispersion relationship for gravity waves: A discussion. 
Coastal engineering, 78, 21-22.

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
