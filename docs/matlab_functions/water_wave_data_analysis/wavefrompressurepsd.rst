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

wavefrompressurepsd
===================

.. code:: MATLAB

    [Hm0, fp, Tp, Tm01, Tm02, f, Syy, m0, ftail] = wavefrompressurepsd(P, fs, h, heightfrombed, fmaxpcorr, fminpcorr, fcL, fcH, fmaxpcorrCalcMethod, Kpafterfmaxpcorr, kCalcMethod, Rho, nfft, SegmentSize, OverlapSize, dispout)

Description
-----------

Calculate wave properties from water pressure by converting it to water surface elevation power spectral density

Inputs
------

P
    Water pressure time series data in (N/m^2)
fs
    Sampling frequency that data collected at in (Hz)
h
    Water depth in (m)
heightfrombed=0;
    Height from bed that data collected at in (m)
fmaxpcorr=fs/2;
    | Maximum frequency that a pressure attenuation factor applies up on that (Hz)
    | If fmaxpcorrCalcMethod='user', then the smaller of calculated and user defined fmaxpcorr will be chosen
fminpcorr=0;
    | Minimum frequency that is used for defining fmaxpcorr if fmaxpcorrCalcMethod='auto' (Hz)
    | fminpcorr should be smaller than fp 
    | If swell energy exists, fminpcorr should be smaller than fp of wind sea (fpsea) and larger than fp of swell (fpswell) if there swell 
fcL=0;
    Low cut-off frequency, between 0*fs to 0.5*fs (Hz)
fcH=fs/2;
    High cut-off frequency, between 0*fs to 0.5*fs (Hz)
fmaxpcorrCalcMethod='auto';
    | Define if to calculate fmaxpcorr and ftail or to use user defined
    | 'user': use user defined value for fmaxpcorr
    | 'auto': automatically define value for fmaxpcorr
Kpafterfmaxpcorr='constant';
    | Define a apressure response factor, Kp, value for frequency larger than fmaxpcorr
    | 'nochange': Kp is not changed for frequency larger than fKuvmin 
    | 'one': Kp=1 for frequency larger than fmaxpcorr 
    | 'constant': Kp for f larger than fmaxpcorr stays equal to Kp at fmaxpcorr (constant)
kCalcMethod='beji';
    | Wave number calculation method 
    | 'hunt': Hunt (1979), 'beji': Beji (2013), 'vatankhah': Vatankhah and Aghashariatmadari (2013) 
    | 'goda': Goda (2010), 'exact': calculate exact value 
Rho=1000;
    Water density (kg/m^3)
nfft=length(Eta);
    Total number of points between 0 and fs that spectrum reports at is (nfft+1)
SegmentSize=256;
    Segment size, data are divided into the segments each has a total element equal to SegmentSize
OverlapSize=128;
    Number of data points that are overlapped with data in previous segments 
dispout='no';
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
ftail
    Frequency that diagnostic tail apply after that (typically: ftail=2.5fm where fm=1/Tm01)

Examples
--------

.. code:: MATLAB

    fs=2; %Sampling frequency
    duration=1024; %Duration of the data
    N=fs*duration; %Total number of points
    df=fs/N; %Frequency difference 
    dt=1/fs; %Time difference, dt=1/fs
    t(:,1)=linspace(0,duration-dt,N); %Time
    Eta(:,1)=detrend(0.5.*cos(2*pi*0.2*t)+(-0.1+(0.1-(-0.1))).*rand(N,1));
    hfrombed=4;
    h=5;
    k=0.2;
    P=Eta.*9.81.*1000.*(cosh(k*hfrombed)/cosh(k*h));
    [Hm0,fp,Tp,Tm01,Tm02,f,Syy,m0,ftail]=wavefrompressurepsd(P,fs,5,4,0.7,0.15,0,fs/2,'auto','constant','beji',1025,N,256,128,'yes');

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
