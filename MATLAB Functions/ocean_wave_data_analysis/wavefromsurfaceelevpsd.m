function [Hm0, fp, Tp, Tm01, Tm02, f, Syy, m0] = wavefromsurfaceelevpsd(Eta, fs, fcL, fcH, nfft, SegmentSize, OverlapSize, dispout)
%{
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

wavefromsurfaceelevpsd
======================

.. code:: MATLAB

    [Hm0, fp, Tp, Tm01, Tm02, f, Syy, m0] = wavefromsurfaceelevpsd(Eta, fs, fcL, fcH, nfft, SegmentSize, OverlapSize, dispout)

Description
-----------

Calculate wave properties from water surface elevation power spectral density

Inputs
------

Eta
    Water surface elevation time series data in (m)
fs
    Sampling frequency that data collected at in (Hz)
fcL=0;
    Low cut-off frequency, between 0*fs to 0.5*fs (Hz)
fcH=fs/2;
    High cut-off frequency, between 0*fs to 0.5*fs (Hz)
nfft=length(Eta);
    Total number of points between 0 and fs that spectrum reports at is (nfft+1)
SegmentSize=256;
    Segment size, data are divided into the segments each has a total element equal to SegmentSize
OverlapSize=128;
    Number of data points that are overlaped with data in previous segments 
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
    [Hm0,fp,Tp,Tm01,Tm02,f,Syy,m0]=wavefromsurfaceelevpsd(Eta,fs,0,fs/2,N,256,128,'yes');

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
%}

%--------------------------------------------------------------------------
%CODE
%--------------------------------------------------------------------------
%Assign default values

switch nargin
    case 2
        fcL=0; fcH=fs/2; nfft=length(Eta(:,1)); SegmentSize=256; OverlapSize=128; dispout='no';
    case 3
        fcH=fs/2; nfft=length(Eta(:,1)); SegmentSize=256; OverlapSize=128; dispout='no';
    case 4
        nfft=length(Eta(:,1)); SegmentSize=256; OverlapSize=128; dispout='no';
    case 5
        SegmentSize=256; OverlapSize=128; dispout='no';
    case 6
        OverlapSize=128; dispout='no';
    case 7
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(Eta)==1
    Eta=Eta';
end

%--------------------------------------------------------------------------

%Deterending input data
EtaDetrended=detrend(Eta,'linear');
fcL(fcL<0)=0;
fcH(fcH>fs/2)=fix(fs/2);
%nfft=2^(nextpow2(length(EtaDetrended)));

%--------------------------------------------------------------------------
%Calculating power spectral density

%Water surface elevation power spectral density and frequency from Welch method
%[Syy,f]=pwelch(EtaDetrended,hamming(SegmentSize),OverlapSize,nfft,fs); 
[Syy,f]=pwelch(EtaDetrended,SegmentSize,[],nfft,fs); 

df=f(2,1)-f(1,1); %Frequency interval

%--------------------------------------------------------------------------
%Cut off spectrum based on fcL and fcH

if fcL>f(1,1) & fcL<fs/2
    %Syy(f<fcL)=[];
    %f(f<fcL)=[];
    Syy(f<fcL)=0;
end

if fcH>f(1,1) & fcH<fs/2
    %Syy(f>fcH)=[];
    %f(f>fcH)=[];
    Syy(f>fcH)=0;
end

%--------------------------------------------------------------------------
%Calculating wave properties

%Calculating spectral moments
m0=sum(Syy.*f.^0*df);
m1=sum(Syy.*f.^1*df);
m2=sum(Syy.*f.^2*df);

%Calculating wave properties
Hm0=4*sqrt(m0); %Zero-Moment wave height
Tm01=m0/m1; %Mean period
Tm02=(m0/m2)^0.5; %Zero crossing period

%Calculation peak period
[Syymax SyymaxIndx]=max(Syy(:,1));
Tp=1/f(SyymaxIndx,1); %Peak period
fp=1/Tp; %Peak frequency

%Calculating peak frequency from weighted integral (Young, 1995)
% fp=(sum(Syy.^5.*f.^1*df))./(sum(Syy.^5.*f.^0*df)); %Peak frequency

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1

    val=[Hm0 fp Tp Tm01 Tm02 m0 m1 m2];
    name={'Hm0','fp','Tp','Tm01','Tm02','m0','m1','m2'};
    for i=1:length(val)
        fprintf('%14s   %g\n',name{i},val(i));
    end

    %Plotting
    loglog(f(f~=0),Syy(f~=0))

    title('Power Spectral Density')
    xlabel('Frequency(Hz)')
    ylabel('Spectral Density(m^2/Hz)')

end

%--------------------------------------------------------------------------
