function [Hm0, fp, Tp, Tm01, Tm02, f, Syy, m0, ftail] = wavefrompressurepsd(P, fs, h, heightfrombed, fmaxpcorr, fminpcorr, fcL, fcH, fmaxpcorrCalcMethod, Kpafterfmaxpcorr, kCalcMethod, Rho, nfft, SegmentSize, OverlapSize, dispout)
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
%}

%--------------------------------------------------------------------------
%CODE
%--------------------------------------------------------------------------
%Assign default values

switch nargin
    case 3
        heightfrombed=0; fmaxpcorr=fs/2; fminpcorr=0; fcL=0; fcH=fs/2; fmaxpcorrCalcMethod='auto'; Kpafterfmaxpcorr='constant'; kCalcMethod='beji'; Rho=1000; nfft=length(P(:,1)); SegmentSize=256; OverlapSize=128; dispout='no';
    case 4
        fmaxpcorr=fs/2; fminpcorr=0; fcL=0; fcH=fs/2; fmaxpcorrCalcMethod='auto'; Kpafterfmaxpcorr='constant'; kCalcMethod='beji'; Rho=1000; nfft=length(P(:,1)); SegmentSize=256; OverlapSize=128; dispout='no';
    case 5
        fminpcorr=0; fcL=0; fcH=fs/2; fmaxpcorrCalcMethod='auto'; Kpafterfmaxpcorr='constant'; kCalcMethod='beji'; Rho=1000; nfft=length(P(:,1)); SegmentSize=256; OverlapSize=128; dispout='no';
    case 6
        fcL=0; fcH=fs/2; fmaxpcorrCalcMethod='auto'; Kpafterfmaxpcorr='constant'; kCalcMethod='beji'; Rho=1000; nfft=length(P(:,1)); SegmentSize=256; OverlapSize=128; dispout='no';
    case 7
        fcH=fs/2; fmaxpcorrCalcMethod='auto'; Kpafterfmaxpcorr='constant'; kCalcMethod='beji'; Rho=1000; nfft=length(P(:,1)); SegmentSize=256; OverlapSize=128; dispout='no';
    case 8
        fmaxpcorrCalcMethod='auto'; Kpafterfmaxpcorr='constant'; kCalcMethod='beji'; Rho=1000; nfft=length(P(:,1)); SegmentSize=256; OverlapSize=128; dispout='no';
    case 9
        Kpafterfmaxpcorr='constant'; kCalcMethod='beji'; Rho=1000; nfft=length(P(:,1)); SegmentSize=256; OverlapSize=128; dispout='no';
    case 10
        kCalcMethod='beji'; Rho=1000; nfft=length(P(:,1)); SegmentSize=256; OverlapSize=128; dispout='no';
    case 11
        Rho=1000; nfft=length(P(:,1)); SegmentSize=256; OverlapSize=128; dispout='no';
    case 12
        nfft=length(P(:,1)); SegmentSize=256; OverlapSize=128; dispout='no';
    case 13
        SegmentSize=256; OverlapSize=128; dispout='no';
    case 14
        OverlapSize=128; dispout='no';
    case 15
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(P)==1
    P=P';
end

%--------------------------------------------------------------------------

%Deterending input data
PDetrended=detrend(P,'linear');
EtaDetrended=PDetrended./(Rho*9.81); %Converting pressure data to water elevation (depth) data, P=Rho.g.h
%h(h<=0)=0.001;
fmaxpcorr(fmaxpcorr>fs/2)=fix(fs/2);
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
%Calculating wave number (k)

w=2.*pi.*f; %Wave angular frequency, w=2*pi/T=2*pi*f

%Deep water
k0=(w.^2)./9.81; %Deep water wave number
k0h=k0.*h;

%Estimation of wave number (k) from Hunt (1979)
if strcmp(kCalcMethod,'hunt')==1
    % c2gh=(k0h+(1.0+0.6522.*k0h+0.4622.*k0h.^2+0.0864.*k0h.^4+0.0675.*k0h.^5).^(-1)).^(-1); %Calculating wave number from Hunt (1979)
    % k=w./(sqrt(c2gh.*9.81*h));
    d1=0.6666666667; d2=0.3555555556; d3=0.1608465608; d4=0.0632098765; d5=0.0217540484; d6=0.0065407983;
    kh=sqrt(k0h.^2+(k0h)./(1+d1.*k0h+d2.*k0h.^2+d3.*k0h.^3+d4.*k0h.^4+d5.*k0h.^5+d6.*k0h.^6)); %Calculating wave number from Hunt (1979)
    k=kh./h; %Calculating wave number from Hunt (1979)

%Estimation of wave number (k) from Beji (2013)
elseif strcmp(kCalcMethod,'beji')==1
    kh=k0h.*(1+k0h.^1.09.*exp(-(1.55+1.3.*k0h+0.216.*k0h.^2)))./sqrt(tanh(k0h)); %Calculating wave number from Beji (2013)
    k=kh./h; %Calculating wave number from Beji (2013)
    k(w==0)=0;

%Estimation of wave number (k) from Vatankhah and Aghashariatmadari (2013)
elseif strcmp(kCalcMethod,'vatankhah')==1
    kh=(k0h+k0h.^2.*exp(-(3.2+k0h.^1.65)))./sqrt(tanh(k0h))+k0h.*(1-exp(-(k0h.^0.132))).^(5.0532+2.158.*k0h.^1.505); %Calculating wave number from Vatankhah and Aghashariatmadari (2013)
    k=kh./h; %Calculating wave number from Vatankhah and Aghashariatmadari (2013)
    k(w==0)=0;

%Estimation of wave number (k) from Goad (2010)
elseif strcmp(kCalcMethod,'goda')==1
    kh=zeros(length(k0h(:,1)),1);
    kh(k0h>=1)=k0h(k0h>=1);
    kh(k0h<1)=(k0h(k0h<1)).^0.5;
    for i=1:3
        kh=kh-((kh-k0h.*coth(kh))./(1+k0h.*((coth(kh)).^2-1))); %Calculating wave number from Goda (2010)
    end
    k=kh./h; %Calculating wave number from Goda (2010)
    k(w==0)=0;

%Calculating exact wave number (k)
elseif strcmp(kCalcMethod,'exact')==1

    %Estimation of wave number (k) from Beji (2013)
    kh=k0h.*(1+k0h.^1.09.*exp(-(1.55+1.3.*k0h+0.216.*k0h.^2)))./sqrt(tanh(k0h)); %Calculating wave number from Beji (2013)
    kini=kh./h; %Initial value for k (Wave number from Beji (2013))
    kini(w==0)=0;

    %Calculating exact wave number (k)
    k=zeros(length(f(:,1)),1); %Pre-assigning array to make program run faster
    for i=1:length(f(:,1))
        if length(h(:,1))==1
            fun = @(x)(w(i,1)^2-(9.81*x*tanh(x*h))); %w^2=g*k*tanh(kh)
        else
            fun = @(x)(w(i,1)^2-(9.81*x*tanh(x*h(i,1)))); %w^2=g*k*tanh(kh)
        end
        k(i,1)=fzero(fun,kini(i,1)); %Wave number
    end

end

%--------------------------------------------------------------------------
%Calculating the pressure response factor, Kp
Kp=cosh(k.*heightfrombed)./cosh(k.*h);

%Calculating KmaxL based on linear wave theory
kmaxL=pi/(h-heightfrombed); %Wave number associated with fmaxpcorrL

%Automatically estimating fmaxpcorr and ftailcorrection
if strcmp(fmaxpcorrCalcMethod,'auto')==1

    fminpcorrIndx=max(find(f<=fminpcorr)); %Locating the index of fminpcorr (fmaxpcorr should be larger than fminpcorr)
    [Syymax SyymaxIndx]=max(Syy(fminpcorrIndx:end,1)); %Locating the peak frequency, fp, of original spectrum before applying Kp

    fmaxpcorrL=1/(2*pi)*sqrt(9.81*kmaxL*tanh(kmaxL*h)); %Maximum frequency that Kp can be applied, calculated from linear wave theory
    fmaxpcorrLIndx=max(find(f<=fmaxpcorrL)); %Locating the index of fmaxpcorrL
    fmaxpcorrLIndx(fmaxpcorrLIndx<fminpcorrIndx+(SyymaxIndx-1))=fminpcorrIndx+(SyymaxIndx-1); %Check if fmaxpcorrLIndx locataed after fp

    Syy1=Syy./(Kp.^2); %Applying Kp on spectrum causing an increase to infinity at the tail of spectrum
    [Syymin SyyminIndx]=min(Syy1(fminpcorrIndx+(SyymaxIndx-1):fmaxpcorrLIndx,1)); %Locating the index of minimum value for Syy between fp and fmaxpcorrL

    fmaxpcorr1=f(fminpcorrIndx+(SyymaxIndx-1)+(SyyminIndx-1),1); %Asigning the frequency of the minimum value of Syy between fp and fmaxpcorrL
    fmaxpcorr1(fmaxpcorr1>fmaxpcorrL)=fmaxpcorrL; %Check fmaxpcorr1 be smaller than fmaxpcorrL
    fmaxpcorr1(fmaxpcorr1==f(fminpcorrIndx+(SyymaxIndx-1)) & fmaxpcorrL>f(fminpcorrIndx+(SyymaxIndx-1)))=fmaxpcorrL; %If fmaxpcorrL>fp then fmaxpcorr1 should not be equal to fp
    fmaxpcorr(fmaxpcorr>fmaxpcorr1)=fmaxpcorr1;

    %ftail=f(fminpcorrIndx+(SyymaxIndx-1)+(SyyminIndx-1),1); %Asigning the frequency of the minimum value of Syy between fp and fmaxpcorrL for ftail
    %ftail(ftail>fmaxpcorrL)=fmaxpcorrL;
    %ftail(ftail>fmaxpcorr1)=fmaxpcorr1;
    ftail=fmaxpcorr;

elseif strcmp(fmaxpcorrCalcMethod,'user')==1

    ftail=fmaxpcorr;

end

%Defining a value of Kp for frequency larger than fmaxpcorr
if strcmp(Kpafterfmaxpcorr,'nochange')==1

    Kp;

elseif strcmp(Kpafterfmaxpcorr,'one')==1

    Kp(f>fmaxpcorr)=1; %Kp for f larger than fmaxpcorr should be 1 (no correction)

    %Linear increase of Kp to 1 for f larger than fmaxpcorr with bandwith of 0.1 Hz
    Indx1=max(find(f<=fmaxpcorr-0.05));
    Indx2=max(find(f<=fmaxpcorr+0.05));
    Indx2(Indx2>length(f))=length(f);
    for i=Indx1:Indx2
        Kp(i)=(Kp(Indx2)-Kp(Indx1))/(Indx2-Indx1)*(i-Indx1)+Kp(Indx1);
    end

elseif strcmp(Kpafterfmaxpcorr,'constant')==1
    fmaxpcorrIndx=max(find(f<=fmaxpcorr));
    fmaxpcorrIndx(fmaxpcorrIndx>length(f))=length(f);
    Kp(f>fmaxpcorr)=Kp(fmaxpcorrIndx); %Kp for f larger than fmaxpcorr stays equal to Kp at fmaxpcorr (constant)
end

%Comparing Kp with KpminL calculated based on linear wave theory
KpminL=cosh(kmaxL*heightfrombed)/cosh(kmaxL*h); %Minimum Limit for Kp calculated based on linear wave theory
Kp(Kp<KpminL)=KpminL; %Check to avoid large amplification, Kp should be larger than minimum Kp calculated based on linear wave theory

%--------------------------------------------------------------------------
%Pressure attenuation correction

SyyBeforeKp=Syy; %Power spectral density before applying Kp

%Applying Kp on power spectral density
Syy= Syy./(Kp.^2); %Applies pressure response factor, Kp

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
