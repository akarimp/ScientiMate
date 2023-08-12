function [Wavedir, theta1, theta2, f] = wavedirpuv(P, Ux, Uy, fs, h, Pheightfrombed, UVheightfrombed, dirCalcMethod, coordinatesys, fmaxpcorr, fminpcorr, fKuvmin, fcL, fcH, fmaxpcorrCalcMethod, Kpafterfmaxpcorr, KuvafterfKuvmin, kCalcMethod, Rho, nfft, SegmentSize, OverlapSize, dispout)
%{
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

wavedirpuv
==========

.. code:: MATLAB

    [Wavedir, theta1, theta2, f] = wavedirpuv(P, Ux, Uy, fs, h, Pheightfrombed, UVheightfrombed, dirCalcMethod, coordinatesys, fmaxpcorr, fminpcorr, fKuvmin, fcL, fcH, fmaxpcorrCalcMethod, Kpafterfmaxpcorr, KuvafterfKuvmin, kCalcMethod, Rho, nfft, SegmentSize, OverlapSize, dispout)

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
Pheightfrombed=0;
    Pressure sensor height from bed that data collected at in (m)
UVheightfrombed=0;
    Velocity sensor height from bed that data collected at in (m)
dirCalcMethod='puv1';
    | Wave number calculation method 
    | 'puv1': PUV Method 1, 'puv2': PUV Method 2, 'puv3': PUV Method 3 
coordinatesys='xyz';
    | Define the coordinate system 
    | 'xyz': XYZ coordinate system, 'enu': ENU (East North Up) coordinate system 
    | If coordinatesys='enu', then x is East and y is North  
    | If coordinatesys='enu', results are reported with respect to true north  
    | In true north coordinate system, wave comes from as:
    | 0 degree: from north, 90 degree: from east, 180 degree: from south, 270 degree: from west  
fmaxpcorr=fs/2;
    | Maximum frequency that a pressure attenuation factor applies up on that (Hz)
    | If fmaxpcorrCalcMethod='user', then the smaller of calculated and user defined fmaxpcorr will be chosen
fminpcorr=0;
    | Minimum frequency that is used for defining fmaxpcorr if fmaxpcorrCalcMethod='auto' (Hz)
    | fminpcorr should be smaller than fp 
    | If swell energy exists, fminpcorr should be smaller than fp of wind sea (fpsea) and larger than fp of swell (fpswell) if there swell 
fKuvmin=fs/2;
    Frequency that a velocity conversion factor (Kuv) at that frequency is considered as a minimum limit for Kuv
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
KuvafterfKuvmin='constant';
    | Define conversion factor, Kuv, value for frequency larger than fKuvmin
    | 'nochange': Kuv is not changed for frequency larger than fKuvmin 
    | 'one': Kuv=1 for frequency larger than fKuvmin 
    | 'constant': Kuv for f larger than fKuvmin stays equal to Kuv at fKuvmin (constant)
kCalcMethod='beji';
    | Wave number calculation method 
    | 'hunt': Hunt (1979), 'beji': Beji (2013), 'vatankhah': Vatankhah and Aghashariatmadari (2013) 
    | 'goda': Goda (2010), 'exact': calculate exact value 
Rho=1000;
    Water density (kg/m^3)
nfft=length(P);
    Total number of points between 0 and fs that spectrum reports at is (nfft+1)
SegmentSize=256;
    Segment size, data are divided into the segments each has a total element equal to SegmentSize
OverlapSize=128;
    Number of data points that are overlaped with data in previous segments 
dispout='no';
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
    Ux=(pi/5).*(2.*Eta).*(cosh(k*hfrombed)/sinh(k*h)); 
    Uy=0.2.*Ux;
    [Wavedir,theta1,theta2,f]=wavedirpuv(P,Ux,Uy,fs,h,4,4,'puv1','xyz',0.7,0,0.7,0,fs/2,'auto','constant','constant','beji',1025,N,256,128,'yes');

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
%}

%--------------------------------------------------------------------------
%CODE
%--------------------------------------------------------------------------
%Assign default values

switch nargin
    case 5
        Pheightfrombed=0; UVheightfrombed=0; dirCalcMethod='puv1'; coordinatesys='xyz'; fmaxpcorr=fs/2; fminpcorr=0; fKuvmin=fs/2; fcL=0; fcH=fs/2; fmaxpcorrCalcMethod='auto'; Kpafterfmaxpcorr='constant'; KuvafterfKuvmin='constant'; kCalcMethod='beji'; Rho=1000; nfft=length(P); SegmentSize=256; OverlapSize=128; dispout='no';
    case 6
        UVheightfrombed=0; dirCalcMethod='puv1'; coordinatesys='xyz'; fmaxpcorr=fs/2; fminpcorr=0; fKuvmin=fs/2; fcL=0; fcH=fs/2; fmaxpcorrCalcMethod='auto'; Kpafterfmaxpcorr='constant'; KuvafterfKuvmin='constant'; kCalcMethod='beji'; Rho=1000; nfft=length(P); SegmentSize=256; OverlapSize=128; dispout='no';
    case 7
        dirCalcMethod='puv1'; coordinatesys='xyz'; fmaxpcorr=fs/2; fminpcorr=0; fKuvmin=fs/2; fcL=0; fcH=fs/2; fmaxpcorrCalcMethod='auto'; Kpafterfmaxpcorr='constant'; KuvafterfKuvmin='constant'; kCalcMethod='beji'; Rho=1000; nfft=length(P); SegmentSize=256; OverlapSize=128; dispout='no';
    case 8
        coordinatesys='xyz'; fmaxpcorr=fs/2; fminpcorr=0; fKuvmin=fs/2; fcL=0; fcH=fs/2; fmaxpcorrCalcMethod='auto'; Kpafterfmaxpcorr='constant'; KuvafterfKuvmin='constant'; kCalcMethod='beji'; Rho=1000; nfft=length(P); SegmentSize=256; OverlapSize=128; dispout='no';
    case 9
        fmaxpcorr=fs/2; fminpcorr=0; fKuvmin=fs/2; fcL=0; fcH=fs/2; fmaxpcorrCalcMethod='auto'; Kpafterfmaxpcorr='constant'; KuvafterfKuvmin='constant'; kCalcMethod='beji'; Rho=1000; nfft=length(P); SegmentSize=256; OverlapSize=128; dispout='no';
    case 10
        fminpcorr=0; fKuvmin=fs/2; fcL=0; fcH=fs/2; fmaxpcorrCalcMethod='auto'; Kpafterfmaxpcorr='constant'; KuvafterfKuvmin='constant'; kCalcMethod='beji'; Rho=1000; nfft=length(P); SegmentSize=256; OverlapSize=128; dispout='no';
    case 11
        fKuvmin=fs/2; fcL=0; fcH=fs/2; fmaxpcorrCalcMethod='auto'; Kpafterfmaxpcorr='constant'; KuvafterfKuvmin='constant'; kCalcMethod='beji'; Rho=1000; nfft=length(P); SegmentSize=256; OverlapSize=128; dispout='no';
    case 12
        fcL=0; fcH=fs/2; fmaxpcorrCalcMethod='auto'; Kpafterfmaxpcorr='constant'; KuvafterfKuvmin='constant'; kCalcMethod='beji'; Rho=1000; nfft=length(P); SegmentSize=256; OverlapSize=128; dispout='no';
    case 13
        fcH=fs/2; fmaxpcorrCalcMethod='auto'; Kpafterfmaxpcorr='constant'; KuvafterfKuvmin='constant'; kCalcMethod='beji'; Rho=1000; nfft=length(P); SegmentSize=256; OverlapSize=128; dispout='no';
    case 14
        fmaxpcorrCalcMethod='auto'; Kpafterfmaxpcorr='constant'; KuvafterfKuvmin='constant'; kCalcMethod='beji'; Rho=1000; nfft=length(P); SegmentSize=256; OverlapSize=128; dispout='no';
    case 15
        Kpafterfmaxpcorr='constant'; KuvafterfKuvmin='constant'; kCalcMethod='beji'; Rho=1000; nfft=length(P); SegmentSize=256; OverlapSize=128; dispout='no';
    case 16
        KuvafterfKuvmin='constant'; kCalcMethod='beji'; Rho=1000; nfft=length(P); SegmentSize=256; OverlapSize=128; dispout='no';
    case 17
        kCalcMethod='beji'; Rho=1000; nfft=length(P); SegmentSize=256; OverlapSize=128; dispout='no';
    case 18
        Rho=1000; nfft=length(P); SegmentSize=256; OverlapSize=128; dispout='no';
    case 19
        nfft=length(P); SegmentSize=256; OverlapSize=128; dispout='no';
    case 20
        SegmentSize=256; OverlapSize=128; dispout='no';
    case 21
        OverlapSize=128; dispout='no';
    case 22
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(P)==1
    P=P';
end

if isrow(Ux)==1
    Ux=Ux';
end

if isrow(Uy)==1
    Uy=Uy';
end

%--------------------------------------------------------------------------

%Deterending input data
PDetrended=detrend(P,'linear');
EtaDetrended=PDetrended./(Rho*9.81); %Converting pressure data to water elevation (depth) data, P=Rho.g.h
UxDetrended=detrend(Ux,'linear');
UyDetrended=detrend(Uy,'linear');

%h(h<=0)=0.001;
fmaxpcorr(fmaxpcorr>fs/2)=fix(fs/2);
fKuvmin(fKuvmin>fs/2)=fix(fs/2);
fcL(fcL<0)=0;
fcH(fcH>fs/2)=fix(fs/2);
%nfft=2^(nextpow2(length(EtaDetrended)));

%--------------------------------------------------------------------------
%Calculating power spectral density

%Power spectral density and frequency from Welch method
[Syy,f]=pwelch(EtaDetrended,SegmentSize,[],nfft,fs); 
[Suu,f]=pwelch(UxDetrended,SegmentSize,[],nfft,fs); 
[Svv,f]=pwelch(UyDetrended,SegmentSize,[],nfft,fs); 

%Cross power spectral density and frequency from Welch method
[Syu,f]=cpsd(EtaDetrended,UxDetrended,SegmentSize,[],nfft,fs);
[Syv,f]=cpsd(EtaDetrended,UyDetrended,SegmentSize,[],nfft,fs);
[Suv,f]=cpsd(UxDetrended,UyDetrended,SegmentSize,[],nfft,fs);

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
Kp=cosh(k.*Pheightfrombed)./cosh(k.*h);

%Calculating KmaxL based on linear wave theory
kmaxL=pi/(h-Pheightfrombed); %Wave number associated with fmaxpcorrL

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
KpminL=cosh(kmaxL*Pheightfrombed)/cosh(kmaxL*h); %Minimum Limit for Kp calculated based on linear wave theory
Kp(Kp<KpminL)=KpminL; %Check to avoid large amplification, Kp should be larger than minimum Kp calculated based on linear wave theory

%--------------------------------------------------------------------------
%Pressure attenuation correction

SyyBeforeKp=Syy; %Power spectral density before applying Kp

%Applying Kp on power spectral density
Syy= Syy./(Kp.^2); %Applies pressure response factor, Kp

%--------------------------------------------------------------------------
%Calculating the velocity response factor, Kuv, (Wiberg, 2008)
Kuv=(2*pi.*f).*cosh(k.*UVheightfrombed)./sinh(k.*h);
Kuv(1,1)=Kuv(2,1);

%Defining a value of Kuv for frequency larger than fKuvmin
if strcmp(KuvafterfKuvmin,'nochange')==1

    Kuv;

elseif strcmp(KuvafterfKuvmin,'one')==1

    Kuv(f>fKuvmin)=1; %Kuv for f larger than fKuvmin should be 1 (no correction)

    %Linear change of Kuv to 1 for f larger than fKuvmin with bandwith of 0.1 Hz
    Indx1=max(find(f<=fKuvmin-0.05));
    Indx2=max(find(f<=fKuvmin+0.05));
    Indx2(Indx2>length(f))=length(f);
    for i=Indx1:Indx2
        Kuv(i)=(Kuv(Indx2)-Kuv(Indx1))/(Indx2-Indx1)*(i-Indx1)+Kuv(Indx1);
    end

elseif strcmp(KuvafterfKuvmin,'constant')==1 %Apply minimum velocity conversion factor after fKuvmin
    %Kuv(Kuv<1)=1; %Conversion factor, Kuv, should be larger than 1
    fKuvminIndx=max(find(f(f<=fKuvmin)));
    fKuvminIndx(fKuvminIndx>length(f))=length(f);
    Kuv(f>fKuvmin)=Kuv(fKuvminIndx,1); %Kuv for f larger than fKuvmin stays equal to Kuv at fKuvmin (constant)
end

%Calculating KmaxL based on linear wave theory
kmaxL=pi/(h-UVheightfrombed); %Wave number associated with fmaxuvcorrL

%Calculating fmaxuvcorrL based on linear wave theory
fmaxuvcorrL=1/(2*pi)*sqrt(9.81*pi/(h-UVheightfrombed)*tanh(pi*h/(h-UVheightfrombed))); %fmaxuvcorrL based on linear wave theory

%Comparing Kuv with KuvminL calculated based on linear wave theory
KuvminL=(2*pi*fmaxuvcorrL)*cosh(kmaxL*UVheightfrombed)/sinh(kmaxL*h); %Minimum Limit for Kuv calculated based on linear wave theory
Kuv(Kuv<KuvminL)=KuvminL; %Check to avoid large amplification, Kuv should be larger than minimum Kuv calculated based on linear wave theory

%--------------------------------------------------------------------------
%Applying Kuv and Kp on power spectral density

Suu=Suu./(Kuv.^2);
Svv=Svv./(Kuv.^2);
Suv=Suv./(Kuv.^2);
Syu=Syu./(Kp.*Kuv);
Syv=Syv./(Kp.*Kuv);

%--------------------------------------------------------------------------
%Checking power spectral density for negative values

Syy(Syy<0)=0; %Syy can not be negative
Suu(Syy<0)=0; %Suu can not be negative
Svv(Syy<0)=0; %Svv can not be negative
Syu(Syy<0)=0; %Syu can not be negative
Syv(Syy<0)=0; %Syv can not be negative
Suv(Syy<0)=0; %Suv can not be negative

%--------------------------------------------------------------------------
%Cut off spectrum based on fcL and fcH

if fcL>f(1,1) & fcL<fs/2
    Syy(f<fcL)=[];
    Suu(f<fcL)=[];
    Svv(f<fcL)=[];
    Syu(f<fcL)=[];
    Syv(f<fcL)=[];
    Suv(f<fcL)=[];
    f(f<fcL)=[];
end

if fcH>f(1,1) & fcH<fs/2
    Syy(f>fcH)=[];
    Suu(f>fcH)=[];
    Svv(f>fcH)=[];
    Syu(f>fcH)=[];
    Syv(f>fcH)=[];
    Suv(f>fcH)=[];
    f(f>fcH)=[];
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
%Calculating Fourier coefficients
%(e.g. Grosskopf (1983); Earle (1995); Ewans (1998); Herbers et al. (1999); Deo et al. (2002))

if strcmp(dirCalcMethod,'puv1')==1

    %Note: power spectra are already devided by Kp and Kuv
    %Note: Kuv already multiplied by (2*pi*f)
    a0=real(Syy./(2*pi));
    a1=real(Syu./(Syy.*(Suu+Svv)).^0.5);
    b1=real(Syv./(Syy.*(Suu+Svv)).^0.5);
    a2=real((Suu-Svv)./(Suu+Svv));
    b2=real(2*Suv./(Suu+Svv));

elseif strcmp(dirCalcMethod,'puv2')==1

    %Note: Power spectra are already devided by Kp and Kuv
    %Note: Kuv already multiplied by (2*pi*f)
    a0=real(Syy./(2*pi));  
    a1=real(Syu./(pi));  
    b1=real(Syv./(pi));
    a2=real((Suu-Svv)./(pi));
    b2=real(2*Suv./(pi));

end

%Checking for NaN and Inf
if strcmp(dirCalcMethod,'puv1')==1 | strcmp(dirCalcMethod,'puv2')==1 

    a0(isnan(a0)==1 | isinf(a0)==1)=0;
    a1(isnan(a1)==1 | isinf(a1)==1)=0;
    b1(isnan(b1)==1 | isinf(b1)==1)=0;
    a2(isnan(a2)==1 | isinf(a2)==1)=0;
    b2(isnan(b2)==1 | isinf(b2)==1)=0;

end

%--------------------------------------------------------------------------
%Calculating mean wave direction

if strcmp(dirCalcMethod,'puv1')==1 | strcmp(dirCalcMethod,'puv2')==1 

    theta1=atan2(real(b1),real(a1)); %Mean wave direction as a function of frequency
    theta2=0.5*atan2(real(b2),real(a2)); %Principal wave direction as a function of frequency
    %theta_m_f=atan2(real(b1),real(a1)); %Mean wave direction as a function of frequency
    %theta_m_rad=atan2(sum(sin(theta_m_f).*(Syy).*df),sum(cos(theta_m_f).*(Syy).*df)); %Mean wave direction in radian, Note: Syy already devided by (Kp)
    theta_m_rad=atan2(sum(sin(theta1).*(Syy).*df),sum(cos(theta1).*(Syy).*df)); %Mean wave direction in radian, Note: Syy already devided by (Kp)
    Wavedir=rad2deg(theta_m_rad); %Mean wave direction in degree

elseif strcmp(dirCalcMethod,'puv3')==1 %Calculating mean wave direction using FFT

    %Calculating Fast Fourier transform
    EtaFFT=fft(EtaDetrended);
    UxFFT=fft(UxDetrended);
    UyFFT=fft(UyDetrended);

    %Generate symmetrical Kp and Kuv respect with fs/2
    KpSym(:,1)=[Kp;flipud(Kp(2:end-1,1))];
    KuvSym(:,1)=[Kuv;flipud(Kuv(2:end-1,1))];

    %Applying Kp and Kuv
    EtaFFT=EtaFFT./(KpSym);
    UxFFT=UxFFT./KuvSym;
    UyFFT=UyFFT./KuvSym;

    %Inverse Fast Fourier transform
    EtaIFFT=real(ifft(EtaFFT));	
    UxIFFT=real(ifft(UxFFT));	
    UyIFFT=real(ifft(UyFFT));

    %Choosing data related to wave motion direction
    UxIFFT=UxIFFT(EtaDetrended>0);
    UyIFFT=UyIFFT(EtaDetrended>0);
    EtaIFFT=EtaIFFT(EtaDetrended>0);

    %Calculating mean wave direction
    theta1=atan2(abs(UyFFT(1:length(f(:,1)))),abs(UxFFT(1:length(f(:,1))))); %Mean wave direction as a function of frequency
    theta2=NaN; %Principal wave direction as a function of frequency
    theta_rad=atan2(UyIFFT,UxIFFT); %Time series of wave direction in radian
    theta_m_rad=(atan2(sum(sin(theta_rad).*EtaIFFT.^1),sum(cos(theta_rad).*EtaIFFT.^1))); %Mean wave direction in radian
    Wavedir=rad2deg(theta_m_rad); %Mean wave direction in degree

end

%--------------------------------------------------------------------------
%Converting mathematical angle direction to compass direction with respect to true north

if strcmp(coordinatesys,'enu')==1
    theta1=270-theta1; %Converting mathematical direction to compass direction
    theta1(theta1<0)=theta1(theta1<0)+360;

    theta2=270-theta2; %Converting mathematical direction to compass direction
    theta2(theta2<0)=theta2(theta2<0)+360;

    Wavedir=270-Wavedir; %Converting mathematical direction to compass direction
    Wavedir(Wavedir<0)=Wavedir(Wavedir<0)+360;
    Wavedir(Wavedir>360)=Wavedir(Wavedir>360)-360;
end

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1 

    val=[Hm0 fp Tp Tm01 Tm02 m0 m1 m2];
    name={'Hm0','fp','Tp','Tm01','Tm02','m0','m1','m2'};
    for i=1:length(val)
        fprintf('%14s   %g\n',name{i},val(i));
    end

    %Plotting
    [Wavedirx,Wavediry]=pol2cart(deg2rad(Wavedir),1);
    compass(Wavedirx,Wavediry)
    %polar([deg2rad(Wavedir),deg2rad(Wavedir)],[0,1])

    title('Mean Wave Direction')

end

%--------------------------------------------------------------------------
