function [Hm0, fp, Tp, Tm01, Tm02, f, Syy, m0] = wavefromvelocitypsd(Ux, Uy, fs, h, heightfrombed, fKuvmin, fcL, fcH, KuvafterfKuvmin, kCalcMethod, nfft, SegmentSize, OverlapSize, dispout)
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

wavefromvelocitypsd
===================

.. code:: MATLAB

    [Hm0, fp, Tp, Tm01, Tm02, f, Syy, m0] = wavefromvelocitypsd(Ux, Uy, fs, h, heightfrombed, fKuvmin, fcL, fcH, KuvafterfKuvmin, kCalcMethod, nfft, SegmentSize, OverlapSize, dispout)

Description
-----------

Calculate wave properties from wave orbital velocity by converting it to water surface elevation power spectral density

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
heightfrombed=0;
    Height from bed that data collected at in (m)
fKuvmin=fs/2;
    Frequency that a velocity conversion factor (Kuv) at that frequency is considered as a minimum limit for Kuv
fcL=0;
    Low cut-off frequency, between 0*fs to 0.5*fs (Hz)
fcH=fs/2;
    High cut-off frequency, between 0*fs to 0.5*fs (Hz)
KuvafterfKuvmin='constant';
    | Define conversion factor, Kuv, value for frequency larger than fKuvmin
    | 'nochange': Kuv is not changed for frequency larger than fKuvmin 
    | 'one': Kuv=1 for frequency larger than fKuvmin 
    | 'constant': Kuv for f larger than fKuvmin stays equal to Kuv at fKuvmin (constant)
kCalcMethod='beji';
    | Wave number calculation method 
    | 'hunt': Hunt (1979), 'beji': Beji (2013), 'vatankhah': Vatankhah and Aghashariatmadari (2013) 
    | 'goda': Goda (2010), 'exact': calculate exact value 
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
    Ux=(pi/5).*(2.*Eta).*(cosh(k*hfrombed)/sinh(k*h)); 
    Uy=0.2.*Ux;
    [Hm0,fp,Tp,Tm01,Tm02,f,Syy,m0]=wavefromvelocitypsd(Ux,Uy,fs,5,4,0.6,0,fs/2,'constant','beji',N,256,128,'yes');

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
%}

%--------------------------------------------------------------------------
%CODE
%--------------------------------------------------------------------------
%Assign default values

switch nargin
    case 4
        heightfrombed=0; fKuvmin=fs/2; fcL=0; fcH=fs/2; KuvafterfKuvmin='constant'; kCalcMethod='beji'; nfft=length(Ux(:,1)); SegmentSize=256; OverlapSize=128; dispout='no';
    case 5
        fKuvmin=fs/2; fcL=0; fcH=fs/2; KuvafterfKuvmin='constant'; kCalcMethod='beji'; nfft=length(Ux(:,1)); SegmentSize=256; OverlapSize=128; dispout='no';
    case 6
        fcL=0; fcH=fs/2; KuvafterfKuvmin='constant'; kCalcMethod='beji'; nfft=length(Ux(:,1)); SegmentSize=256; OverlapSize=128; dispout='no';
    case 7
        KuvafterfKuvmin='constant'; kCalcMethod='beji'; nfft=length(Ux(:,1)); SegmentSize=256; OverlapSize=128; dispout='no';
    case 8
        kCalcMethod='beji'; nfft=length(Ux(:,1)); SegmentSize=256; OverlapSize=128; dispout='no';
    case 9
        nfft=length(Ux(:,1)); SegmentSize=256; OverlapSize=128; dispout='no';
    case 10
        SegmentSize=256; OverlapSize=128; dispout='no';
    case 11
        SegmentSize=256; OverlapSize=128; dispout='no';
    case 12
        OverlapSize=128; dispout='no';
    case 13
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(Ux)==1
    Ux=Ux';
end

if isrow(Uy)==1
    Uy=Uy';
end

%--------------------------------------------------------------------------

%Deterending input data
UxDetrended=detrend(Ux,'linear');
UyDetrended=detrend(Uy,'linear');
%h(h<=0)=0.001;
fKuvmin(fKuvmin>fs/2)=fix(fs/2);
fcL(fcL<0)=0;
fcH(fcH>fs/2)=fix(fs/2);
%nfft=2^(nextpow2(length(UxDetrended)));

%--------------------------------------------------------------------------
%Calculating power spectral density

%Velocity power spectral density and frequency from Welch method
%[Suu,f] = pwelch(UxDetrended,hamming(SegmentSize),OverlapSize,nfft,fs);
%[Svv,f] = pwelch(UyDetrended,hamming(SegmentSize),OverlapSize,nfft,fs);
[Suu,f] = pwelch(UxDetrended,SegmentSize,[],nfft,fs);
[Svv,f] = pwelch(UyDetrended,SegmentSize,[],nfft,fs);

Suv=Suu+Svv; %Combined horizontal spectrum (Wiberg, 2008)

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
%Calculating the velocity response factor, Kuv, (Wiberg, 2008)
Kuv=(2*pi.*f).*cosh(k.*heightfrombed)./sinh(k.*h);
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
kmaxL=pi/(h-heightfrombed); %Wave number associated with fmaxuvcorrL

%Calculating fmaxuvcorrL based on linear wave theory
fmaxuvcorrL=1/(2*pi)*sqrt(9.81*pi/(h-heightfrombed)*tanh(pi*h/(h-heightfrombed))); %fmaxuvcorrL based on linear wave theory

%Comparing Kuv with KuvminL calculated based on linear wave theory
KuvminL=(2*pi*fmaxuvcorrL)*cosh(kmaxL*heightfrombed)/sinh(kmaxL*h); %Minimum Limit for Kuv calculated based on linear wave theory
Kuv(Kuv<KuvminL)=KuvminL; %Check to avoid large amplification, Kuv should be larger than minimum Kuv calculated based on linear wave theory

%--------------------------------------------------------------------------
%Converting velocity power spectral density to water level power spectral density

Syy=Suv./(Kuv.^2); %Calculating water level spectrum from orbital velocity spectrum

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
