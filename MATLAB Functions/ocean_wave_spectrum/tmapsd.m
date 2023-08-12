function [f, Syy, Hm0, fp, Tp, Tm01, Tm02, PHI] = tmapsd(U10, F, h, fp, fs, N, CalSpectralSP, transfCalcMethod, kCalcMethod, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-08-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

tmapsd
======

.. code:: MATLAB

    [f, Syy, Hm0, fp, Tp, Tm01, Tm02, PHI] = tmapsd(U10, F, h, fp, fs, N, CalSpectralSP, transfCalcMethod, kCalcMethod, dispout)

Description
-----------

Calculate TMA spectrum (power spectral density), (Bouws et al. 1985)

Inputs
------

U10=10;
    Wind velocity at 10 meter above surface level in (m/s)
F=10000;
    Wind fetch length in (m)
h=5;
    Mean water depth in (m)
fp=0.33;
    | Peak wave frequency (fp=1/Tp) in (Hz)
    | If CalSpectralSP='yes'; then fp is calculated from U10 and F
fs=2;
    Sampling frequency that data collected at in (Hz)
N=256;
    Total number of points between 0 and fs that spectrum reports at is (N+1)
CalSpectralSP='yes';
    Define to calculate spectral shape parameters or not ('yes': calculate, 'no': use given parameters by user)
transfCalcMethod='approx';
    | Transformation function from JONSWAP into TMA calculation method 
    | 'approx': approximated method, 'tucker': Tucker (1994), 'kitaigordskii': Kitaigordskii et al. (1975) 
kCalcMethod='beji';
    | Wave number calculation method 
    | 'hunt': Hunt (1979), 'beji': Beji (2013), 'vatankhah': Vatankhah and Aghashariatmadari (2013) 
    | 'goda': Goda (2010), 'exact': calculate exact value 
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

f
    Frequency (Hz)
Syy
    Wave Energy Power Spectrum (m^2/Hz)
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
PHI
    Transformation function from JONSWAP into TMA

Examples
--------

.. code:: MATLAB

    [f,Syy,Hm0,fp,Tp,Tm01,Tm02,PHI]=tmapsd(10,10000,5,0.33,2,256,'yes','approx','beji','yes');

References
----------

Beji, S. (2013). 
Improved explicit approximation of linear dispersion relationship for gravity waves. 
Coastal Engineering, 73, 11-12.

Bouws, E.; Günther, H.; Rosenthal, W., and Vincent, C.L., (1985). 
Similarity of the wind wave spectrum in finite depth water: 1. Spectral form. 
Journal of Geophysical Research: Oceans, 90(C1), 975-986.

Goda, Y. (2010). 
Random seas and design of maritime structures. 
World scientific.

Hunt, J. N. (1979). 
Direct solution of wave dispersion equation. 
Journal of the Waterway Port Coastal and Ocean Division, 105(4), 457-459.

Kitaigordskii, S. A., Krasitskii, V. P., & Zaslavskii, M. M. (1975). 
On Phillips' theory of equilibrium range in the spectra of wind-generated gravity waves. 
Journal of Physical Oceanography, 5(3), 410-420.

Vatankhah, A. R., & Aghashariatmadari, Z. (2013). 
Improved explicit approximation of linear dispersion relationship for gravity waves: A discussion. 
Coastal engineering, 78, 21-22.

Tucker, M. J. (1994). 
Nearshore waveheight during storms. 
Coastal Engineering, 24(1-2), 111-136.

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
    case 0
        U10=10; F=10000; h=5; fp=0.33; fs=2; N=256; CalSpectralSP='yes'; transfCalcMethod='approx'; kCalcMethod='beji'; dispout='no';
    case 1
        F=10000; h=5; fp=0.33; fs=2; N=256; CalSpectralSP='yes'; transfCalcMethod='approx'; kCalcMethod='beji'; dispout='no';
    case 2
        h=5; fp=0.33; fs=2; N=256; CalSpectralSP='yes'; transfCalcMethod='approx'; kCalcMethod='beji'; dispout='no';
    case 3
        fp=0.33; fs=2; N=256; CalSpectralSP='yes'; transfCalcMethod='approx'; kCalcMethod='beji'; dispout='no';
    case 4
        fs=2; N=256; CalSpectralSP='yes'; transfCalcMethod='approx'; kCalcMethod='beji'; dispout='no';
    case 5
        N=256; CalSpectralSP='yes'; transfCalcMethod='approx'; kCalcMethod='beji'; dispout='no';
    case 6
        CalSpectralSP='yes'; transfCalcMethod='approx'; kCalcMethod='beji'; dispout='no';
    case 7
        transfCalcMethod='approx'; kCalcMethod='beji'; dispout='no';
    case 8
        kCalcMethod='beji'; dispout='no';
    case 9
        dispout='no';
end

%--------------------------------------------------------------------------
%Calculating wave number (k)

%Calculating spectral shape parameters
if strcmp(CalSpectralSP,'yes')==1
    fp=3.5*9.81/U10*(9.81*F/(U10^2))^-0.33; %Peak frequency
end

wp=2.*pi.*fp; %Wave angular frequency, w=2*pi/T=2*pi*f

%Deep water
k0=(wp.^2)./9.81; %Deep water wave number
k0h=k0.*h;

%Estimation of wave number (k) from Hunt (1979)
if strcmp(kCalcMethod,'hunt')==1
    % c2gh=(k0h+(1.0+0.6522.*k0h+0.4622.*k0h.^2+0.0864.*k0h.^4+0.0675.*k0h.^5).^(-1)).^(-1); %Calculating wave number from Hunt (1979)
    % k=wp./(sqrt(c2gh.*9.81*h));
    d1=0.6666666667; d2=0.3555555556; d3=0.1608465608; d4=0.0632098765; d5=0.0217540484; d6=0.0065407983;
    kh=sqrt(k0h.^2+(k0h)./(1+d1.*k0h+d2.*k0h.^2+d3.*k0h.^3+d4.*k0h.^4+d5.*k0h.^5+d6.*k0h.^6)); %Calculating wave number from Hunt (1979)
    kp=kh./h; %Calculating wave number from Hunt (1979)

%Estimation of wave number (k) from Beji (2013)
elseif strcmp(kCalcMethod,'beji')==1
    kh=k0h.*(1+k0h.^1.09.*exp(-(1.55+1.3.*k0h+0.216.*k0h.^2)))./sqrt(tanh(k0h)); %Calculating wave number from Beji (2013)
    kp=kh./h; %Calculating wave number from Beji (2013)
    kp(wp==0)=0;

%Estimation of wave number (k) from Vatankhah and Aghashariatmadari (2013)
elseif strcmp(kCalcMethod,'vatankhah')==1
    kh=(k0h+k0h.^2.*exp(-(3.2+k0h.^1.65)))./sqrt(tanh(k0h))+k0h.*(1-exp(-(k0h.^0.132))).^(5.0532+2.158.*k0h.^1.505); %Calculating wave number from Vatankhah and Aghashariatmadari (2013)
    kp=kh./h; %Calculating wave number from Vatankhah and Aghashariatmadari (2013)
    kp(wp==0)=0;

%Estimation of wave number (k) from Goad (2010)
elseif strcmp(kCalcMethod,'goda')==1
    kh=zeros(length(k0h(:,1)),1);
    kh(k0h>=1)=k0h(k0h>=1);
    kh(k0h<1)=(k0h(k0h<1)).^0.5;
    for i=1:3
        kh=kh-((kh-k0h.*coth(kh))./(1+k0h.*((coth(kh)).^2-1))); %Calculating wave number from Goda (2010)
    end
    kp=kh./h; %Calculating wave number from Goda (2010)
    kp(wp==0)=0;

%Calculating exact wave number (k)
elseif strcmp(kCalcMethod,'exact')==1

    %Estimation of wave number (k) from Beji (2013)
    kh=k0h.*(1+k0h.^1.09.*exp(-(1.55+1.3.*k0h+0.216.*k0h.^2)))./sqrt(tanh(k0h)); %Calculating wave number from Beji (2013)
    kini=kh./h; %Initial value for k (Wave number from Beji (2013))
    kini(wp==0)=0;

    %Calculating exact wave number (k)
    kp=zeros(length(fp(:,1)),1); %Pre-assigning array to make program run faster
    for i=1:length(fp(:,1))
        if length(h(:,1))==1
            fun = @(x)(wp(i,1)^2-(9.81*x*tanh(x*h))); %w^2=g*k*tanh(kh)
        else
            fun = @(x)(wp(i,1)^2-(9.81*x*tanh(x*h(i,1)))); %w^2=g*k*tanh(kh)
        end
        kp(i,1)=fzero(fun,kini(i,1)); %Wave number
    end

end

%--------------------------------------------------------------------------
%Calculating wave length (L), wave celereity (C)

Tp1=1./fp; %Peak wave period
Lp=2*pi/kp; %Peak wave length
Cp=9.81*Tp1/(2*pi)*tanh(kp*h); %Peak wave celerity

%--------------------------------------------------------------------------
%Calculating TMA Spectrum
fs(fs<0)=-fs;
df=fs/N; %Frequency difference between consecutive samples, df=fs/N

f(:,1)=[0:df:fs/2]; %Frequency vector with (N/2+1) elements from 0 Hz to fNy=fs/2 Hz, equally spaced at df

%Calculating spectral shape parameters
alpha=0.0078*(2*pi*(U10^2)/(9.81*Lp))^0.49; %Phillip's constant
gama=2.47*(2*pi*(U10^2)/(9.81*Lp))^0.39; %peak enhancement parameter

%Spectral width parameter
sigma_a=0.07; %Lower limit for spectral width coefficient
sigma_b=0.09; %Upper limit for spectral width coefficient
sigma=ones(length(f(:,1)),1);
sigma(f<=fp)=sigma_a; %spectral width coefficient
sigma(f>fp)=sigma_b; %spectral width coefficient

q=exp(-((f-fp).^2)./(2*sigma.^2.*fp^2));

omega=2*pi.*f.*sqrt(h/9.81);

%Transformation function from JONSWAP into TMA, approximated method
if strcmp(transfCalcMethod,'approx')==1
    PHI=ones(length(omega(:,1)),1);
    PHI(omega<=1)=omega(omega<=1).^2./2;
    PHI(omega>1 & omega<2)=1-0.5*(2-omega(omega>1 & omega<2)).^2;
    PHI(omega>=2)=1;

%Transformation function from JONSWAP into TMA, exact method (Tucker, 1994)
elseif strcmp(transfCalcMethod,'tucker')==1

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

    PHI=(tanh(k.*h)).^2./(1+(2.*k.*h)./(sinh(2.*k.*h)));
    PHI(1,1)=0;

%Transformation function from JONSWAP into TMA, exact method (Kitaigordskii et al., 1975)
elseif strcmp(transfCalcMethod,'kitaigordskii')==1
    for i=length(f(:,1)):-1:2
        fun = @(x)(x*tanh(omega(i,1)^2*x)-1); %R(wh).tanh(wh^2*R(wh))=1
        if i==length(f(:,1))
            Rwh(i,1)=fzero(fun,0);
        else    
            Rwh(i,1)=fzero(fun,Rwh(i+1,1));
        end 
    end
    Rwh(1,1)=1;
    PHI=Rwh.^-2.*(1+(2.*omega.^2.*Rwh)./(sinh(2.*omega.^2.*Rwh))).^-1;
    PHI(1,1)=0;

end

Sj=alpha*9.81^2./((2*pi).^4.*(f.^5)).*exp(-1.25.*(fp./f).^4).*gama.^q; %Calculating JONSWAP Spectrum (Hasselmann et al. 1973)

Syy=Sj.*PHI; %Calculating TMA Spectrum (Bouws et al. 1985)

%Assign zero value to fisrt element of Syy
Syy(1,1)=0;

%--------------------------------------------------------------------------
%Calculating wave properties

%Calculating spectral moments
m0=sum(Syy.*f.^0*df);
m1=sum(Syy.*f.^1*df);
m2=sum(Syy.*f.^2*df);

%Calculating wave properties
Hm0=4*sqrt(m0); %Zero-Moment wave height
Tm01=m0/m1; %mean period
Tm02=(m0/m2)^0.5; %zero crossing period

%Calculation peak period
[Syymax Syymaxloc]=max(Syy(:,1));
Tp=1/f(Syymaxloc,1); %peak period
fp=1/Tp; %peak frequency

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
    plot(f(f~=0),Syy(f~=0))

    title('TMA Power Spectral Density')
    xlabel('Frequency(Hz)')
    ylabel('Spectral Density(m^2/Hz)')

end

%--------------------------------------------------------------------------
