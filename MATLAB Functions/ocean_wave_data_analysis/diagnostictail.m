function [fOut, SyyOut, Hm0, fp, Tp, Tm01, Tm02] = diagnostictail(fIn, SyyIn, ftail, tailtype, tailpower, h, transfCalcMethod, kCalcMethod, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-03-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

diagnostictail
==============

.. code:: MATLAB

    [fOut, SyyOut, Hm0, fp, Tp, Tm01, Tm02] = diagnostictail(fIn, SyyIn, ftail, tailtype, tailpower, h, transfCalcMethod, kCalcMethod, dispout)

Description
-----------

Replace a spectrum tail with JONSWAP (Hasselmann et al.,  1973) or TMA Spectrum (Bouws et al., 1985)

Inputs
------

fIn
    Frequency (Hz), Input
SyyIn
    Power spectral density (m^2/Hz), Input
ftail
    Frequency that diagnostic tail apply after that (typically: ftail=2.5fm where fm=1/Tm01)
tailtype='jonswap';
    | Define type of the diagnostic tail to be applied 
    | 'jonswap': JONSWAP Spectrum tail, 'tma': TMA Spectrum tail
tailpower=-5;
    | Tail power that diagnostic tail apply based on that 
    | tailpower=-3 for shallow water, tailpower=-5 for deep water
h=0;
    Mean water depth in (m)
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

fOut
    Frequency (Hz), Output
SyyOut
    Power spectral density (m^2/Hz) with diagnostic tail, Output
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

Examples
--------

.. code:: MATLAB

    N=2^11; %Total number of points
    fs=8; %Sampling frequency
    df=fs/N; %Frequency difference 
    f(:,1)=[0:df:fs/2]; %Frequency vector 
    f(1,1)=f(2,1)/2; %Assign temporarily non-zero value to fisrt element of f to prevent division by zero
    Syy=0.016.*9.81.^2./((2.*pi).^4.*(f.^5)).*exp(-1.25.*(0.33./f).^4); %Calculating Spectrum 
    f(1,1)=0;
    Syy(1,1)=0;

    [fOut,SyyOut,Hm0,fp,Tp,Tm01,Tm02]=diagnostictail(f,Syy,0.5,'jonswap',-5,5,'approx','beji','yes');

    [fOut,SyyOut,Hm0,fp,Tp,Tm01,Tm02]=diagnostictail(f,Syy,0.5,'tma',-3,5,'approx','beji','yes');

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

Hasselmann, K.; Barnett, T. P.; Bouws, E.; Carlson, H.; Cartwright, D. E.; Enke, K.; Ewing, J. A.; 
Gienapp, H.; Hasselmann, D. E.; Kruseman, P.; Meerbrug, A.; Muller, P.; Olbers, D. J.; Richter, K.; 
Sell, W., and Walden, H., (1973). 
Measurements of wind-wave growth and swell decay during the Joint North Sea Wave Project (JONSWAP). 
Deutsche Hydrographische Zeitschrift A80(12), 95p.

Hunt, J. N. (1979). 
Direct solution of wave dispersion equation. 
Journal of the Waterway Port Coastal and Ocean Division, 105(4), 457-459.

Vatankhah, A. R., & Aghashariatmadari, Z. (2013). 
Improved explicit approximation of linear dispersion relationship for gravity waves: A discussion. 
Coastal engineering, 78, 21-22.

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
        tailtype='jonswap'; tailpower=-5; h=0; transfCalcMethod='approx'; kCalcMethod='beji'; dispout='no';
    case 4
        tailpower=-5; h=0; transfCalcMethod='approx'; kCalcMethod='beji'; dispout='no';
    case 5
        h=0; transfCalcMethod='approx'; kCalcMethod='beji'; dispout='no';
    case 6
        transfCalcMethod='approx'; kCalcMethod='beji'; dispout='no';
    case 7
        kCalcMethod='beji'; dispout='no';
    case 8
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(fIn)==1
    fIn=fIn';
end

if isrow(SyyIn)==1
    SyyIn=SyyIn';
end

%--------------------------------------------------------------------------

Syy=SyyIn; %Preserve input Spectrum
f=fIn; %Preserve input Spectrum
df=f(2,1)-f(1,1); %Frequency difference between consecutive samples, df=fs/N

%--------------------------------------------------------------------------
%Calculating PHI for TMA Spectrum

if strcmp(tailtype,'tma')==1

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

end

%--------------------------------------------------------------------------
%Applying diagnostic tail

%Index of ftail
Indxftail=min(find(f>=ftail));

%Applying JONSWAP Spectrum tail after ftail
if strcmp(tailtype,'jonswap')==1

    Syy(f>ftail)=Syy(Indxftail,1).*(f(f>ftail)./ftail).^tailpower; %Adding diagnostic tail
    Syy(Syy<0)=0; %Syy can not be negative

%Applying TMA Spectrum tail after ftail
elseif strcmp(tailtype,'tma')==1

    Syy(f>ftail)=Syy(Indxftail,1).*(PHI(f>ftail)./PHI(Indxftail,1)).*(f(f>ftail)./ftail).^tailpower; %Adding TMA Spectrum tail
    Syy(Syy<0)=0; %Syy can not be negative

end

SyyOut=Syy; %Assigning Syy to SyyOut
fOut=f; %Assigning f to fOut

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
[Syymax SyymaxIndx]=max(Syy(:,1));
Tp=1/f(SyymaxIndx,1); %peak period
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
    plot(f,SyyIn)
    hold on
    plot(f,SyyOut)

    title('Power Spectral Density')
    xlabel('Frequency(Hz)')
    ylabel('Spectral Density(m^2/Hz)')
    legend('Input','Output')

end

%--------------------------------------------------------------------------
