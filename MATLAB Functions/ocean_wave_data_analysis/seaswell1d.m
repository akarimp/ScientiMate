function [fseparation, Hm0, Hm0sea, Hm0swell, Tp, Tpsea, Tpswell, m0, m0sea, m0swell, fp, Tm01, Tm02] = seaswell1d(f, Syy, fsepCalcMethod, fu, fmaxswell, fpminswell, Windvel, dispout)
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

seaswell1d
==========

.. code:: MATLAB

    [fseparation, Hm0, Hm0sea, Hm0swell, Tp, Tpsea, Tpswell, m0, m0sea, m0swell, fp, Tm01, Tm02] = seaswell1d(f, Syy, fsepCalcMethod, fu, fmaxswell, fpminswell, Windvel, dispout)

Description
-----------

Partition (separate) wind sea from swell in a power spectral density using an one dimensional method

Inputs
------

f
    Frequency (Hz)
Syy
    Power spectral density (m^2/Hz)
fsepCalcMethod='hwang';
    | Wind sea swell separating calculation method 
    | 'celerity': using deep water wave celerity, 'gilhousen': Gilhousen and Hervey (2001), 
    | 'hwang': Hwang et al. (2012), 'exact': calculate exact value 
fu=0.5;
    An upper bound of a spectrum integration frequency (Hz)
fmaxswell=0.25;
    Maximum frequency that swell can have, It is about 0.2 in Gulf of Mexico
fpminswell=0.1;
    A lower bound of a spectrum (minimum frequency) that is used for Tpswell calculation
Windvel=10;
    Wind velocity (m/s), only required for Gilhousen and Hervey (2001) method
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

fseparation
    Wind sea and Swell Separation Frequency (Hz)
Hm0
    Zero-Moment Wave Height (m)
Hm0sea
    Sea Zero-Moment Wave Height (m)
Hm0swell
    Swell Zero-Moment Wave Height (m)
Tp
    Peak wave period (second)
Tpsea
    Peak Sea period (second)
Tpswell
    Peak Swell Period (second)
m0
    Zero-Moment of the power spectral density (m^2)
m0sea
    Zero-Moment of the wind sea spectrum (m^2)
m0swell
    Zero-Moment of the swell spectrum (m^2)
fp
    Peak wave frequency (Hz)
Tm01
    Wave Period from m01 (second), Mean Wave Period
Tm02
    Wave Period from m02 (second), Mean Zero Crossing Period

Examples
--------

.. code:: MATLAB

    N=2^11; %Total number of points
    fs=2; %Sampling frequency
    df=fs/N; %Frequency difference 
    f(:,1)=[0:df:fs/2]; %Frequency vector 
    f(1,1)=f(2,1)/2; %Assign temporarily non-zero value to fisrt element of f to prevent division by zero
    SyySea=0.016.*9.81.^2./((2.*pi).^4.*(f.^5)).*exp(-1.25.*(0.33./f).^4); %Calculating Spectrum 
    SyySwell=0.016.*9.81.^2./((2.*pi).^4.*(f.^5)).*exp(-1.25.*(0.15./f).^4).*0.005; %Calculating Spectrum 
    Syy=SyySea+SyySwell;
    f(1,1)=0;
    Syy(1,1)=0;
    [fseparation,Hm0,Hm0sea,Hm0swell,Tp,Tpsea,Tpswell,m0,m0sea,m0swell,fp,Tm01,Tm02]=seaswell1d(f,Syy,'exact',0.5,0.3,0,10,'yes');

References
----------

Gilhousen, D. B., & Hervey, R. (2002). 
Improved estimates of swell from moored buoys. 
In Ocean Wave Measurement and Analysis (2001) (pp. 387-393).

Hwang, P. A., Ocampo-Torres, F. J., & García-Nava, H. (2012). 
Wind sea and swell separation of 1D wave spectrum by a spectrum integration method. 
Journal of Atmospheric and Oceanic Technology, 29(1), 116-128.

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
        fsepCalcMethod='hwang'; fu=0.5; fmaxswell=0.2; fpminswell=0; Windvel=10; dispout='no';
    case 3
        fu=0.5; fmaxswell=0.2; fpminswell=0; Windvel=10; dispout='no';
    case 4
        fmaxswell=0.2; fpminswell=0; Windvel=10; dispout='no';
    case 5
        fpminswell=0; Windvel=10; dispout='no';
    case 6
        Windvel=10; dispout='no';
    case 7
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(f)==1
    f=f';
end

if isrow(Syy)==1
    Syy=Syy';
end

%--------------------------------------------------------------------------

df=f(2,1)-f(1,1); %Frequency difference between consecutive samples, df=fs/N

%--------------------------------------------------------------------------
%Calculating 1D separation frequency of wind sea from swell

%Calculating 1D separation frequency of wind sea from swell using wave celerity (phase speed)
if strcmp(fsepCalcMethod,'celerity')==1

    % fseparation is associated with a frequency where wave celerity becomes equal to wind velocity 
    fseparation=9.81/(2*pi*Windvel);  

%Calculating 1D separation frequency of wind sea from swell using Gilhousen and Hervey (2001) method
elseif strcmp(fsepCalcMethod,'gilhousen')==1
    
    fup=f(f<=fu); %Frequency up to fu
    for i=1:length(fup(:,1))
        fstar(i,1)=fup(i,1);
        m0fstar(i,1)=sum(Syy(i:length(fup(:,1)),1).*f(i:length(fup(:,1)),1).^(0)*df);
        m2fstar(i,1)=sum(Syy(i:length(fup(:,1)),1).*f(i:length(fup(:,1)),1).^(2)*df);
    end
    
    alfafstar=(8.*pi.*m2fstar)./(9.81.*sqrt(m0fstar)); %Steepness function
    
    [alfafstarmax alfafstarmaxIndx]=max(alfafstar(:,1));
    fx=fstar(alfafstarmaxIndx,1);
    fseparation=max(0.75*fx,0.9*1.25/Windvel);
    
    fseparation(fseparation>fmaxswell | fseparation==inf | fseparation==NaN | fseparation==0)=fmaxswell; %fseperation is about 0.2 in Gulf of Mexico
    
%Calculating 1D separation frequency of wind sea from swell using Hwang et al. (2012) method
elseif strcmp(fsepCalcMethod,'hwang')==1
    
    fup=f(f<=fu); %Frequency up to fu
    for i=1:length(fup(:,1))
        fstar(i,1)=fup(i,1);
        m1fstar(i,1)=sum(Syy(i:length(fup(:,1)),1).*f(i:length(fup(:,1)),1).^1*df);
        mminus1fstar(i,1)=sum(Syy(i:length(fup(:,1)),1).*f(i:length(fup(:,1)),1).^(-1)*df);
    end
    
    alfafstar=(m1fstar)./sqrt(mminus1fstar);
    
    [alfafstarmax alfafstarmaxIndx]=max(alfafstar(:,1));
    fm=fstar(alfafstarmaxIndx,1);
    fseparation=24.2084*fm^3-9.2021*fm^2+1.8906*fm-0.04286;
    
    fseparation(fseparation>fmaxswell | fseparation==inf | fseparation==NaN | fseparation==0)=fmaxswell; %fseperation is about 0.2 in Gulf of Mexico
    

%Calculating the exact separation frequency of wind sea from swell
elseif strcmp(fsepCalcMethod,'exact')==1

    %Estimating 1D separation frequency of wind sea from swell using Hwang et al. (2012) method
    fup=f(f<=fu); %Frequency up to fu
    for i=1:length(fup(:,1))
        fstar(i,1)=fup(i,1);
        m1fstar(i,1)=sum(Syy(i:length(fup(:,1)),1).*f(i:length(fup(:,1)),1).^1*df);
        mminus1fstar(i,1)=sum(Syy(i:length(fup(:,1)),1).*f(i:length(fup(:,1)),1).^(-1)*df);
    end
    
    alfafstar=(m1fstar)./sqrt(mminus1fstar);
    
    [alfafstarmax alfafstarmaxIndx]=max(alfafstar(:,1));
    fm=fstar(alfafstarmaxIndx,1);
    fseparation=24.2084*fm^3-9.2021*fm^2+1.8906*fm-0.04286;
    
    fseparation(fseparation>fmaxswell | fseparation==inf | fseparation==NaN | fseparation==0)=fmaxswell; %fseperation is about 0.2 in Gulf of Mexico

    %Calculating the exact separation frequency of wind sea from swell
    fsepIndx=length(find(f<=fseparation)); %location of fseperation
    fpsea1=(sum(Syy(fsepIndx:end,1).^5.*f(fsepIndx:end,1).^1*df))./(sum(Syy(fsepIndx:end,1).^5.*f(fsepIndx:end,1).^0*df)); %Wind sea peak frequency based on fseparation from previous step
    fpswell1=(sum(Syy(1:fsepIndx,1).^5.*f(1:fsepIndx,1).^1*df))./(sum(Syy(1:fsepIndx,1).^5.*f(1:fsepIndx,1).^0*df)); %Swell peak frequency based on fseparation from previous step
    [Syyminswell Indx1]=min(Syy(f>fpswell1 & f<fpsea1));
    Indx2=length(find(f<=fpswell1)); 
    fseparation=f(Indx1+Indx2,1);
    fseparation(fseparation>fmaxswell | fseparation==inf | fseparation==NaN | fseparation==0)=fmaxswell; %fseperation is about 0.2 in Gulf of Mexico
    fseparation(length(fseparation)==0)=fmaxswell; %fseperation is about 0.2 in Gulf of Mexico

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
%Calculating wind sea and swell wave properties

fsepIndx=length(find(f<=fseparation)); %Index of fseperation

m0sea=sum(Syy(fsepIndx+1:end,1).*f(fsepIndx+1:end,1).^0*df); %Zero-Moment of wind sea 
m0swell=sum(Syy(1:fsepIndx,1).*f(1:fsepIndx,1).^0*df); %Zero-Moment of swell

Hm0sea=4*sqrt(m0sea); %Wind sea zero-Moment wave height
Hm0swell=4*sqrt(m0swell); %Swell zero-Moment wave height

[Syymaxsea SyymaxseaIndx]=max(Syy(fsepIndx:end,:));
Tpsea=1/(f(fsepIndx+SyymaxseaIndx-1,1)); %Wind sea peak period, fsepIndx+SyymaxseaIndx-1 is the location of sea peak

fpminswellIndx=length(find(f<=fpminswell)); %Index of fpminswell
fpminswellIndx(fpminswellIndx>=fsepIndx)=1;
[Syymaxswell SyymaxswellIndx]=max(Syy(fpminswellIndx:fsepIndx,:));
Tpswell=1/(f(fpminswellIndx+SyymaxswellIndx-1,1)); %Swell peak period, fpminswellIndx+SyymaxswellIndx-1 is the location of swell peak

%Calculating peak frequency from weighted integral (Young, 1995)
fpsea=(sum(Syy(fsepIndx:end,1).^5.*f(fsepIndx:end,1).^1*df))./(sum(Syy(fsepIndx:end,1).^5.*f(fsepIndx:end,1).^0*df)); %Wind sea peak frequency
fpswell=(sum(Syy(1:fsepIndx,1).^5.*f(1:fsepIndx,1).^1*df))./(sum(Syy(1:fsepIndx,1).^5.*f(1:fsepIndx,1).^0*df));%Swell peak frequency
    
%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1

    val=[fseparation Hm0 Hm0sea Hm0swell Tp Tpsea Tpswell m0 m0sea m0swell fp Tm01 Tm02];
    name={'fseparation','Hm0','Hm0sea','Hm0swell','Tp','Tpsea','Tpswell','m0','m0sea','m0swell','fp','Tm01','Tm02'};
    for i=1:length(val)
        fprintf('%14s   %g\n',name{i},val(i));
    end

    %Plotting
    plot(f,Syy)
    hold on
    plot([fseparation fseparation],[min(Syy) max(Syy)])

    title('Power Spectral Density')
    xlabel('Frequency(Hz)')
    ylabel('Spectral Density(m^2/Hz)')
    legend('Spectrum','Separation Frequency')

end

%--------------------------------------------------------------------------
