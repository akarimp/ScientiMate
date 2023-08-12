function [f, Syy, Hm0, fp, Tp, Tm01, Tm02] = jonswappsd(U10, F, fp, fs, N, gama, CalSpectralSP, dispout)
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

jonswappsd
==========

.. code:: MATLAB

    [f, Syy, Hm0, fp, Tp, Tm01, Tm02] = jonswappsd(U10, F, fp, fs, N, gama, CalSpectralSP, dispout)

Description
-----------

Calculate JONSWAP spectrum (power spectral density), (Hasselmann et al. 1973)

Inputs
------

U10=10;
    Wind velocity at 10 meter above surface level in (m/s)
F=10000;
    Wind fetch length in (m)
fp=0.33;
    | Peak wave frequency (fp=1/Tp) in (Hz)
    | If CalSpectralSP='yes'; then fp is calculated from U10 and F
gama=3.3;
    Peak enhancement parameter (between 1 and 7)
fs=2;
    Sampling frequency that data collected at in (Hz)
N=256;
    Total number of points between 0 and fs that spectrum reports at is (N+1)
CalSpectralSP='yes';
    Define to calculate spectral shape parameters or not ('yes': calculate, 'no': use given parameters by user)
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

Examples
--------

.. code:: MATLAB

    [f,Syy,Hm0,fp,Tp,Tm01,Tm02]=jonswappsd(10,10000,0.33,2,256,3.3,'yes','yes');

References
----------

Hasselmann, K.; Barnett, T. P.; Bouws, E.; Carlson, H.; Cartwright, D. E.; Enke, K.; Ewing, J. A.; 
Gienapp, H.; Hasselmann, D. E.; Kruseman, P.; Meerbrug, A.; Muller, P.; Olbers, D. J.; Richter, K.; 
Sell, W., and Walden, H., (1973). 
Measurements of wind-wave growth and swell decay during the Joint North Sea Wave Project (JONSWAP). 
Deutsche Hydrographische Zeitschrift A80(12), 95p.

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
        U10=10; F=10000; fp=0.33; fs=2; N=256; gama=3.3; CalSpectralSP='yes'; dispout='no';
    case 1
        F=10000; fp=0.33; fs=2; N=256; gama=3.3; CalSpectralSP='yes'; dispout='no';
    case 2
        fp=0.33; fs=2; N=256; gama=3.3; CalSpectralSP='yes'; dispout='no';
    case 3
        fs=2; N=256; gama=3.3; CalSpectralSP='yes'; dispout='no';
    case 4
        N=256; gama=3.3; CalSpectralSP='yes'; dispout='no';
    case 5
        gama=3.3; CalSpectralSP='yes'; dispout='no';
    case 6
        CalSpectralSP='yes'; dispout='no';
    case 7
        dispout='no';
end

%--------------------------------------------------------------------------
%Calculating JONSWAP Spectrum
fs(fs<0)=-fs;
df=fs/N; %Frequency difference between consecutive samples, df=fs/N

f(:,1)=[0:df:fs/2]; %Frequency vector with (N/2+1) elements from 0 Hz to fNy=fs/2 Hz, equally spaced at df

%Calculating spectral shape parameters
if strcmp(CalSpectralSP,'yes')==1
    fp=3.5*9.81/U10*(9.81*F/(U10^2))^-0.33; %Peak frequency
end

alpha=0.076*(9.81*F/(U10^2))^-0.22; %Phillip's constant

%Spectral width parameter
sigma_a=0.07; %Lower limit for spectral width coefficient
sigma_b=0.09; %Upper limit for spectral width coefficient
sigma=ones(length(f(:,1)),1);
sigma(f<=fp)=sigma_a; %spectral width coefficient
sigma(f>fp)=sigma_b;  %spectral width coefficient

q=exp(-((f-fp).^2)./(2*sigma.^2.*fp^2));

Syy=alpha*9.81^2./((2*pi).^4.*(f.^5)).*exp(-1.25.*(fp./f).^4).*gama.^q; %Calculating JONSWAP Spectrum (Hasselmann et al. 1973) 

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
    
    title('JONSWAP Power Spectral Density')
    xlabel('Frequency(Hz)')
    ylabel('Spectral Density(m^2/Hz)')
    
end

%--------------------------------------------------------------------------
