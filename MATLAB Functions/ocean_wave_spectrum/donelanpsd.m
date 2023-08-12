function [f, Syy, Hm0, fp, Tp, Tm01, Tm02] = donelanpsd(U10, F, fp, fs, N, CalSpectralSP, dispout)
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

donelanpsd
==========

.. code:: MATLAB

    [f, Syy, Hm0, fp, Tp, Tm01, Tm02] = donelanpsd(U10, F, fp, fs, N, CalSpectralSP, dispout)

Description
-----------

Calculate Donelan spectrum (power spectral density), (Donelan et al. 1985)

Inputs
------

U10=10;
    Wind velocity at 10 meter above surface level in (m/s)
F=10000;
    Wind fetch length in (m)
fp=0.33;
    | Peak wave frequency (fp=1/Tp) in (Hz)
    | If CalSpectralSP='yes'; then fp is calculated from U10 and F
fs=2;
    Sampling frequency that data collected at in (Hz)
N=256;
    | Total number of points between 0 and fs-df, where df=fs/N
    | Spectrum is reported between 0 and fs/2 with total number of points equal to (N/2+1)
    | Total number of points between 0 and fs is (N+1)
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

    [f,Syy,Hm0,fp,Tp,Tm01,Tm02]=donelanpsd(10,10000,0.33,2,256,'yes','yes');

References
----------

Donelan, M.A.; Hamilton, J., and Hui, W., 1985. Directional spectra of wind-generated waves. 
Philosophical Transactions of the Royal Society of London A: Mathematical, 
Physical and Engineering Sciences, 315(1534), 509–562.

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
        U10=10; F=10000; fp=0.33; fs=2; N=256; CalSpectralSP='yes'; dispout='no';
    case 1
        F=10000; fp=0.33; fs=2; N=256; CalSpectralSP='yes'; dispout='no';
    case 2
        fp=0.33; fs=2; N=256; CalSpectralSP='yes'; dispout='no';
    case 3
        fs=2; N=256; CalSpectralSP='yes'; dispout='no';
    case 4
        N=256; CalSpectralSP='yes'; dispout='no';
    case 5
        CalSpectralSP='yes'; dispout='no';
    case 6
        dispout='no';
end

%--------------------------------------------------------------------------
%Calculating Donelan Spectrum
fs(fs<0)=-fs;
df=fs/N; %Frequency difference between consecutive samples, df=fs/N

f(:,1)=[0:df:fs/2]; %Frequency vector with (N/2+1) elements from 0 Hz to fNy=fs/2 Hz, equally spaced at df

%Calculating spectral shape parameters
if strcmp(CalSpectralSP,'yes')==1
    fp=3.5*9.81/U10*(9.81*F/(U10^2))^-0.33; %Peak frequency
end

Cp=9.81/(2*pi*fp); %Peak wave celerity (Cp) for deep water

alpha=0.006*(U10/Cp)^0.55; %Phillip's constant

if (U10/Cp)>0.83 & (U10/Cp)<1
    gama=1.7; %Peak enhancement parameter
elseif (U10/Cp)>=1 & (U10/Cp)<5
    gama=1.7+6*log10(U10/Cp); %Peak enhancement parameter
end

%Spectral width parameter
sigma=ones(length(f(:,1)),1);
sigma=0.08.*(1+4./(U10/Cp).^3);

q=exp(-((f-fp).^2)./(2*sigma.^2.*fp^2));

Syy=alpha*9.81^2./((2*pi).^4).*fp.^-1.*f.^-4.*exp(-1.*(f./fp).^-4).*gama.^q; %Calculating Donelan Spectrum (Donelan et al. 1985)

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

    title('Donelan Power Spectral Density')
    xlabel('Frequency(Hz)')
    ylabel('Spectral Density(m^2/Hz)')

end

%--------------------------------------------------------------------------
