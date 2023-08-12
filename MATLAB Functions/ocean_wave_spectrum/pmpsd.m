function [f, Syy, Hm0, fp, Tp, Tm01, Tm02] = pmpsd(U195, fp, fs, N, CalSpectralSP, dispout)
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

pmpsd
=====

.. code:: MATLAB

    [f, Syy, Hm0, fp, Tp, Tm01, Tm02] = pmpsd(U195, fp, fs, N, CalSpectralSP, dispout)

Description
-----------

Calculate Pierson-Moskowitz spectrum (power spectral density), (Pierson and Moskowitz 1964) 

Inputs
------

U195=10;
    Wind velocity at 19.5 meter above surface level in (m/s)
fp=0.33;
    | Peak wave frequency (fp=1/Tp) in (Hz)
    | If CalSpectralSP='yes'; then fp is calculated from U195
fs=8;
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

    [f,Syy,Hm0,fp,Tp,Tm01,Tm02]=pmpsd(10,0.33,2,256,'yes','yes');

References
----------

Pierson, W. J., & Moskowitz, L. (1964). 
A proposed spectral form for fully developed wind seas based on the similarity theory of SA Kitaigorodskii. 
Journal of geophysical research, 69(24), 5181-5190.

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
        U195=10; fp=0.33; fs=2; N=256; CalSpectralSP='yes'; dispout='no';
    case 1
        fp=0.33; fs=2; N=256; CalSpectralSP='yes'; dispout='no';
    case 2
        fs=2; N=256; CalSpectralSP='yes'; dispout='no';
    case 3
        N=256; CalSpectralSP='yes'; dispout='no';
    case 4
        CalSpectralSP='yes'; dispout='no';
    case 5
        dispout='no';
end

%--------------------------------------------------------------------------
%Calculating Pierson-Moskowitz Spectrum
fs(fs<0)=-fs;
df=fs/N; %Frequency difference between consecutive samples, df=fs/N

f(:,1)=[0:df:fs/2]; %Frequency vector with (N/2+1) elements from 0 Hz to fNy=fs/2 Hz, equally spaced at df

%Calculating spectral shape parameters
if strcmp(CalSpectralSP,'yes')==1
    fp=0.87*9.81/(2*pi*U195); %Peak frequency
end

alpha=8.1e-3; %Phillip's constant

Syy=alpha*9.81^2./((2*pi).^4.*(f.^5)).*exp(-1.25.*(fp./f).^4); %Calculating Pierson-Moskowitz Spectrum 

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
    
    title('Pierson-Moskowitz Power Spectral Density')
    xlabel('Frequency(Hz)')
    ylabel('Spectral Density(m^2/Hz)')
    
end

%--------------------------------------------------------------------------
