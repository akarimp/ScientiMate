function [Hm0, fp, Tp, Tm01, Tm02, m0, m1, m2] = wavepropfrompsd(Syy, f, fcL, fcH, dispout)
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

wavepropfrompsd
===============

.. code:: MATLAB

    [Hm0, fp, Tp, Tm01, Tm02, m0, m1, m2] = wavepropfrompsd(Syy, f, fcL, fcH, dispout)

Description
-----------

Calculate wave properties from a power spectral density

Inputs
------

Syy
    Power spectral density (m^2/Hz)
f
    Frequency (Hz)
fcL=0;
    Low cut-off frequency, between 0*fs to 0.5*fs (Hz)
fcH=max(f);
    High cut-off frequency, between 0*fs to 0.5*fs (Hz)
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
m0
    Zero-Moment of the power spectral density (m^2)
m1
    First-Moment of the power spectral density (m^2)
m2
    Second-Moment of the power spectral density (m^2)

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
    [Hm0,fp,Tp,Tm01,Tm02,m0,m1,m2]=wavepropfrompsd(Syy,f,0,8/2,'yes');

References
----------


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
        fcL=0; fcH=max(f); dispout='no';
    case 3
        fcH=max(f); dispout='no';
    case 4
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(Syy)==1
    Syy=Syy';
end

if isrow(f)==1
    f=f';
end

%--------------------------------------------------------------------------

fcL(fcL<0)=0;
fcH(fcH>max(f))=max(f);
df=f(2,1)-f(1,1); %Frequency interval

%--------------------------------------------------------------------------
%Cut off spectrum based on fcL and fcH

if fcL>f(1,1) & fcL<max(f)
    %Syy(f<fcL)=[];
    %f(f<fcL)=[];
    Syy(f<fcL)=0;
end

if fcH>f(1,1) & fcH<max(f)
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
    plot(f,Syy)

    title('Power Spectral Density')
    xlabel('Frequency(Hz)')
    ylabel('Spectral Density(m^2/Hz)')

end

%--------------------------------------------------------------------------
