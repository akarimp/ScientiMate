function [Eta, t, Hm0, fp, fEta, SxxEta, a, w, Phi] = wavespectrum2timeseries(f, Sxx, fs, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2018-05-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

wavespectrum2timeseries
=======================

.. code:: MATLAB

    [Eta, t, Hm0, fp, fEta, SxxEta, a, w, Phi] = wavespectrum2timeseries(f, Sxx, fs, dispout)

Description
-----------

| Generate random water wave data from a given water wave spectrum using wave superposition
| For more options use psd2timeseries

Inputs
------

f
    Frequency (Hz)
Sxx
    Wave power spectral density (m^2s)
fs=2;
    Sampling frequency that data are collected at in (Hz)
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

Eta
    Water surface level time series in (m)
t
    Time in (s)
Hm0
    Zero moment wave height (m)
fp
    Peak wave frequency (Hz), fp=1/Tp (Tp: Peak wave period (s))
fEta
    Frequency from generated time series(Hz)
SxxEta
    Power spectral density from generated time series (m^2s)
a
    Wave amplitude for for one-sided spectrum (0<fEta<fs/2) from generated time series (m)
w
    Wave angular frequency for for one-sided spectrum (0<fEta<fs/2) from generated time series (rad/s)
Phi
    Wave random phase for for one-sided spectrum (0<fEta<fs/2) from generated time series (rad)

Examples
--------

.. code:: MATLAB

    N=2^11; %Total number of points
    fs=8; %Sampling frequency
    df=fs/N; %Frequency difference 
    f(:,1)=[0:df:fs/2]; %Frequency vector 
    f(1,1)=f(2,1)/2; %Assign temporarily non-zero value to fisrt element of f to prevent division by zero
    Sxx=0.016.*9.81.^2./((2.*pi).^4.*(f.^5)).*exp(-1.25.*(0.33./f).^4); %Calculating Spectrum 
    f(1,1)=0;
    Sxx(1,1)=0;
    [Eta,t,Hm0,fp,fEta,SxxEta,a,w,Phi]=wavespectrum2timeseries(f,Sxx,fs,'yes');

References
----------

Branlard, E. (2010).
Generation of time series from a spectrum.
Technical University Denmark. National Laboratory for Sustainable Energy.

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
    fs=2; dispout='no';
case 3
    dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(f)==1
    f=f';
end

if isrow(Sxx)==1
    Sxx=Sxx';
end

%--------------------------------------------------------------------------
%Projecting f to f1(0:fs/2,1) and Sxx to Sxx1(0:fs/2,1)
df=f(2,1)-f(1,1); %Delta f

lenf=length(f(:,1));

if f(1,1)>df
    f11(:,1)=[f(1,1)-df:-df:df];
    f11=flipud(f11);
    lenf11=length(f11(:,1));
    f1(1:lenf11,1)=f11;
    Sxx1(1:lenf11,1)=0;
else
    lenf11=0;
end

f1(lenf11+1:lenf11+lenf,1)=f;
Sxx1(lenf11+1:lenf11+lenf,1)=Sxx;

if f(end,1)<fs/2
    f12(:,1)=[f(end,1)+df:df:fs/2];
    lenf12=length(f12(:,1));
    f1(lenf11+lenf+1:lenf11+lenf+lenf12,1)=f12;
    Sxx1(lenf11+lenf+1:lenf11+lenf+lenf12,1)=0;
else
    lenf12=0;
end

fEta=f1;
sample=length(f1(:,1)); %Number of sample for 0<f<fs/2 which is equal to N/2+1
N=2*(sample-1); %Total number of points between 0<f<fs is N+1 where N/2+1=sample is a total number of points between 0<f<fs/2
dt=1/fs; %Time difference between consecutive samples, dt=1/fs
duration=N*dt; %Total time of time series
t(:,1)=[0:dt:(N-1)*dt]; %Time from 0 to T-dt, equally spaced at dt

%--------------------------------------------------------------------------
%Calculating random time series, using wave superposition

%----------------------------------------------------------------------
%Calculating random phase
mu=0; %mean=0
sigma=1; %standard deviation=1
RandNum=sigma*randn(N/2+1,1)+mu;  % Random number with mean=0 and standard deviation=1 (normal distribution)
Phi(:,1)=2.*pi*RandNum; % Random phase

%----------------------------------------------------------------------
%Calculating wave properties
%Note: Hm0=4.*sqrt(Sxx.*df), Hrms=Hm0/sqrt(2)=4/sqrt(2).*sqrt(Sxx.*df), a=Hrms/2=4/(2*sqrt(2)).*sqrt(Sxx.*df)=sqrt(2.*Sxx.*df)

Hm01(:,1)=4.*(Sxx1.*df).^0.5; %Wave height for each deltaf from one-sided spectrum
H=Hm01./sqrt(2); %Wave height: Hrms=Hm0/sqrt(2)
w(:,1)=2.*pi.*f1; %Wave angular frequency

%----------------------------------------------------------------------
%Calculating wave amplitude for one-sided spectrum
%Note: Hm0=4.*sqrt(Sxx.*df), Hrms=Hm0/sqrt(2)=4/sqrt(2).*sqrt(Sxx.*df), a=Hrms/2=4/(2*sqrt(2)).*sqrt(Sxx.*df)=sqrt(2.*Sxx.*df)
%a(:,1)=sqrt(2.*Sxx1.*df); %Wave amplitude for 0<f<fs/2
a=H(1:N/2+1,1)./2; %Wave amplitude: a=H/2 for 0<f<fs/2

%----------------------------------------------------------------------
%Calculating random time series using wave superposition
Etaij=zeros(length(t(:,1)),N/2+1); %Pre-assigning array to make program run faster
for i=1:N/2+1
    Etaij(:,i)=(a(i,1)).*cos(-w(i,1).*t+Phi(i,1));
end

Eta=sum(Etaij,2);

%----------------------------------------------------------------------
%Calculating spectrum
SxxEta=a.^2./2./df; % Wave energy spectrum Sxx=1/2*a^2/deltaf
    
%--------------------------------------------------------------------------
%Calculating wave properties

% Calculating spectrum using Welch method
% [Sxxwelch,fwelch] = pwelch(Eta,hanning(256),[],N,fs); %Wave power spectrum and Frequency
% [SxxPG,fPG] = periodogram(Eta,hamming(length(Eta(:,1))),length(Eta(:,1)),fs); %Wave power spectrum and Frequency
% SxxPGSmooth=smooth(SxxPG,0.01,'lowess');

%Calculating spectral moments
m0=sum(SxxEta.*fEta.^0*df);
m1=sum(SxxEta.*fEta.^1*df);
m2=sum(SxxEta.*fEta.^2*df);

%Calculating wave properties
Hm0=4*sqrt(m0); %Zero-Moment wave height
Tm01=m0/m1; %mean period
Tm02=(m0/m2)^0.5; %zero crossing period

%Calculation peak period
[Sxxmax loc6]=max(SxxEta(:,1));
Tp=1/fEta(loc6,1); %peak period
fp=1/Tp; %peak frequency

%Calculating peak frequency from weighted integral (Young, 1995)
% fp=(sum(Sxx.^5.*f.^1*df))./(sum(Sxx.^5.*f.^0*df)); %Peak frequency

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1
    
    val=[Hm0 fp Tp Tm01 Tm02 m0 m1 m2];
    name={'Hm0','fp','Tp','Tm01','Tm02','m0','m1','m2'};
    for i=1:length(val)
        fprintf('%14s   %g\n',name{i},val(i));
    end
    
    %plotting
    subplot(2,1,1)
    plot(f1(f1~=0),Sxx1(f1~=0))
    hold on
    plot(fEta(fEta~=0),SxxEta(fEta~=0),'--')
    % plot(fwelch(fwelch~=0),Sxxwelch(fwelch~=0))
    
    title('Power Spectral Density')
    xlabel('Frequency(Hz)')
    ylabel('Spectral Density(m^2/Hz)')
    
    legend('Input PSD','Output PSD')
    
    subplot(2,1,2)
    plot(t,Eta)
    hold on
    
    title('Water Level')
    xlabel('Time(s)')
    ylabel('Water Level(m)')
    
end

%--------------------------------------------------------------------------
