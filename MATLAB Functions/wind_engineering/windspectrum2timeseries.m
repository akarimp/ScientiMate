function [U, t] = windspectrum2timeseries(f, Sxx, fs, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2020-08-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

windspectrum2timeseries
=======================

.. code:: MATLAB

    [Eta, t] = windspectrum2timeseries(f, Sxx, fs, dispout)

Description
-----------

| Generate zero-mean wind velocity time series from a given spectrum
| Add mean of the wind velocity to generated time series to get desired time series

Inputs
------

f
    Frequency (Hz)
Sxx
    | Power spectral density in ((m/s)^2/Hz)
    | Length of Sxx and f should be odd number
fs=2*max(f);
    Sampling frequency that data collected at in (Hz)
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

U
    Wind velocity time series in (m/s)
t
    Time in (s)

Examples
--------

.. code:: MATLAB

    N=2048+1; %Total number of points
    fs=2; %Sampling frequency
    df=fs/N; %Frequency difference 
    f(:,1)=[0.002:df:fs/2]; %Frequency vector 
    Sxx=1.0^2.*(6.868*100/10)./(1+10.32.*f.*100./10).^(5/3); %Calculate Spectrum
    [U, t]=windspectrum2timeseries(f, Sxx, fs, 'yes');

References
----------

Branlard, E. (2010).
Generation of time series from a spectrum.
Technical University Denmark. National Laboratory for Sustainable Energy.

Rose, S., & Apt, J. (2012). 
Generating wind time series as a hybrid of measured and simulated data. 
Wind Energy, 15(5), 699-715.

Shinozuka, M., & Jan, C. M. (1972). 
Digital simulation of random processes and its applications. 
Journal of sound and vibration, 25(1), 111-128.

Veers, P. (1984). 
Modeling stochastic wind loads on vertical axis wind turbines. 
In 25th Structures, Structural Dynamics and Materials Conference (p. 910).

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
        fs=2*max(f); dispout='no';
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
%%
%Projecte f to f1(0:fs/2,1) and Sxx to Sxx1(0:fs/2,1)
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

f_U=f1;
%sample=length(f1(:,1)); %Number of sample for 0<f<fs/2 which is equal to N/2+1
sample=length(f(:,1)); %Number of sample for 0<f<fs/2 which is equal to N/2+1
N=2*(sample-1); %Total number of points between 0<f<fs is N+1 where N/2+1=sample is a total number of points between 0<f<fs/2
dt=1/fs; %Time difference between consecutive samples, dt=1/fs
duration=N*dt; %Total time of time series
t(:,1)=[0:dt:(N-1)*dt]; %Time from 0 to T-dt, equally spaced at dt

%--------------------------------------------------------------------------
%Calculating random time series, using wave superposition

%----------------------------------------------------------------------
%Calculate random phase
mu=0; %mean=0
sigma=1; %standard deviation=1
%RandNum=sigma*randn(N/2+1,1)+mu;  %Random number with mean=0 and standard deviation=1 (normal distribution)
RandNum=rand(N/2+1,1);  %Random number (uniform distribution)
Phi(:,1)=2.*pi.*RandNum; % Random phase

%----------------------------------------------------------------------
%Calculate wave properties

dw=2*pi*df;
an=sqrt(1/2*Sxx.*df).*sin(Phi); %Veers (1984)
bn=sqrt(1/2*Sxx.*df).*cos(Phi); %Veers (1984)
cn=sqrt(Sxx.*df); %Shinozuka & Jan (1972)
%cn(1,1)=0;
w(:,1)=2.*pi.*f; %Wave angular frequency

%----------------------------------------------------------------------
%Calculating random time series using wave superposition

%Shinozuka, M., & Jan, C. M. (1972). 
%Digital simulation of random processes and its applications. 
%Journal of sound and vibration, 25(1), 111-128.
%Eq(15)

Uij=zeros(length(t(:,1)),N/2+1); %Pre-assigning array to make program run faster
for i=1:N/2+1
    %Uij(:,i)=an(i,1).*sin(w(i,1).*t)+bn(i,1)*cos(w(i,1).*t); %Veers (1984)
    Uij(:,i)=sqrt(2)*cn(i,1)*cos(w(i,1).*t+Phi(i,1)); %Shinozuka & Jan (1972)
end

U=sum(Uij,2);
%U=2*U; %Veers (1984)

%----------------------------------------------------------------------
%Calculate spectrum

%[Sxxwelch,fwelch] = pwelch(U,[],[],N,fs); %Wave power spectrum and Frequency
%trapz(f,Sxx) %Input spectrum
%trapz(fwelch,Sxxwelch) %Spectrum from generated time series

%plot(sqrt(Sxxwelch./Sxx))
%figure(2)

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1
    
    subplot(2,1,1)
    plot(f(f~=0),Sxx(f~=0))
    %hold on
    %plot(fwelch(fwelch~=0),Sxxwelch(fwelch~=0))
    
    title('Power Spectral Density')
    xlabel('Frequency(Hz)')
    ylabel('Spectral Density ((m/s)^2/Hz)')
    
    %legend('Input PSD','Output PSD')
    
    subplot(2,1,2)
    plot(t,U)
    
    %title('Wind Velocity')
    xlabel('Time(s)')
    ylabel('Wind Velocity (m/s)')
    
end

%--------------------------------------------------------------------------
