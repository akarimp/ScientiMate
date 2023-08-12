function [Eta, t, Hm0, fp, fEta, SxxEta, a] = psd2timeseries(f, Sxx, fs, CalcMethod, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-01-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

psd2timeseries
==============

.. code:: MATLAB

    [Eta, t, Hm0, fp, fEta, SxxEta, a] = psd2timeseries(f, Sxx, fs, CalcMethod, dispout)

Description
-----------

Generate random wave data from a given spectrum

Inputs
------

f
    Frequency (Hz)
Sxx
    Power spectral density (m^2s)
fs=8;
    Sampling frequency that data collected at in (Hz)
CalcMethod='fft';
    | Method for Calculating random time series, 
    | 'fft': using Fast Fourier Transform, 'sp': using wave superposition
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

Eta
    Water Surface Level Time Series in (m)
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
    [Eta,t,Hm0,fp,fEta,SxxEta,a]=psd2timeseries(f,Sxx,fs,'fft','yes');

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
        fs=2; CalcMethod='fft'; dispout='no';
    case 3
        CalcMethod='fft'; dispout='no';
    case 4
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
%projecting f to f1(0:fs/2,1) and Sxx to Sxx1(0:fs/2,1)
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
%duration=N/fs; %Total time of time series
duration=N*dt; %Total time of time series
%t(:,1)=linspace(dt,N/fs,N); %Time from 0 to T-dt, equally spaced at dt
%t(:,1)=[0:dt:duration-dt]; %Time from 0 to T-dt, equally spaced at dt
t(:,1)=[0:dt:(N-1)*dt]; %Time from 0 to T-dt, equally spaced at dt

%--------------------------------------------------------------------------
% converts one-sided spectrum into two-sided
Sxx2=Sxx1; 
Sxx2(2:end-1,1)=Sxx1(2:end-1,1)./2; %division by 2 converts one-sided spectrum into two-sided, 
                                    %by reducing the energy in two-sides to half to make it equal with one-sided spectrum
Sxx3(1:N/2+1,1)=Sxx2; %make Sxx symetric around fs/2
Sxx3(N/2+2:N,1)=flipud(Sxx2(2:end-1,1)); %make Sxx symetric around fs/2

%--------------------------------------------------------------------------
%%
% Calculating random time series, First method, using FFT

if strcmp(CalcMethod,'fft')==1
    
    %----------------------------------------------------------------------
    %Calculating FFT from Spectrum
    
    % Note: Hm0=4.*sqrt(Sxx.*df), Hrms=Hm0/sqrt(2)=4/sqrt(2).*sqrt(Sxx.*df), a=Hrms/2=4/(2*sqrt(2)).*sqrt(Sxx.*df)=sqrt(2.*Sxx.*df)
    % Note: total number of points=N
    
    % method 1, 1st approach
    EtaFFT1=sqrt(2.*df.*Sxx3(1:N,1)); %Wave amplitude from two-sided spectrum, Eta=Sigma(dn*cos(nwt+phi)), Sxx.df=1/(2)*(dn)^2, dn=sqrt(2.Sxx.df)
    EtaFFT2=EtaFFT1.*(N); % in Matlab Y=fft(y,N)/length(y)
    
    % method 1, 2nd approach
    % EtaFFT2=sqrt(2.*(N+1).*fs.*Sxx3(1:N+1,1)); %Wave amplitude from two-sided spectrum, Sxx=2.*(1./(fs.*N)).*abs(fft)^2
    
    %----------------------------------------------------------------------
    %Calculating wave amplitude for one-sided spectrum
    dn=sqrt(2.*df.*Sxx1); %Wave amplitude from one-sided spectrum, Eta=Sigma(dn*cos(nwt+phi)), Sxx.df=1/(2)*(dn)^2, dn=sqrt(2.Sxx.df)
    a=dn; %Wave amplitude: a=dn=H/2 for 0<f<fs/2
    
    %----------------------------------------------------------------------
    %Calculating white noise
    mu=0; %mean=0
    sigma=1; %standard deviation=1
    RandNum=sigma*randn(N/2+1,1)+mu;  %Random number with mean=0 and standard deviation=1 (normal distribution)
    %RandNum=(-1)+(1-(-1)).*rand(N/2+1,1); %Random number between -1 and 1 (uniform distribution)
    Phi(:,1)=2.*pi.*RandNum; % Random phase
    WhiteNoise1=sigma.*exp(1i*Phi); % Generating white noise
    
    WhiteNoise2=flipud(WhiteNoise1(2:end-1,1));
    %WhiteNoise=WhiteNoise1; 
    %WhiteNoise(N/2+2:N+1,1)=WhiteNoise2; %make White Noise symetric around fs/2
    
    WhiteNoise=[WhiteNoise1;WhiteNoise2]; %make White Noise symetric around fs/2

    %----------------------------------------------------------------------
    %Calculating random time series
    
    %Adding white noise to spectrum (double-sided) in frequency-domain
    EtaFFT=EtaFFT2.*WhiteNoise; % (meters)
    
    %Calculating time series, method 1
    Eta=real(ifft(EtaFFT(1:N,1),N));	% corected water surface levels time series, total number of points=N+1
    Eta=Eta(1:N,1); %change the size Eta equal to size of t
    
    %Calculating random time series, method 1, 3rd approach 
    % w(:,1)=2.*pi.*f1; % Angular frequency
    
    % for i=1:N/2
    %     Etaij(:,i)=(EtaFFT1(i,1)./sqrt(2)).*cos(w(i,1).*t+Phi(i,1)); %Hrms=Hm0/sqrt(2)
    % end
    
    % for i=N/2+1:N
    %     Etaij(:,i)=Etaij(:,N/2-(i-(N/2+1))); %make Wij symetric around fs/2
    % end
    
    % for i=1:N
    %     Eta1(i,1)=sum(Etaij(i,1:end));
    % end
    
    %----------------------------------------------------------------------
    %Calculating spectrum
    
    % First method: Sxx=2.*(1./(dt^2./duration)).*abs(fft).^2 or Sxx=2.*(1./(fs.*N)).*abs(fft)^2
    Y=fft(Eta,N); %calculating Fast Fourier transform
    YMagnitude=abs(Y); %magnitude of complex fft Y=sqrt(an^2+bn^2)
    % psd=(dt^2/duration)*(YMagnitude).^2; %calculating two-side power density spectrum, %Sxx=T/dt^2*(cn)^2
    psd=(1/(N*fs))*(YMagnitude).^2; %Calculating psd using fs
    SxxEta(1,1)=psd(1,1); %assigning one-sided power density 
    SxxEta(N/2+1,1)=psd(N/2+1,1); %assigning one-sided power density 
    SxxEta(2:(N/2+1)-1,1)=2.*psd(2:(N/2+1)-1,1); %calculating one-sided power density (2:(N/2+1)-1)
                                                 %(multiplying by 2 converts two sides spectrum into one side, 
                                                 %by double the energy in one-sided spectrum to compansete for energy from eleminated side of two-sided spectrum)
    
    %Creating Hamming window function
    Nw=ceil((length(SxxEta(:,1)))/32); %Number of elements in a window function
    Nw(Nw<3)=3; %Window size should be larger than 3
    nwin(:,1)=[0:1:Nw-1]; %Counter from 0 to Nw-1

    WnHamm=0.54-0.46*cos(2*pi*nwin/(Nw-1)); %Hamming window function
    WnNorm=WnHamm./sum(WnHamm(:,1)); %Normalizing a filter
    
    % smoothing data
    SxxEta=conv(SxxEta,WnNorm,'same'); %Smoothing a power density spectra
    Hm0_Smooth=4*sqrt(sum(SxxEta.*fEta.^0*df));
    
end

%--------------------------------------------------------------------------
%%
% Calculating random time series, Second method, using wave superposition

if strcmp(CalcMethod,'sp')==1
    %----------------------------------------------------------------------
    %Calculating random phase
    mu=0; %mean=0
    sigma=1; %standard deviation=1
    RandNum=sigma*randn(N/2+1,1)+mu;  % Random number with mean=0 and standard deviation=1 (normal distribution)
    %RandNum=(-1)+(1-(-1)).*rand(N/2,1); %Random number between -1 and 1 (uniform distribution)
    Phi(:,1)=2.*pi*RandNum; % Random phase
    
    % Phi1=Phi;
    % Phi1=flipud(Phi1);
    % Phi(N/2+1:2*N/2,1)=Phi1; %make Random phase symetric around fs/2
    
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

    %for i=1:length(t(:,1))
    %    Eta(i,1)=sum(Etaij(i,1:end));
    %end
    
    %----------------------------------------------------------------------
    %Calculating spectrum
    SxxEta=a.^2./2./df; % Wave energy spectrum Sxx=1/2*a^2/deltaf
    
end

%--------------------------------------------------------------------------
%%
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
    plot(fEta(fEta~=0),SxxEta(fEta~=0))
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
