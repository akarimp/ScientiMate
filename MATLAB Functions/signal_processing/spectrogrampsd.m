function [f, t, Sxx] = spectrogrampsd(x, fs, SegmentSize, OverlapSize, WindowName, nfft, outtype, OutputSmoothSize, dispout)
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

spectrogrampsd
==============

.. code:: MATLAB

    [f, t, Sxx] = spectrogrampsd(x, fs, SegmentSize, OverlapSize, WindowName, nfft, outtype, OutputSmoothSize, dispout)

Description
-----------

Calculate spectrogram following Welch's method without averaging

Inputs
------

x
    Input data
fs=2;
    Sampling frequency that data collected at in (Hz)
SegmentSize=256;
    Segment size, data are divided into the segments each has a total element equal to SegmentSize
OverlapSize=0;
    | Number of data points that are overlaped with data in previous segments 
    | OverlapSize is recomneded to be half of the SegmentSize
WindowName='hamming';
    | Window name, define if multiplying input data by window function or not ('none': not multiplying)
    | 'none','rectangular','triangular','welch','hanning','hamming','gaussian','blackman','nuttall','blackmanharris'
nfft=length(x);
    Total number of points between 0 and fs that spectrum reports at is (nfft+1)
outtype='psd';
    | Define output type
    | 'psd': power spectral density, 'db': decibel   
OutputSmoothSize=0;
    Window size for smoothing calculated spectrum (0, 1 or 2: not smoothing, reports original Welch spectrum)
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

f
    Frequency in (Hz)
Sxx
    Spectrogram (m^2/Hz) or (dB/Hz)
t
    Time at midpoint of each section (s)

Examples
--------

.. code:: MATLAB

    x(:,1)=detrend(0.5.*cos(2*pi*0.2*linspace(0,1023.5,2048))+(-0.1+(0.1-(-0.1))).*rand(1,1024*2));
    [f,t,Sxx]=spectrogrampsd(x,2,256,128,'hamming',2048,'psd',0,'yes');

    x(:,1)= chirp([0:0.001:4],50,2,150,'quadratic');
    [f,t,Sxx]=spectrogrampsd(x,1000,128,64,'hamming',128,'db',0,'yes');

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
    case 1
        SegmentSize=256; fs=2; OverlapSize=0; WindowName='hamming'; nfft=256; outtype='psd'; OutputSmoothSize=0; dispout='no';
    case 2
        fs=2; OverlapSize=0; WindowName='hamming'; nfft=256; outtype='psd'; OutputSmoothSize=0; dispout='no';
    case 3
        OverlapSize=0; WindowName='hamming'; nfft=256; outtype='psd'; OutputSmoothSize=0; dispout='no';
    case 4
        WindowName='hamming'; nfft=256; outtype='psd'; OutputSmoothSize=0; dispout='no';
    case 5
        nfft=256; outtype='psd'; OutputSmoothSize=0; dispout='no';
    case 6
        outtype='psd'; OutputSmoothSize=0; dispout='no';
    case 7
        OutputSmoothSize=0; dispout='no';
    case 8
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(x)==1
    x=x';
end

%--------------------------------------------------------------------------

%Detrending the input data
x=detrend(x);

N=length(x(:,1)); %Length of the input data

%Checking size of the window
Nseg=SegmentSize; %Segment size
Nseg(Nseg>N)=N;
Nseg(Nseg<=0)=N;

%Checking size of the overlap window
Nlap=OverlapSize; %Overlap size
Nlap(Nlap>=Nseg)=Nseg-1;
Nlap(Nlap<0)=Nseg-1;

%--------------------------------------------------------------------------
%Creating window function in time domain

if strcmp(WindowName,'none')~=1
    
    %Creating window function in time domain
    Nwseg=Nseg; %Number of elements in a window function
    nwinseg(:,1)=[0:1:Nwseg-1]; %Counter from 0 to Nw-1
    
    %Calculating window function
    WnRec=ones(Nwseg,1); %Rectangular window function (moving average filter)
    WnTri=1-abs((nwinseg-((Nwseg-1)/2))/((Nwseg-1)/2)); %Triangular or Bartlett window function
    WnWelch=1-((nwinseg-((Nwseg-1)/2))/((Nwseg-1)/2)).^2; %Welch window function
    WnHann=0.5-0.5*cos(2*pi*nwinseg/(Nwseg-1)); %Hanning window function
    WnHamm=0.54-0.46*cos(2*pi*nwinseg/(Nwseg-1)); %Hamming window function
    WnGauss=exp(-0.5*((nwinseg-((Nwseg-1)/2))/(0.4*(Nwseg-1)/2)).^2); %Gaussian window function
    WnBlack=0.42659-0.49656*cos(2*pi*nwinseg/(Nwseg-1))+0.076849*cos(4*pi*nwinseg/(Nwseg-1)); %Blackman window function
    WnNutt=0.355768-0.487396*cos(2*pi*nwinseg/(Nwseg-1))+0.144232*cos(4*pi*nwinseg/(Nwseg-1))-0.012604*cos(6*pi*nwinseg/(Nwseg-1)); %Nuttall window function
    WnBlHarr=0.35875-0.48829*cos(2*pi*nwinseg/(Nwseg-1))+0.14128*cos(4*pi*nwinseg/(Nwseg-1))-0.01168*cos(6*pi*nwinseg/(Nwseg-1)); %Blackman-Harris window function
    
    %Assigning window function
    if strcmp(WindowName,'rectangular')==1
        Wn=WnRec;
    elseif strcmp(WindowName,'triangular')==1
        Wn=WnTri;
    elseif strcmp(WindowName,'welch')==1
        Wn=WnWelch;
    elseif strcmp(WindowName,'hanning')==1
        Wn=WnHann;
    elseif strcmp(WindowName,'hamming')==1
        Wn=WnHamm;
    elseif strcmp(WindowName,'gaussian')==1
        Wn=WnGauss;
    elseif strcmp(WindowName,'blackman')==1
        Wn=WnBlack;
    elseif strcmp(WindowName,'nuttall')==1
        Wn=WnNutt;
    elseif strcmp(WindowName,'blackmanharris')==1
        Wn=WnBlHarr;
    end

    Wnseg=Wn; %Window Function
    
    WnNormseg=Wnseg./sum(Wnseg(:,1)); %Normalizing a window function in time domain
    WnSqrSumseg=sum((Wnseg(:,1)).^2); %Sum of the window function in time domain
    U=WnSqrSumseg/Nwseg; %Energy loss from multiplying data by window function in time domain is compansated by 1/U
    
end

%--------------------------------------------------------------------------
%Calculating power spectral density using Welch method

%Generating frequency vector
dfsegWithZeroPadding=fs/N; %Frequency difference between consecutive samples, df=fs/N
fsegWithZeroPadding(:,1)=[0:dfsegWithZeroPadding:fs/2]; %Frequency vector from 0 Hz to fNy=fs/2 Hz, equally spaced at df
% Note: in general f(:,1)=[0:N-1]*df from 0 Hz to (fs-df) Hz, equally spaced at df

dt=1/fs; %Time difference between consecutive samples, dt=1/fs
duration=N/fs; %total time of time series
tTotal(:,1)=[0:dt:duration-dt]; %Time from 0 to T-dt, equally spaced at dt

i=1; %Counter
StartIndx=1; %Initating the first index of the fisrt data segment
EndIndx=Nseg; %Initating last index of the fisrt data segment

%Calculating total number segments
K=0; %Total number of segments
i=1; %Counter
StartIndx=0; %Initating the first index of the fisrt data segment
EndIndx=Nseg-1; %Initating last index of the fisrt data segment
while ((EndIndx-Nlap+1)+(Nseg-1)<=N | i==1)
    if i==1
        StartIndx=(i-1)*Nseg+1; %Segment start index
        EndIndx=i*Nseg; %Segment end index
    else
        StartIndx=(EndIndx-Nlap)+1; %Segment start index
        EndIndx=StartIndx+Nseg-1; %Segment end index
    end
    i=i+1; %Counter
    K=K+1; %Total number of segments
end

%Calculating power spectral density
Sxx_Segment=zeros((N/2+1),K); %Initilizing Sxx_Segment
FFT_Segment=zeros((N/2+1),K); %Initilizing FFT_Segment
t_Segment=zeros((Nseg),K); %Initilizing t_Segment
t=zeros(K,1); %Initilizing t
for i=1:K
    
    if i==1
        StartIndx=(i-1)*Nseg+1; %Segment start index
        EndIndx=i*Nseg; %Segment end index
    else
        StartIndx=(EndIndx-Nlap)+1; %Segment start index
        EndIndx=StartIndx+Nseg-1; %Segment end index
    end
    
    xseg=x(StartIndx:EndIndx,1); %Time seris segment with length of Nseg data points
    
    %Appling window function on each segment
    if strcmp(WindowName,'none')~=1
        xseg=Wnseg.*xseg; %Multiplying time data by window function
    end
    
    %Padding zero to the end of x to make its length equal to N
    xsegWithZeroPadding(1:Nseg,1)=xseg;
    xsegWithZeroPadding(Nseg+1:N,1)=0;
    
    %Fast Fourier Transform (FFT) of the water surface elevation from f=0 Hz to f=fs Hz
    xFFT=fft(xsegWithZeroPadding);
    
    %Fast Fourier Transform (FFT) of the water surface elevation from f=0 Hz to fNy=fs/2 Hz
    xFFT=xFFT(1:length(fsegWithZeroPadding(:,1))); %(N/2+1)elements if N is even and (N+1)/2 elements if N is odd
    
    %Half of the two-sided power spectral density (psd) from f=0 Hz to fNy=fs/2 Hz
    Sxx_2Sided=(1/(Nseg*fs))*(abs(xFFT)).^2; %Calculating psd using fs in (m^2/Hz)
    % Sxx_2Sided=(dt^2/(N*dt))*(abs(xFFT)).^2; %Calculating psd using dt
    
    %One-sided power spectral density (psd) from f=0 Hz to fNy=fs/2 Hz
    Sxx_1Sided=Sxx_2Sided;
    Sxx_1Sided(2:end-1,1)=2.*Sxx_1Sided(2:end-1,1); %one-side-spectrum=2*two-sided-spectrum in (m^2/Hz)
    
    %Compensating energy loss by multiplying by 1/U
    if strcmp(WindowName,'none')~=1
        Sxx_Segment(:,i)=Sxx_1Sided./U; %Multiplying by 1/U to compensate energy loss from multiplying time data by window function
    else
        Sxx_Segment(:,i)=Sxx_1Sided; 
    end        
    
    FFT_Segment(:,i)=xFFT; %FFT of each segment

    t_Segment(:,i)=tTotal(StartIndx:EndIndx,1); %Time segment with length of Nseg data points
    t(i,1)=mean(t_Segment(:,i)); %Time at midpoint of each section
    
end

SxxAvg=mean(Sxx_Segment,2); %Averaging value for each frequency to get mean power spectral density

%--------------------------------------------------------------------------
% Converting length of Sxx from SegmentSize to nfft elements

%Checking size of the nfft
nfft(nfft>N)=N;
nfft(nfft<=0)=N;

if nfft==N
    f=fsegWithZeroPadding;
    Sxx=Sxx_Segment;
    FFT=FFT_Segment;
else
    %Generating frequency vector with nfft elements
    dfnfft=fs/nfft; %Frequency difference between consecutive samples, df=fs/N
    fnfft(:,1)=[0:dfnfft:fs/2]; %Frequency vector from 0 Hz to fNy=fs/2 Hz, equally spaced at df
    % Note: in general f(:,1)=[0:N-1]*df from 0 Hz to (fs-df) Hz, equally spaced at df        
    f=fnfft;
    
    Sxx=zeros(length(f(:,1)),K); %Initilizing Sxx with size ((nfft/2+1),K)
    FFT=zeros(length(f(:,1)),K); %Initilizing FFT with size ((nfft/2+1),K)
    for i=1:K
        Sxx(:,i)=interp1(fsegWithZeroPadding,Sxx_Segment(:,i),fnfft); %Linear interpolation of (fwin,Syywin) to (f,Syy)
        FFT(:,i)=interp1(fsegWithZeroPadding,FFT_Segment(:,i),fnfft); %Linear interpolation of (fwin,FFTwin) to (f,FFT)
    end
end

%--------------------------------------------------------------------------
%Set output type

if strcmp(outtype,'psd')==1
    Sxx=Sxx; %Power spectral density (m^2/Hz)
elseif strcmp(outtype,'db')==1
    %Sxx=10.*log10(abs(Sxx)./max(abs(Sxx(:)))); %Power spectral density, Converting to decibel, (dB/Hz)
    Sxx=20.*log10(abs(FFT)./max(abs(FFT(:)))); %Power spectral density, Converting to decibel, (dB/Hz)
end    

%--------------------------------------------------------------------------
%Smoothing a power spectral density using convolution

if OutputSmoothSize>2 %Window size at least should be 3
    
    %Checking size of the window
    OutputSmoothSize(OutputSmoothSize>N)=N;
    
    %Creating window function in time domain
    Nw=OutputSmoothSize; %Number of elements in a window function
    nwin(:,1)=[0:1:Nw-1]; %Counter from 0 to Nw-1
    
    %Calculating window function
    WnRec=ones(Nw,1); %Rectangular window function (moving average filter)
    WnTri=1-abs((nwin-((Nw-1)/2))/((Nw-1)/2)); %Triangular or Bartlett window function
    WnWelch=1-((nwin-((Nw-1)/2))/((Nw-1)/2)).^2; %Welch window function
    WnHann=0.5-0.5*cos(2*pi*nwin/(Nw-1)); %Hanning window function
    WnHamm=0.54-0.46*cos(2*pi*nwin/(Nw-1)); %Hamming window function
    WnGauss=exp(-0.5*((nwin-((Nw-1)/2))/(0.4*(Nw-1)/2)).^2); %Gaussian window function
    WnBlack=0.42659-0.49656*cos(2*pi*nwin/(Nw-1))+0.076849*cos(4*pi*nwin/(Nw-1)); %Blackman window function
    WnNutt=0.355768-0.487396*cos(2*pi*nwin/(Nw-1))+0.144232*cos(4*pi*nwin/(Nw-1))-0.012604*cos(6*pi*nwin/(Nw-1)); %Nuttall window function
    WnBlHarr=0.35875-0.48829*cos(2*pi*nwin/(Nw-1))+0.14128*cos(4*pi*nwin/(Nw-1))-0.01168*cos(6*pi*nwin/(Nw-1)); %Blackman-Harris window function
    
    %Assigning window function
    if strcmp(WindowName,'rectangular')==1
        Wn=WnRec;
    elseif strcmp(WindowName,'triangular')==1
        Wn=WnTri;
    elseif strcmp(WindowName,'welch')==1
        Wn=WnWelch;
    elseif strcmp(WindowName,'hanning')==1
        Wn=WnHann;
    elseif strcmp(WindowName,'hamming')==1
        Wn=WnHamm;
    elseif strcmp(WindowName,'gaussian')==1
        Wn=WnGauss;
    elseif strcmp(WindowName,'blackman')==1
        Wn=WnBlack;
    elseif strcmp(WindowName,'nuttall')==1
        Wn=WnNutt;
    elseif strcmp(WindowName,'blackmanharris')==1
        Wn=WnBlHarr;
    elseif strcmp(WindowName,'none')==1
        Wn=WnHamm;
    end

    WnNorm=Wn./sum(Wn(:,1)); %Normalizing a window function in time domain
    
    for i=1:K
        Sxx(:,i)=conv(Sxx(:,i),WnNorm,'same'); %Smoothing a power spectral density
    end    
end

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1
    
    if length(Sxx(1,:))>1
        %[X,Y]=meshgrid(t,f);
        %[Cont,hCont]=contourf(X,Y,Sxx);
        %set(hCont,'LineColor','none');

        %pcol=pcolor(t,f,Sxx);
        %set(pcol,'edgecolor','none')

        imagesc(t,f,Sxx);
        set(gca,'YDir','normal');
        
        if strcmp(outtype,'psd')==1
            title('Power Spectral Density (m^2/Hz)')
        elseif strcmp(outtype,'db')==1
            title('Power Spectral Density (dB/Hz)')
        end
        xlabel('Time (s)')
        ylabel('Frequency (Hz)')
        colorbar
    end
        
end

%--------------------------------------------------------------------------
