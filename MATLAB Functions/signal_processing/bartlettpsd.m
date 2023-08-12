function [f, Sxx] = bartlettpsd(x, fs, SegmentSize, WindowName, nfft, OutputSmoothSize, dispout)
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

bartlettpsd
===========

.. code:: MATLAB

    [f, Sxx] = bartlettpsd(x, fs, SegmentSize, WindowName, nfft, OutputSmoothSize, dispout)

Description
-----------

Calculate power spectral density using Bartlett's method

Inputs
------

x
    Input data
fs=2;
    Sampling frequency that data collected at in (Hz)
SegmentSize=256
    Segment size, data are divided into the segments each has a total element equal to SegmentSize
WindowName='hamming';
    | Window name, define if multiplying input data by window function or not ('none': not multiplying)
    | 'none','rectangular','triangular','welch','hanning','hamming','gaussian','blackman','nuttall','blackmanharris'
nfft=length(x);
    Total number of points between 0 and fs that spectrum reports at is (nfft+1)
OutputSmoothSize=0;
    | Window size for smoothing calculated spectrum (0, 1 or 2: not smoothing, reports original periodogram)
    | if WindowName='none' and OutputSmoothSize>2, then WindowName='hamming'
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

f
    Frequency in (Hz)
Sxx
    Power spectral density using Bartlett's method (m^2/Hz)

Examples
--------

.. code:: MATLAB

    x(:,1)=(0.5.*cos(2*pi*0.2*linspace(0,1023.5,2048))+(-0.1+(0.1-(-0.1))).*rand(1,1024*2));
    [f,Sxx]=bartlettpsd(x,2,256,'hamming',2048,0,'yes');

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
        fs=2; SegmentSize=256; WindowName='hamming'; nfft=265; OutputSmoothSize=0; dispout='no';
    case 2
        SegmentSize=256; WindowName='hamming'; nfft=265; OutputSmoothSize=0; dispout='no';
    case 3
        WindowName='hamming'; nfft=265; OutputSmoothSize=0; dispout='no';
    case 4
        nfft=265; OutputSmoothSize=0; dispout='no';
    case 5
        OutputSmoothSize=0; dispout='no';
    case 6
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
%Calculating power spectral density using Bartlett method

K=fix(N/Nseg); %Total number of segments

%Generating frequency vector
dfsegWithZeroPadding=fs/N; %Frequency difference between consecutive samples, df=fs/N
fsegWithZeroPadding(:,1)=[0:dfsegWithZeroPadding:fs/2]; %Frequency vector from 0 Hz to fNy=fs/2 Hz, equally spaced at df
% Note: in general f(:,1)=[0:N-1]*df from 0 Hz to (fs-df) Hz, equally spaced at df

Sxx_Segment=zeros((N/2+1),K); #Initilizing Sxx_Segment
for i=1:K
    
    StartIndx=(i-1)*Nseg+1; %Segment start index
    EndIndx=i*Nseg; %Segment end index
    xseg=x(StartIndx:EndIndx,1); %Time seris segment with length of Nseg data points
    
    %Multiplying time data by window function
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
    
    Sxx_Segment(:,i)=Sxx_1Sided;
    
end

SxxAvg=mean(Sxx_Segment,2); %Averaging value for each frequency to get mean power spectral density

%Compensating energy loss by multiplying time data by window function
if strcmp(WindowName,'none')~=1
    SxxAvg=SxxAvg./U; %Compensating energy loss by multiplying by 1/U
end

%--------------------------------------------------------------------------
% Converting length of Syy from SegmentSize to nfft elements

%Checking size of the nfft
nfft(nfft>N)=N;
nfft(nfft<=0)=N;

if nfft==N
    f=fsegWithZeroPadding;
    Sxx=SxxAvg;
else
    %Generating frequency vector with nfft elements
    dfnfft=fs/nfft; %Frequency difference between consecutive samples, df=fs/N
    fnfft(:,1)=[0:dfnfft:fs/2]; %Frequency vector from 0 Hz to fNy=fs/2 Hz, equally spaced at df
    % Note: in general f(:,1)=[0:N-1]*df from 0 Hz to (fs-df) Hz, equally spaced at df
    
    f=fnfft;
    Sxx=interp1(fsegWithZeroPadding,SxxAvg,fnfft); %Linear interpolation of (fwin,Syywin) to (f,Syy)
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
    
    Sxx=conv(Sxx,WnNorm,'same'); %Smoothing a power spectral density
    
end

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1
    
    loglog(f(f~=0),Sxx(f~=0))
    
    title('Power Spectral Density')
    xlabel('Frequency (Hz)')
    ylabel('Spectral Density (m^2/Hz)')
    
end

%--------------------------------------------------------------------------
