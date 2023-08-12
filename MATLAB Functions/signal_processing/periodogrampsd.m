function [f, Sxx] = periodogrampsd(x, fs, WindowName, OutputSmoothSize, dispout)
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

periodogrampsd
==============

.. code:: MATLAB

    [f, Sxx] = periodogrampsd(x, fs, WindowName, OutputSmoothSize, dispout)

Description
-----------

Calculate power spectral density using periodogram method

Inputs
------

x
    Input data
fs=2;
    Sampling frequency that data collected at in (Hz)
WindowName='none';
    | Window name, define if multiplying input data by window function or not ('none': not multiplying)
    | 'none','rectangular','triangular','welch','hanning','hamming','gaussian','blackman','nuttall','blackmanharris'
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
    Power spectral density using periodogram method (m^2/Hz)

Examples
--------

.. code:: MATLAB

    x(:,1)=(0.5.*cos(2*pi*0.2*linspace(0,1023.5,2048))+(-0.1+(0.1-(-0.1))).*rand(1,1024*2));
    [f,Sxx]=periodogrampsd(x,2,'none',0,'yes');

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
        fs=2; WindowName='none'; OutputSmoothSize=0; dispout='no';
    case 2
        WindowName='none'; OutputSmoothSize=0; dispout='no';
    case 3
        OutputSmoothSize=0; dispout='no';
    case 4
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

%--------------------------------------------------------------------------
%Creating window function in time domain

if strcmp(WindowName,'none')~=1

    %Creating window function in time domain
    NwTime=N; %Window size
    nwinTime(:,1)=[0:1:NwTime-1]; %Counter from 0 to Nw-1

    %Calculating window function
    WnRec=ones(NwTime,1); %Rectangular window function (moving average filter)
    WnTri=1-abs((nwinTime-((NwTime-1)/2))/((NwTime-1)/2)); %Triangular or Bartlett window function
    WnWelch=1-((nwinTime-((NwTime-1)/2))/((NwTime-1)/2)).^2; %Welch window function
    WnHann=0.5-0.5*cos(2*pi*nwinTime/(NwTime-1)); %Hanning window function
    WnHamm=0.54-0.46*cos(2*pi*nwinTime/(NwTime-1)); %Hamming window function
    WnGauss=exp(-0.5*((nwinTime-((NwTime-1)/2))/(0.4*(NwTime-1)/2)).^2); %Gaussian window function
    WnBlack=0.42659-0.49656*cos(2*pi*nwinTime/(NwTime-1))+0.076849*cos(4*pi*nwinTime/(NwTime-1)); %Blackman window function
    WnNutt=0.355768-0.487396*cos(2*pi*nwinTime/(NwTime-1))+0.144232*cos(4*pi*nwinTime/(NwTime-1))-0.012604*cos(6*pi*nwinTime/(NwTime-1)); %Nuttall window function
    WnBlHarr=0.35875-0.48829*cos(2*pi*nwinTime/(NwTime-1))+0.14128*cos(4*pi*nwinTime/(NwTime-1))-0.01168*cos(6*pi*nwinTime/(NwTime-1)); %Blackman-Harris window function

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

    WnTime=Wn; %Window Function
    WnNormTime=WnTime./sum(WnTime(:,1)); %Normalizing a window function in time domain
    WnSqrSumTime=sum((WnTime(:,1)).^2); %Sum of the window function in time domain
    U=WnSqrSumTime/NwTime; %Energy loss from multiplying data by window function in time domain is compansated by 1/U
    
    x=x.*WnTime; %Multiplying time data by window function
end

%--------------------------------------------------------------------------
%Calculating power spectral density

%Generating frequency vector
N=length(x(:,1)); %Length of a time series
df=fs/N; %Frequency difference between consecutive samples, df=fs/N
f(:,1)=[0:df:fs/2]; %Frequency vector from 0 Hz to fNy=fs/2 Hz, equally spaced at df
% Note: in general f(:,1)=[0:N-1]*df from 0 Hz to (fs-df) Hz, equally spaced at df

%Fast Fourier Transform (FFT) of the water surface elevation from f=0 Hz to f=fs Hz
xFFT=fft(x);

%Fast Fourier Transform (FFT) of the water surface elevation from f=0 Hz to fNy=fs/2 Hz
xFFT=xFFT(1:length(f(:,1))); %(N/2+1)elements if N is even and (N+1)/2 elements if N is odd

% if mod(N,2)==0
%     xFFT=xFFT(1:N/2+1); %If N is even
% else
%     xFFT=xFFT(1:(N+1)/2); %If N is odd
% end

%Half of the two-sided power spectral density (psd) from f=0 Hz to fNy=fs/2 Hz
Sxx_2Sided=(1/(N*fs))*(abs(xFFT)).^2; %Calculating psd using fs in (m^2/Hz)
% Sxx_2Sided=(dt^2/(N*dt))*(abs(xFFT)).^2; %Calculating psd using dt

%One-sided power spectral density (psd) from f=0 Hz to fNy=fs/2 Hz
Sxx_1Sided=Sxx_2Sided;
Sxx_1Sided(2:end-1,1)=2.*Sxx_1Sided(2:end-1,1); %one-side-spectrum=2*two-sided-spectrum in (m^2/Hz)

Sxx=Sxx_1Sided; %One-sided power spectral density (psd)

%Compensating energy loss by multiplying time data by window function
if strcmp(WindowName,'none')~=1
    Sxx=Sxx./U; %Compensating energy loss by multiplying by 1/U
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
