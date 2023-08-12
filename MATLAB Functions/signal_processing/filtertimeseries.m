function [xFiltered, t] = filtertimeseries(x, fs, fcL, fcH, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2020-09-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

filtertimeseries
================

.. code:: MATLAB

    [xFiltered, t] = filtertimeseries(x, fs, fcL, fcH, dispout)

Description
-----------

| Filter time-series to retain signals with frequencies of fcL <= f <= fcH
| It assumes time-series is stationary

Inputs
------

x
    | Time-series
fs
    | Sampling frequency that data collected at in (Hz)
    | fs=1/dt where dt is time interval between two successive points
fcL=0;
    | Low cut-off frequency, between 0*fs to 0.5*fs (Hz)
    | Signals with frequencies f<fcL are removed
    | Must be between 0 and fs/2
fcH=fs/2;
    | High cut-off frequency, between 0*fs to 0.5*fs (Hz)
    | Signals with frequencies f>fcH are removed
    | Must be between 0 and fs/2
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

xFiltered
    Filtered time-series
t
    Time

Examples
--------

.. code:: MATLAB

    %Generate time-series from 3 waves with frequencies of 0.5, 2, 4 Hz
    fs=32; %Sampling frequency
    d=20; %Duration
    f1=0.5; %1st wave frequency
    f2=2; %2nd wave frequency
    f3=4; %3rd wave frequency
    a1=1; %1st wave amplitude
    a2=0.2; %2nd wave amplitude
    a3=0.1; %3rd wave amplitude
    dt=1/fs; %Time interval
    t=linspace(0,d-dt,fs*d); %Time data points
    x=a1*sin(2*pi*f1.*t)+a2*sin(2*pi*f2.*t)+a3*sin(2*pi*f3.*t); %Time-series

    %Keep only wave with f1 frequency and remove waves with f2 and f3 frequencies
    fcL=f1-0.2; %Lower cut-off frequency
    fcH=f1+0.2; %Higher cut-off frequency
    [xFiltered, time] = filtertimeseries(x, fs, fcL, fcH, 'yes');

    %Keep waves with f2 and f3 frequencies and remove wave with f1 frequency
    fcL=f2-0.2; %Lower cut-off frequency
    fcH=f3+0.2; %Higher cut-off frequency
    [xFiltered, time] = filtertimeseries(x, fs, fcL, fcH, 'yes');

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
        fcL=0; fcH=fs/2; dispout='no';
    case 3
        fcH=fs/2; dispout='no';
    case 4
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(x)==1
    x=x';
end

%--------------------------------------------------------------------------

N=length(x(:,1)); %Length of the input data
xDetrended=detrend(x,'constant'); %Deterending input data
xMean=mean(x); %Mean
fcL(fcL<0)=0;
fcH(fcH>fs/2)=fix(fs/2);
dt=1/fs; %Time interval
t=[0:1:N-1].*dt; %Time data points

%--------------------------------------------------------------------------
%Filter time-series

%Windowed-sinc filter order (must be even)
if mod(N,2)==0
    nf=N;
else
    nf=N+1; %Make number of datapoints even by padding a single zero at the end of FFT
end

%Band-pass windowed-sinc filter
b = fir1(nf, [fcL/(fs/2),fcH/(fs/2)], 'bandpass', hann(nf+1)); %b is impulse response
[h,f] = freqz(b, 1, length(b)-1, 'whole', fs); %h is frequency response
%plot(f,abs(h)) %Plot frequency response vs frequency

%FFT of time-series
x_fft=fft(xDetrended,nf);

%Apply band-pass windowed-sinc filter to data
x_fft_Filtered = x_fft.*abs(h);

%Inverst FFT
x_ifft_Filtered=real(ifft(x_fft_Filtered));

%Filtered time-series
if mod(N,2)==0
    xFiltered=x_ifft_Filtered;
else
    xFiltered=x_ifft_Filtered(1:end-1,1);
end

%Add mean back to time-series
xFiltered=xFiltered+xMean;

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1
    
    plot(t,x)
    hold on
    plot(t,xFiltered)
    
    title('Filtered Data')
    xlabel('Time (s)')
    ylabel('Data')
    legend('Input','Filtered')
    
end

%--------------------------------------------------------------------------
