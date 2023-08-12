function [xSmoothed] = smoothsignal(x, WindowSize, method, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2020-01-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

smoothsignal
============

.. code:: MATLAB

    [xSmoothed] = smoothsignal(x, WindowSize, method, dispout)

Description
-----------

Smooth input data using a window function

Inputs
------

x
    Input data
WindowSize=33;
    Window size (number of adjacent elements)that is used for smoothing, should be equal or larger than 3
method='moveavg';
    | Smoothing method
    | 'moveavg': moving average
    | 'lowpass': Low-Pass filter
    | 'savgol': Savitzky-Golay filter
    | 'butter': Butterworth filter
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

xSmoothed
    Smoothed data

Examples
--------

.. code:: MATLAB

    x(:,1)=sin(linspace(0,6*pi,200))+rand(1,200);
    [xSmoothed]=smoothsignal(x,33,'moveavg','yes');

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
        WindowSize=33; method='moveavg'; dispout='no';
    case 2
        method='moveavg'; dispout='no';
    case 3
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(x)==1
    x=x';
end

%--------------------------------------------------------------------------

N=length(x(:,1)); %Length of the input data
fs=2; %Sampling frequency
fc=0.25; %Cut-off frequency (Considered as fc/(fs/2)=0.25)

%Checking size of the window
Nw=WindowSize;
Nw(Nw>N)=fix(N/8);
Nw(mod(Nw,2)==0)=Nw+1; %Window and filter length must be odd
Nw(Nw<=3)=3; %Window and filter length should be equal or larger than 3


%--------------------------------------------------------------------------
%Smoothig data

if strcmp(method,'moveavg')==1
    win = ones(Nw,1); %Or win = rectwin(Nw);
    win_norm = win/sum(win(:)); %Normalize 'boxcar' window function
    xSmoothed = conv(x, win_norm, 'same'); %Smoothing input using convolution

elseif strcmp(method,'lowpass')==1
    b = fir1(Nw,fc/(fs/2),'low',hann(Nw+1));
    impulse_response = impz(b,1,Nw,fs);
    xSmoothed = conv(x, impulse_response, 'same');

elseif strcmp(method,'savgol')==1
    poly_order=5; %Polynomial order
    xSmoothed = sgolayfilt(x,poly_order,Nw);

elseif strcmp(method,'butter')==1
    filter_order=4; %Filter order
    [b,a] = butter(filter_order,fc/(fs/2),'low'); %Design Butterworth filter
    xSmoothed = filtfilt(b,a,x);

end


%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1
    
    plot(x)
    hold on
    plot(xSmoothed)
    
    title('Soothed Data')
    xlabel('Sample')
    ylabel('Data')
    legend('Input','Smoothed')
    
end

%--------------------------------------------------------------------------
