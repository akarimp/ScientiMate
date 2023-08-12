function [freq] = fftfrequency(N, fs)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2020-05-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

fftfrequency
============

.. code:: MATLAB

    [freq] = fftfrequency(N, fs)

Description
-----------

Return frequencies for Fast Fourier Transform

Inputs
------

N
    Total number of data points
fs=2;
    Sampling frequency (Hz)

Outputs
-------

freq
    Frequency in (Hz)

Examples
--------

.. code:: MATLAB

    [freq] = fftfrequency(64, 4);
    [freq] = fftfrequency(33, 2);

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
        fs=2;
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

%if isrow(x)==1
%    x=x';
%end

%--------------------------------------------------------------------------
%Calculate frequencies
%https://www.mathworks.com/matlabcentral/answers/141271-what-are-the-frequencies-when-n-in-fft-x-n-is-odd

dt=1/fs; %Time interval
df=fs/N; %Frequency interval

%Frequencies
if mod(N,2)==0
    freq = [0:1:N-1].*df; %All frequencies
    freq(N/2+1:end) = freq(N/2+1:end)-fs; %Map negative frequencies
else
    freq = [0:1:N-1].*df; %All frequencies
    freq((N+1)/2+1:end) = freq((N+1)/2+1:end)-fs; %Map negative frequencies
end

%--------------------------------------------------------------------------
