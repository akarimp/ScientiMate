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
