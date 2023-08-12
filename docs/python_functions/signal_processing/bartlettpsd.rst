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

scientimate.bartlettpsd
=======================

.. code:: python

    f, Sxx = scientimate.bartlettpsd(x, fs=2, SegmentSize=256, WindowName='hamming', nfft=256, OutputSmoothSize=0, dispout='no')

Description
-----------

Calculate power spectral density using Bartlett's method

Inputs
------

x
    Input data
fs=2
    Sampling frequency that data collected at in (Hz)
SegmentSize=256
    Segment size, data are divided into the segments each has a total element equal to SegmentSize
WindowName='hamming'
    | Window name, define if multiplying input data by window function or not ('none': not multiplying)
    | 'none','rectangular','triangular','welch','hanning','hamming','gaussian','blackman','nuttall','blackmanharris'
nfft=length(x)
    Total number of points between 0 and fs that spectrum reports at is (nfft+1)
OutputSmoothSize=0
    | Window size for smoothing calculated spectrum (0, 1 or 2: not smoothing, reports original periodogram)
    | if WindowName='none' and OutputSmoothSize>2, then WindowName='hamming'
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

f
    Frequency in (Hz)
Sxx
    Power spectral density using Bartlett's method (m^2/Hz)

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np
    import scipy as sp
    from scipy import signal

    x=sp.signal.detrend(0.5*np.cos(2*np.pi*0.2*np.arange(0,1024,1/2))+(-0.1+(0.1-(-0.1)))*np.random.rand(1024*2))
    f,Sxx=sm.bartlettpsd(x,2,256,'hamming',2048,0,'yes')

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
