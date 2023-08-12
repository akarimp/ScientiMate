.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-04-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

scientimate.levelcrossing
=========================

.. code:: python

    UpCrossIndx, UpCrossTime, UpCrossValue, CrestIndx, CrestTime, CrestValue, TroughIndx, TroughTime, TroughValue, t = scientimate.levelcrossing(x, CrossingLevel=None, fs=1, dispout='no')

Description
-----------

Calculate crossing point for a given level by using an upward zero crossing method

Inputs
------

x
    Input time series (oscillatory) data
CrossingLevel=mean(x)
    | Level that crossing reported for
    | If x is deterended wave (mean=0), then set CrossingLevel=0
fs=1
    | Sampling frequency that data collected at in (Hz)
    | If fs is not given, then default fs is fs=1;
    | If fs=1, then index of data points represents time as well
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

UpCrossIndx
    Index of up-crossing location
UpCrossTime
    Time of up-crossing location (s)
UpCrossValue
    Value of up-crossing
CrestIndx
    Index of crest location
CrestTime
    Time of wave crest location (s)
CrestValue
    Value of wave crest
TroughIndx
    Index of trough location
TroughTime
    Time of wave trough location (s)
TroughValue
    Value of wave trough
t
    Time (s)

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np
    import scipy as sp
    from scipy import signal

    fs=2 #Sampling frequency
    duration=1024 #Duration of the data
    N=fs*duration #Total number of points
    df=fs/N #Frequency difference 
    dt=1/fs #Time difference, dt=1/fs
    t=np.linspace(0,duration-dt,N) #Time
    x=sp.signal.detrend(0.5*np.cos(2*np.pi*0.2*t)+(-0.1+(0.1-(-0.1)))*np.random.rand(N))
    UpCrossIndx,UpCrossTime,UpCrossValue,CrestIndx,CrestTime,CrestValue,TroughIndx,TroughTime,TroughValue,t=sm.levelcrossing(x,np.mean(x),fs,'yes')

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
