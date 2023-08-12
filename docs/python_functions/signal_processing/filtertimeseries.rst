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

scientimate.filtertimeseries
============================

.. code:: python

    xFiltered, t = scientimate.filtertimeseries(x, fs, fcL, fcH, dispout)

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
fcL=0
    | Low cut-off frequency, between 0*fs to 0.5*fs (Hz)
    | Signals with frequencies f<fcL are removed
    | Must be between 0 and fs/2
fcH=fs/2
    | High cut-off frequency, between 0*fs to 0.5*fs (Hz)
    | Signals with frequencies f>fcH are removed
    | Must be between 0 and fs/2
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

xFiltered
    Filtered time-series
t
    Time

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np

    #Generate time-series from 3 waves with frequencies of 0.5, 2, 4 Hz
    fs=32 #Sampling frequency
    d=20 #Duration
    f1=0.5 #1st wave frequency
    f2=2 #2nd wave frequency
    f3=4 #3rd wave frequency
    a1=1 #1st wave amplitude
    a2=0.2 #2nd wave amplitude
    a3=0.1 #3rd wave amplitude
    dt=1/fs #Time interval
    t=np.linspace(0,d-dt,fs*d) #Time data points
    x=a1*np.sin(2*np.pi*f1*t)+a2*np.sin(2*np.pi*f2*t)+a3*np.sin(2*np.pi*f3*t) #Time-series

    #Keep only wave with f1 frequency and remove waves with f2 and f3 frequencies
    fcL=f1-0.2 #Lower cut-off frequency
    fcH=f1+0.2 #Higher cut-off frequency
    xFiltered, time = sm.filtertimeseries(x, fs, fcL, fcH, 'yes')

    #Keep waves with f2 and f3 frequencies and remove wave with f1 frequency
    fcL=f2-0.2 #Lower cut-off frequency
    fcH=f3+0.2 #Higher cut-off frequency
    xFiltered, time = sm.filtertimeseries(x, fs, fcL, fcH, 'yes')

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
