.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2020-08-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

windspectrum2timeseries
=======================

.. code:: MATLAB

    [Eta, t] = windspectrum2timeseries(f, Sxx, fs, dispout)

Description
-----------

| Generate zero-mean wind velocity time series from a given spectrum
| Add mean of the wind velocity to generated time series to get desired time series

Inputs
------

f
    Frequency (Hz)
Sxx
    | Power spectral density in ((m/s)^2/Hz)
    | Length of Sxx and f should be odd number
fs=2*max(f);
    Sampling frequency that data collected at in (Hz)
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

U
    Wind velocity time series in (m/s)
t
    Time in (s)

Examples
--------

.. code:: MATLAB

    N=2048+1; %Total number of points
    fs=2; %Sampling frequency
    df=fs/N; %Frequency difference 
    f(:,1)=[0.002:df:fs/2]; %Frequency vector 
    Sxx=1.0^2.*(6.868*100/10)./(1+10.32.*f.*100./10).^(5/3); %Calculate Spectrum
    [U, t]=windspectrum2timeseries(f, Sxx, fs, 'yes');

References
----------

Branlard, E. (2010).
Generation of time series from a spectrum.
Technical University Denmark. National Laboratory for Sustainable Energy.

Rose, S., & Apt, J. (2012). 
Generating wind time series as a hybrid of measured and simulated data. 
Wind Energy, 15(5), 699-715.

Shinozuka, M., & Jan, C. M. (1972). 
Digital simulation of random processes and its applications. 
Journal of sound and vibration, 25(1), 111-128.

Veers, P. (1984). 
Modeling stochastic wind loads on vertical axis wind turbines. 
In 25th Structures, Structural Dynamics and Materials Conference (p. 910).

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
