.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2018-05-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

wavespectrum2timeseries
=======================

.. code:: MATLAB

    [Eta, t, Hm0, fp, fEta, SxxEta, a, w, Phi] = wavespectrum2timeseries(f, Sxx, fs, dispout)

Description
-----------

| Generate random water wave data from a given water wave spectrum using wave superposition
| For more options use psd2timeseries

Inputs
------

f
    Frequency (Hz)
Sxx
    Wave power spectral density (m^2s)
fs=2;
    Sampling frequency that data are collected at in (Hz)
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

Eta
    Water surface level time series in (m)
t
    Time in (s)
Hm0
    Zero moment wave height (m)
fp
    Peak wave frequency (Hz), fp=1/Tp (Tp: Peak wave period (s))
fEta
    Frequency from generated time series(Hz)
SxxEta
    Power spectral density from generated time series (m^2s)
a
    Wave amplitude for for one-sided spectrum (0<fEta<fs/2) from generated time series (m)
w
    Wave angular frequency for for one-sided spectrum (0<fEta<fs/2) from generated time series (rad/s)
Phi
    Wave random phase for for one-sided spectrum (0<fEta<fs/2) from generated time series (rad)

Examples
--------

.. code:: MATLAB

    N=2^11; %Total number of points
    fs=8; %Sampling frequency
    df=fs/N; %Frequency difference 
    f(:,1)=[0:df:fs/2]; %Frequency vector 
    f(1,1)=f(2,1)/2; %Assign temporarily non-zero value to fisrt element of f to prevent division by zero
    Sxx=0.016.*9.81.^2./((2.*pi).^4.*(f.^5)).*exp(-1.25.*(0.33./f).^4); %Calculating Spectrum 
    f(1,1)=0;
    Sxx(1,1)=0;
    [Eta,t,Hm0,fp,fEta,SxxEta,a,w,Phi]=wavespectrum2timeseries(f,Sxx,fs,'yes');

References
----------

Branlard, E. (2010).
Generation of time series from a spectrum.
Technical University Denmark. National Laboratory for Sustainable Energy.

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
