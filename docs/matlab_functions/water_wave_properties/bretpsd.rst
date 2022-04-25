.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-08-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bretpsd
=======

.. code:: MATLAB

    [f, Syy, Hm0, fp, Tp, Tm01, Tm02] = bretpsd(Hs, fp, fs, N, dispout)

Description
-----------

Calculate Bretschneider spectrum (power spectral density), (Bretschneider, 1959), ITTC spectrum 

Inputs
------

Hs=1;
    Significant wave height (m)
fp=0.33;
    Peak wave frequency in (Hz)
fs=2;
    Sampling frequency that data collected at in (Hz)
N=256;
    Total number of points between 0 and fs that spectrum reports at is (N+1)
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

f
    Frequency (Hz)
Syy
    Wave Energy Power Spectrum (m^2/Hz)
Hm0
    Zero-Moment Wave Height (m)
fp
    Peak wave frequency (Hz)
Tp
    Peak wave period (second)
Tm01
    Wave Period from m01 (second), Mean Wave Period
Tm02
    Wave Period from m02 (second), Mean Zero Crossing Period

Examples
--------

.. code:: MATLAB

    [f,Syy,Hm0,fp,Tp,Tm01,Tm02]=bretpsd(1,0.33,2,256,'yes');

References
----------

Bretschneider, C. L. (1959). 
Wave variability and wave spectra for wind-generated gravity waves
(No. TM-118). CORPS OF ENGINEERS WASHINGTON DC BEACH EROSION BOARD.

Stansberg, C. T., Contento, G., Hong, S. W., Irani, M., Ishida, S., & Mercier, R. (2002). 
The specialist committee on waves final report and recommendations to the 23rd ITTC. 
In Proceedings of the 23rd ITTC (Vol. 2, pp. 505-551).

Zwolan, P., & Czaplewski, K. (2012). 
Sea waves models used in maritime simulators. 
Zeszyty Naukowe/Akademia Morska w Szczecinie, (32 (104) z. 2), 186-190.

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
