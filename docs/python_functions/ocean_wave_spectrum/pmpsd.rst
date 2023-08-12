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

scientimate.pmpsd
=================

.. code:: python

    f, Syy, Hm0, fp, Tp, Tm01, Tm02 = scientimate.pmpsd(U195=10, fp=0.33, fs=2, N=256, CalSpectralSP='yes', dispout='no')

Description
-----------

Calculate Pierson-Moskowitz spectrum (power spectral density), (Pierson and Moskowitz 1964) 

Inputs
------

U195=10
    Wind velocity at 19.5 meter above surface level in (m/s)
fp=0.33
    | Peak wave frequency (fp=1/Tp) in (Hz)
    | If CalSpectralSP='yes'; then fp is calculated from U195
fs=8
    Sampling frequency that data collected at in (Hz)
N=256
    Total number of points between 0 and fs that spectrum reports at is (N+1)
CalSpectralSP='yes'
    Define to calculate spectral shape parameters or not ('yes': calculate, 'no': use given parameters by user)
dispout='no'
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

.. code:: python

    import scientimate as sm
    f,Syy,Hm0,fp,Tp,Tm01,Tm02=sm.pmpsd(10,0.33,2,256,'yes','yes')

References
----------

Pierson, W. J., & Moskowitz, L. (1964). 
A proposed spectral form for fully developed wind seas based on the similarity theory of SA Kitaigorodskii. 
Journal of geophysical research, 69(24), 5181-5190.

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
