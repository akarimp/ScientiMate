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

tmapsd
======

.. code:: MATLAB

    [f, Syy, Hm0, fp, Tp, Tm01, Tm02, PHI] = tmapsd(U10, F, h, fp, fs, N, CalSpectralSP, transfCalcMethod, kCalcMethod, dispout)

Description
-----------

Calculate TMA spectrum (power spectral density), (Bouws et al. 1985)

Inputs
------

U10=10;
    Wind velocity at 10 meter above surface level in (m/s)
F=10000;
    Wind fetch length in (m)
h=5;
    Mean water depth in (m)
fp=0.33;
    | Peak wave frequency (fp=1/Tp) in (Hz)
    | If CalSpectralSP='yes'; then fp is calculated from U10 and F
fs=2;
    Sampling frequency that data collected at in (Hz)
N=256;
    Total number of points between 0 and fs that spectrum reports at is (N+1)
CalSpectralSP='yes';
    Define to calculate spectral shape parameters or not ('yes': calculate, 'no': use given parameters by user)
transfCalcMethod='approx';
    | Transformation function from JONSWAP into TMA calculation method 
    | 'approx': approximated method, 'tucker': Tucker (1994), 'kitaigordskii': Kitaigordskii et al. (1975) 
kCalcMethod='beji';
    | Wave number calculation method 
    | 'hunt': Hunt (1979), 'beji': Beji (2013), 'vatankhah': Vatankhah and Aghashariatmadari (2013) 
    | 'goda': Goda (2010), 'exact': calculate exact value 
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
PHI
    Transformation function from JONSWAP into TMA

Examples
--------

.. code:: MATLAB

    [f,Syy,Hm0,fp,Tp,Tm01,Tm02,PHI]=tmapsd(10,10000,5,0.33,2,256,'yes','approx','beji','yes');

References
----------

Beji, S. (2013). 
Improved explicit approximation of linear dispersion relationship for gravity waves. 
Coastal Engineering, 73, 11-12.

Bouws, E.; Günther, H.; Rosenthal, W., and Vincent, C.L., (1985). 
Similarity of the wind wave spectrum in finite depth water: 1. Spectral form. 
Journal of Geophysical Research: Oceans, 90(C1), 975-986.

Goda, Y. (2010). 
Random seas and design of maritime structures. 
World scientific.

Hunt, J. N. (1979). 
Direct solution of wave dispersion equation. 
Journal of the Waterway Port Coastal and Ocean Division, 105(4), 457-459.

Kitaigordskii, S. A., Krasitskii, V. P., & Zaslavskii, M. M. (1975). 
On Phillips' theory of equilibrium range in the spectra of wind-generated gravity waves. 
Journal of Physical Oceanography, 5(3), 410-420.

Vatankhah, A. R., & Aghashariatmadari, Z. (2013). 
Improved explicit approximation of linear dispersion relationship for gravity waves: A discussion. 
Coastal engineering, 78, 21-22.

Tucker, M. J. (1994). 
Nearshore waveheight during storms. 
Coastal Engineering, 24(1-2), 111-136.

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
