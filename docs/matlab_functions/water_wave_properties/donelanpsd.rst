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

donelanpsd
==========

.. code:: MATLAB

    [f, Syy, Hm0, fp, Tp, Tm01, Tm02] = donelanpsd(U10, F, fp, fs, N, CalSpectralSP, dispout)

Description
-----------

Calculate Donelan spectrum (power spectral density), (Donelan et al. 1985)

Inputs
------

U10=10;
    Wind velocity at 10 meter above surface level in (m/s)
F=10000;
    Wind fetch length in (m)
fp=0.33;
    | Peak wave frequency (fp=1/Tp) in (Hz)
    | If CalSpectralSP='yes'; then fp is calculated from U10 and F
fs=2;
    Sampling frequency that data collected at in (Hz)
N=256;
    | Total number of points between 0 and fs-df, where df=fs/N
    | Spectrum is reported between 0 and fs/2 with total number of points equal to (N/2+1)
    | Total number of points between 0 and fs is (N+1)
CalSpectralSP='yes';
    Define to calculate spectral shape parameters or not ('yes': calculate, 'no': use given parameters by user)
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

    [f,Syy,Hm0,fp,Tp,Tm01,Tm02]=donelanpsd(10,10000,0.33,2,256,'yes','yes');

References
----------

Donelan, M.A.; Hamilton, J., and Hui, W., 1985. Directional spectra of wind-generated waves. 
Philosophical Transactions of the Royal Society of London A: Mathematical, 
Physical and Engineering Sciences, 315(1534), 509–562.

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
