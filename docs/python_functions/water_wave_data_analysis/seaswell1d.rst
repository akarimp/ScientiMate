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

scientimate.seaswell1d
======================

.. code:: python

    fseparation, Hm0, Hm0sea, Hm0swell, Tp, Tpsea, Tpswell, m0, m0sea, m0swell, fp, Tm01, Tm02 = scientimate.seaswell1d(f, Syy, fsepCalcMethod='hwang', fu=0.5, fmaxswell=0.2, fpminswell=0, Windvel=10, dispout='no')

Description
-----------

Partition (separate) wind sea from swell in a power spectral density using an one dimensional method

Inputs
------

f
    Frequency (Hz)
Syy
    Power spectral density (m^2/Hz)
fsepCalcMethod='hwang'
    | Wind sea swell separating calculation method 
    | 'celerity': using deep water wave celerity, 'gilhousen': Gilhousen and Hervey (2001), 
    | 'hwang': Hwang et al. (2012), 'exact': calculate exact value 
fu=0.5
    An upper bound of a spectrum integration frequency (Hz)
fmaxswell=0.25
    Maximum frequency that swell can have, It is about 0.2 in Gulf of Mexico
fpminswell=0.1
    A lower bound of a spectrum (minimum frequency) that is used for Tpswell calculation
Windvel=10
    Wind velocity (m/s), only required for Gilhousen and Hervey (2001) method
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

fseparation
    Wind sea and Swell Separation Frequency (Hz)
Hm0
    Zero-Moment Wave Height (m)
Hm0sea
    Sea Zero-Moment Wave Height (m)
Hm0swell
    Swell Zero-Moment Wave Height (m)
Tp
    Peak wave period (second)
Tpsea
    Peak Sea period (second)
Tpswell
    Peak Swell Period (second)
m0
    Zero-Moment of the power spectral density (m^2)
m0sea
    Zero-Moment of the wind sea spectrum (m^2)
m0swell
    Zero-Moment of the swell spectrum (m^2)
fp
    Peak wave frequency (Hz)
Tm01
    Wave Period from m01 (second), Mean Wave Period
Tm02
    Wave Period from m02 (second), Mean Zero Crossing Period

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np

    N=2**11 #Total number of points
    fs=2 #Sampling frequency
    df=fs/N #Frequency difference 
    f=np.arange(0,fs/2+df,df) #Frequency vector 
    f[0]=f[1]/2 #Assign temporarily non-zero value to fisrt element of f to prevent division by zero
    SyySea=0.016*9.81**2/((2*np.pi)**4*(f**5))*np.exp(-1.25*(0.33/f)**4) #Calculating Spectrum 
    SyySwell=0.016*9.81**2/((2*np.pi)**4*(f**5))*np.exp(-1.25*(0.15/f)**4)*0.005 #Calculating Spectrum 
    Syy=SyySea+SyySwell
    f[0]=0
    Syy[0]=0
    fseparation,Hm0,Hm0sea,Hm0swell,Tp,Tpsea,Tpswell,m0,m0sea,m0swell,fp,Tm01,Tm02=sm.seaswell1d(f,Syy,'exact',0.5,0.3,0,10,'yes')

References
----------

Gilhousen, D. B., & Hervey, R. (2002). 
Improved estimates of swell from moored buoys. 
In Ocean Wave Measurement and Analysis (2001) (pp. 387-393).

Hwang, P. A., Ocampo-Torres, F. J., & García-Nava, H. (2012). 
Wind sea and swell separation of 1D wave spectrum by a spectrum integration method. 
Journal of Atmospheric and Oceanic Technology, 29(1), 116-128.

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
