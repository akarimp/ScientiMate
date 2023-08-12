.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-03-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

scientimate.diagnostictail
==========================

.. code:: python

    fOut, SyyOut, Hm0, fp, Tp, Tm01, Tm02 = scientimate.diagnostictail(fIn, SyyIn, ftail, tailtype='jonswap', tailpower=-5, h=0, transfCalcMethod='approx', kCalcMethod='beji', dispout='no')

Description
-----------

Replace a spectrum tail with JONSWAP (Hasselmann et al., 1973) or TMA Spectrum (Bouws et al., 1985)

Inputs
------

fIn
    Frequency (Hz), Input
SyyIn
    Power spectral density (m^2/Hz), Input
ftail
    Frequency that diagnostic tail apply after that (typically: ftail=2.5fm where fm=1/Tm01)
tailtype='jonswap'
    | Define type of the diagnostic tail to be applied 
    | 'jonswap': JONSWAP Spectrum tail, 'tma': TMA Spectrum tail
tailpower=-5
    | Tail power that diagnostic tail apply based on that 
    | tailpower=-3 for shallow water, tailpower=-5 for deep water
h=0
    Mean water depth in (m)
transfCalcMethod='approx'
    | Transformation function from JONSWAP into TMA calculation method 
    | 'approx': approximated method, 'tucker': Tucker (1994), 'kitaigordskii': Kitaigordskii et al. (1975) 
kCalcMethod='beji'
    | Wave number calculation method 
    | 'hunt': Hunt (1979), 'beji': Beji (2013), 'vatankhah': Vatankhah and Aghashariatmadari (2013) 
    | 'goda': Goda (2010), 'exact': calculate exact value 
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

fOut
    Frequency (Hz), Output
SyyOut
    Power spectral density (m^2/Hz) with diagnostic tail, Output
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
    import numpy as np

    N=2**11 #Total number of points
    fs=8 #Sampling frequency
    df=fs/N #Frequency difference 
    f=np.arange(0,fs/2+df,df) #Frequency vector 
    f[0]=f[1]/2 #Assign temporarily non-zero value to fisrt element of f to prevent division by zero
    Syy=0.016*9.81**2/((2*np.pi)**4*(f**5))*np.exp(-1.25*(0.33/f)**4) #Calculating Spectrum 
    f[0]=0
    Syy[0]=0
    
    fOut,SyyOut,Hm0,fp,Tp,Tm01,Tm02=sm.diagnostictail(f,Syy,0.5,'jonswap',-5,5,'approx','beji','yes')
    
    fOut,SyyOut,Hm0,fp,Tp,Tm01,Tm02=sm.diagnostictail(f,Syy,0.5,'tma',-3,5,'approx','beji','yes')

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

Hasselmann, K.; Barnett, T. P.; Bouws, E.; Carlson, H.; Cartwright, D. E.; Enke, K.; Ewing, J. A.; 
Gienapp, H.; Hasselmann, D. E.; Kruseman, P.; Meerbrug, A.; Muller, P.; Olbers, D. J.; Richter, K.; 
Sell, W., and Walden, H., (1973). 
Measurements of wind-wave growth and swell decay during the Joint North Sea Wave Project (JONSWAP). 
Deutsche Hydrographische Zeitschrift A80(12), 95p.

Hunt, J. N. (1979). 
Direct solution of wave dispersion equation. 
Journal of the Waterway Port Coastal and Ocean Division, 105(4), 457-459.

Vatankhah, A. R., & Aghashariatmadari, Z. (2013). 
Improved explicit approximation of linear dispersion relationship for gravity waves: A discussion. 
Coastal engineering, 78, 21-22.

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
