.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-02-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

scientimate.incidentreflectedwave
=================================

.. code:: python

    Kr, EtaInc, EtaRef, aInc, aRef, t, f = scientimate.incidentreflectedwave(Eta1, Eta2, dx, h, fs=2, fmin=0, fmax=None, SepMethod='goda', kCalcMethod='beji', dispout='no')

Description
-----------

Separate incident and reflected waves

Inputs
------

Eta1
    Wave signal at position x1 in (m), should have same size as Eta2  
Eta2
    Wave signal at position x2 in (m), should have same size as Eta1 
dx
    Distance between x1 and x2 in (m), dx=x2-x1 
h
    Mean water depth in (m) 
fs=2
    Sampling frequency that data collected at in (Hz), if fs=1 then output is equivalent to normalized filter
fmin=0
    Minimum frequency to be considered in (Hz)
fmax=fs/2
    Maximum frequency to be considered in (Hz)
SepMethod='goda'
    | Incident and reflected waves Separation method 
    | 'goda': Goda and Suzuki (1977) 
    | 'ma': Ma et al. (2010)
    | 'frigaard': Frigaard and Brorsen (1995) 
kCalcMethod='beji'
    | Wave number calculation method 
    | 'hunt': Hunt (1979), 'beji': Beji (2013), 'vatankhah': Vatankhah and Aghashariatmadari (2013) 
    | 'goda': Goda (2010), 'exact': calculate exact value 
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

Kr
    Reflection coefficient
EtaInc
    Incident wave (m)
EtaRef
    Reflected wave (m)
aInc
    Amplitude of incident wave (m)
aRef
    Amplitude of reflected wave (m)
t
    Time (s)
f
    Frequency (Hz)

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np

    h=1
    fs=2
    dt=1/fs
    duration=1024
    t=np.linspace(0,duration-dt,duration*fs)
    x1=1 
    x2=3.7
    dx=x2-x1
    W1=0.5*np.cos(0.412*x1-2*np.pi*0.2*t)
    W2=0.5*np.cos(0.739*x1-2*np.pi*0.34*t)
    EtaIncGauge1=(W1+W2)
    W3=0.5*np.cos(0.412*x2-2*np.pi*0.2*t)
    W4=0.5*np.cos(0.739*x2-2*np.pi*0.34*t)
    EtaIncGauge2=(W3+W4)
    W5=0.1*np.cos(0.412*x1+2*np.pi*0.2*t)
    W6=0.1*np.cos(0.739*x1+2*np.pi*0.34*t)
    EtaRefGauge1=(W5+W6)
    W7=0.1*np.cos(0.412*x2+2*np.pi*0.2*t)
    W8=0.1*np.cos(0.739*x2+2*np.pi*0.34*t)
    EtaRefGauge2=(W7+W8)
    Eta1=EtaIncGauge1+EtaRefGauge1
    Eta2=EtaIncGauge2+EtaRefGauge2

    Kr,EtaInc,EtaRef,aInc,aRef,t,f=sm.incidentreflectedwave(Eta1,Eta2,dx,h,fs,0,fs/2,'goda','beji','yes')

References
----------

Beji, S. (2013). 
Improved explicit approximation of linear dispersion relationship for gravity waves. 
Coastal Engineering, 73, 11-12.

Baldock, T. E., & Simmonds, D. J. (1999). Separation of incident and reflected waves over sloping bathymetry. 
Coastal Engineering, 38(3), 167-176.

Frigaard, P., Brorsen, M., 1995. A time domain method for separating incident and reflected irregular waves.
Coastal Eng. 24, 205–215

Goda, Y., & Suzuki, Y. (1977). Estimation of incident and reflected waves in random wave experiments. 
In Coastal Engineering 1976 (pp. 828-845).

Goda, Y. (2010). 
Random seas and design of maritime structures. 
World scientific.

Hunt, J. N. (1979). 
Direct solution of wave dispersion equation. 
Journal of the Waterway Port Coastal and Ocean Division, 105(4), 457-459.

Ma, Y., Dong, G., Ma, X., & Wang, G. (2010). 
A new method for separation of 2D incident and reflected waves by the Morlet wavelet transform. 
Coastal Engineering, 57(6), 597-603.

Mansard, E. P., & Funke, E. R. (1980). 
The measurement of incident and reflected spectra using a least squares method. 
In Coastal Engineering 1980 (pp. 154-172).

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
