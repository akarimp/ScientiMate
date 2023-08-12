.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-01-0                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

scientimate.stokeswavegenerator
===============================

.. code:: python

    Eta, t, Etaij = scientimate.stokeswavegenerator(h, amin, amax, Tmin, Tmax, Phimin=0, Phimax=2*3.1416, fs=32, duration=10, NoOfWave=2, kCalcMethod='beji', dispout='no')

Description
-----------

Generate second order stokes' waves

Inputs
------

h
    Mean water depth in (m)
amin
    Min wave amplitude in (m)
amax
    Max wave amplitude in (m)
Tmin
    Min wave mean period in (s)
Tmax
    Max wave mean period in (s)
Phimin=0
    Min Phase (radian)
Phimax=2*pi
    Max Phase (radian) 
fs=32
    Sample generation frequency (Hz), number of data points in one second
duration=10
    Duration time that data will be generated in (s)
NoOfWave=2
    Number of waves to be combined with each other
kCalcMethod='beji'
    | Wave number calculation method 
    | 'hunt': Hunt (1979), 'beji': Beji (2013), 'vatankhah': Vatankhah and Aghashariatmadari (2013) 
    | 'goda': Goda (2010), 'exact': calculate exact value 
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

Eta
    Water Surface Level Time Series in (m)
t
    Time in (s)
Etaij
    Separated Water Surface Level Time Series in (m)

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np

    Eta,t,Etaij=sm.stokeswavegenerator(5,0.2,0.4,1,3,0,2*np.pi,32,10,2,'beji','yes')

References
----------

Beji, S. (2013). 
Improved explicit approximation of linear dispersion relationship for gravity waves. 
Coastal Engineering, 73, 11-12.

Goda, Y. (2010). 
Random seas and design of maritime structures. 
World scientific.

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
