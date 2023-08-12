.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-01-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

scientimate.wavebedstress
=========================

.. code:: python

    Tau, Tauc, Tauw, Ust, z0 = scientimate.wavebedstress(h, heightfrombed, d50, Uc=0, H=0, T=0, Rho=1000, kCalcMethod='beji')

Description
-----------

Calculate the bottom shear velocity and shear stress from current velocity and wave 

Inputs
------

h
    Water depth in (m) 
heightfrombed
    Sensor height from bed in (m)
d50
    Median bed particle diameter in (m)
Uc=0
    Current velocity at sensor height in (m/s), set equal to 0 if not exist 
H=0
    | Wave height in (m), H can be approximated and replaced by Hrms for random waves
    | Hrms: root mean square wave height, Hrms=Hm0/sqrt(2)=Hs/sqrt(2) 
    | Hm0: zero moment wave height, Hs: significant wave height
    | set equal to 0 if not exist  
T=0
    | Wave period in (s), T can be approximated and replaced by mean wave period for random waves 
    | if peak wave frequency (Tp) is used, calculated values represent peak wave 
    | set equal to 0 if not exist  
Rho=1000
    Water density in (kg/m^3)
kCalcMethod='beji'
    | Wave number calculation method 
    | 'hunt': Hunt (1979), 'beji': Beji (2013), 'vatankhah': Vatankhah and Aghashariatmadari (2013) 
    | 'goda': Goda (2010), 'exact': calculate exact value 
    | Note: inputs can be as a single value or a 1-D vertical array

Outputs
-------

Tau
    Total bottom shear stress from current velocity and wave
Tauc
    Bottom shear stress from current velocity
Tauw
    Bottom shear stress from wave
Ust
    Bottom shear velocity (m/s)
z0
    Surface roughness in (N/m^2)

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np

    Tau,Tauc,Tauw,Ust,z0=sm.wavebedstress(2,1.1,0.0000245,1.5,0.5,3,1000,'beji')

    Tau,Tauc,Tauw,Ust,z0=sm.wavebedstress(2,1.1,0.0000245,[1.5,2],[0.5,0.6],[3,3.1],1000,'exact')

    Tau,Tauc,Tauw,Ust,z0=sm.wavebedstress(2,1.1,0.0000245,np.array([1.5,2]),np.array([0.5,0.6]),np.array([3,3.1]),1000,'exact')

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
