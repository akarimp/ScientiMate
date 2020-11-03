.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-09-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

scientimate.parametricwavedeep
==============================

.. code:: python

    H, T, Ehat, fphat, m0, Fetchhat = parametricwavedeep(windvel, Fetch, CalcMethod='jonswap', dispout='no')

Description
-----------

Calculate wave properties using parametric wave models in deep water

Inputs
------

windvel
    | Wind velocity in (m/s)
    | Wind velocity should be measured (or represents velocity) at 10 m above surface
    | Wind velocity should be a 10 minutes averaged wind for all methods, except for for 'cem' and 'spmdeep' methods
    | For 'cem' and 'spmdeep' methods, wind velocity should be converted to duration of sustained wind by using gust factor
Fetch
    Wind fetch in (m)
CalcMethod='jonswap'
    | Parametric wave model to be used for wave calculation 
    | 'wilson': Use method by Wislon (1965)
    | 'jonswap': Use method by Hasselmann et al. (1973) known as JONSWAP
    | 'spmdeep': Use method by Shore Protection Manual (SPM),
    |     U.S. Army Corps of Engineers (1984) in deep water
    | 'kahma': Use method by Kahma and Calkoen (1992)
    | 'hwang': Use method by Hwang and Wang (2004)
    | 'cem': Use method by Coastal Engineering Manual (CEM),
    |     U.S. Army Corps of Engineers (2015)
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

H
    | Predicted wave height in (m) 
    | For all methods excepth for 'wilson': H=Hm0 where, Hm0 is a zero-moment wave height
    | For 'wilson' method: H=Hs, where Hs is a significant wave height
T
    | Predicted wave period in (s) 
    | For all methods excepth for 'wilson': T=Tp where, Tp is a peak wave period
    | For 'wilson' method: T=Ts, where Ts is a significant wave period (Ts=0.95Tp)
Ehat
    Predicted dimensionless wave energy, Ehat=g^2*m0/U10^4
fphat
    Predicted dimensionless peak wave frequency, fphat=fp*U10/g
m0
    Zero-moment of water surface elevation power spectral density in (m^2)
Fetchhat
    | Dimensionless Fetch: Fetchhat=g*Fetch/U10^2
    | Note, g=9.81: gravitational acceleration
    |     U10: wind velocity
    |     fp: peak wave frequency

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np

    windvel=10*np.random.rand(100)
    Fetch=10000*np.random.rand(100)
    H,T,Ehat,fphat,m0,Fetchhat=sm.parametricwavedeep(windvel,Fetch,'jonswap','no')

    windvel=10
    Fetch=np.arange(1e3,1e6,1000)
    H,T,Ehat,fphat,m0,Fetchhat=sm.parametricwavedeep(windvel,Fetch,'jonswap','yes')

References
----------

Department of the Army, Waterways Experiment Station, Corps of Engineers, 
and Coastal Engineering Research Center (1984), 
Shore Protection Manual, Washington, 
D.C., vol. 1, 4th ed., 532 pp.

Hasselmann, K.; Barnett, T. P.; Bouws, E.; Carlson, H.; Cartwright, D. E.; Enke, K.; Ewing, J. A.; 
Gienapp, H.; Hasselmann, D. E.; Kruseman, P.; Meerbrug, A.; Muller, P.; Olbers, D. J.; Richter, K.; 
Sell, W., and Walden, H., (1973). 
Measurements of wind-wave growth and swell decay during the Joint North Sea Wave Project (JONSWAP). 
Deutsche Hydrographische Zeitschrift A80(12), 95p.

Holthuijsen, L. H. (2007). 
Waves in oceanic and coastal waters. 
Cambridge university press.

Hwang, P. A., & Wang, D. W. (2004). 
Field measurements of duration-limited growth of wind-generated ocean surface waves at young stage of development. 
Journal of Physical Oceanography, 34(10), 2316-2326.

Kahma, K. K., & Calkoen, C. J. (1992). 
Reconciling discrepancies in the observed growth of wind-generated waves. 
Journal of Physical Oceanography, 22(12), 1389-1405.

Pierson, W. J., & Moskowitz, L. (1964). 
A proposed spectral form for fully developed wind seas based on the similarity theory of SA Kitaigorodskii. 
Journal of geophysical research, 69(24), 5181-5190.

U.S. Army Corps of Engineers (2015). 
Coastal Engineering Manual. 
Engineer Manual 1110-2-1100, Washington, D.C.: U.S. Army Corps of Engineers.

Wilson, B. W. (1965). 
Numerical prediction of ocean waves in the North Atlantic for December, 1959. 
Ocean Dynamics, 18(3), 114-130.

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
