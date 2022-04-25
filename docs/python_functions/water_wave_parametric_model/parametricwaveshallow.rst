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

scientimate.parametricwaveshallow
=================================

.. code:: python

    H, T, Ehat, fphat, m0, Fetchhat, hhat = scientimate.parametricwaveshallow(windvel, Fetch, hmean, CalcMethod='karimpour', dispout='no')

Description
-----------

Calculate wave properties using parametric wave models in shallow and intermediate water

Inputs
------

windvel
    | Wind velocity in (m/s)
    | Wind velocity should be measured (or represents velocity) at 10 m above surface
    | Wind velocity should be a 10 minutes averaged wind for all methods, except for for 'spmshallow' methods
    | For 'spmshallow' methods, wind velocity should be converted to duration of sustained wind by using gust factor
Fetch
    Wind fetch in (m)
hmean
    Mean water depth along a wind fetch in (m)
CalcMethod='karimpour'
    | Parametric wave model to be used for wave calculation 
    | 'spmshallow': Use method by Shore Protection Manual (SPM),
    |     U.S. Army Corps of Engineers (1984) in shallow and intermediate water water
    | 'young': Use method by Young and Verhagen (1996) and Young and Babanin (2006)
    | 'karimpour': Use method by Karimpour et al. (2017)
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

H
    | Predicted wave height in (m) 
    | For all methods: H=Hm0 where, Hm0 is a zero-moment wave height
T
    | Predicted wave period in (s) 
    | For all methods excepth for 'spmshallow': T=Tp where, Tp is a peak wave period
    | For 'spmshallow' method: T=Ts, where Ts is a significant wave period (Ts=0.95Tp)
Ehat
    Predicted dimensionless wave energy, Ehat=g^2*m0/U10^4
fphat
    Predicted dimensionless peak wave frequency, fphat=fp*U10/g
m0
    Zero-moment of water surface elevation power spectral density in (m^2)
Fetchhat
    Dimensionless wind fetch: Fetchhat=g*Fetch/U10^2
hhat
    | Dimensionless mean water depth along a wind fetch: hhat=g*hmean/U10^2
    | Note, g=9.81: gravitational acceleration
    |     U10: wind velocity
    |     fp: peak wave frequency
    |     hmean: mean depth along a wind fetch

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np

    windvel=10*np.random.rand(100)
    Fetch=10000*np.random.rand(100)
    hmean=3*np.random.rand(100)
    H,T,Ehat,fphat,m0,Fetchhat,hhat=sm.parametricwaveshallow(windvel,Fetch,hmean,'karimpour','no')

    windvel=10
    Fetch=np.arange(1e3,1e6,1000)
    hmean=3
    H,T,Ehat,fphat,m0,Fetchhat,hhat=sm.parametricwaveshallow(windvel,Fetch,hmean,'karimpour','yes')

References
----------

Bretschneider, C. L. (1952). 
Revised wave forecasting relationships. 
Coastal Engineering Proceedings, 1(2), 1.

Bretschneider, C. L. (1958). 
Revisions in wave forecasting: deep and shallow water. 
Coastal Engineering Proceedings, 1(6), 3.

Department of the Army, Waterways Experiment Station, Corps of Engineers, 
and Coastal Engineering Research Center (1984), 
Shore Protection Manual, Washington, 
D.C., vol. 1, 4th ed., 532 pp.

Karimpour, A., Chen, Q., & Twilley, R. R. (2017). 
Wind Wave Behavior in Fetch and Depth Limited Estuaries. 
Scientific reports, 7, 40654.

Pierson, W. J., & Moskowitz, L. (1964). 
A proposed spectral form for fully developed wind seas based on the similarity theory of SA Kitaigorodskii. 
Journal of geophysical research, 69(24), 5181-5190.

Sverdrup, H. U., & Munk, W. H. (1947). 
Wind, sea, and swell: theory of relations for forecasting. 
U.S. Navy Department, Hydrographic Office, Publication No. 601, 44 pp. 

Young, I. R., & Verhagen, L. A. (1996). 
The growth of fetch limited waves in water of finite depth. Part 1. Total energy and peak frequency. 
Coastal Engineering, 29(1-2), 47-78.

Young, I. R., & Babanin, A. V. (2006). 
The form of the asymptotic depth‐limited wind wave frequency spectrum. 
Journal of Geophysical Research: Oceans, 111(C6).

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
