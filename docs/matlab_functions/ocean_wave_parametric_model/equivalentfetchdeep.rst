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

equivalentfetchdeep
===================

.. code:: MATLAB

    [FetchEq, FetchhatEq] = equivalentfetchdeep(windvel, SustWindDur, CalcMethod)

Description
-----------

Calculate an equivalent wind fetch for duration limited wave growth in deep water

Inputs
------

windvel
    | Wind velocity in (m/s)
    | Wind velocity should be measured (or represents velocity) at 10 m above surface
    | Wind velocity should be a 10 minutes averaged wind for all methods, except for for 'cem' and 'spmdeep' methods
    | For 'cem' and 'spmdeep' methods, wind velocity should be converted to duration of sustained wind by using gust factor
SustWindDur
    Sustained wind duration in (second)
CalcMethod='cem';
    | Parametric wave model to be used for calculation 
    | 'wilson': Use method by Wislon (1965)
    | 'jonswap': Use method by Hasselmann et al. (1973) known as JONSWAP
    | 'carter': Use method by Carter (1982)
    | 'spmdeep': Use method by Shore Protection Manual (SPM),
    |     U.S. Army Corps of Engineers (1984) in deep water
    | 'cem': Use method by Coastal Engineering Manual (CEM),
    |     U.S. Army Corps of Engineers (2015)

Outputs
-------

FetchEq
    Equivalent wind fetch values for a duration limited wave growth in (m)
FetchhatEq
    | Dimensionless equivalent wind fetch values for a duration limited wave growth: Fetchhat=g*Fetch/U10^2
    | Note, g=9.81: gravitational acceleration
    |       U10: wind velocity
    |       Fetch: wind fetch

Examples
--------

.. code:: MATLAB

    windvel(:,1)=10.*rand(100,1);
    SustWindDur(:,1)=10000.*rand(100,1);
    [FetchEq,FetchhatEq]=equivalentfetchdeep(windvel,SustWindDur,'cem');

References
----------

Carter, D. J. T. (1982). 
Prediction of wave height and period for a constant wind velocity using the JONSWAP results. 
Ocean Engineering, 9(1), 17-33.

Department of the Army, Waterways Experiment Station, Corps of Engineers, 
and Coastal Engineering Research Center (1984), 
Shore Protection Manual, Washington, 
D.C., vol. 1, 4th ed., 532 pp.

Hasselmann, K.; Barnett, T. P.; Bouws, E.; Carlson, H.; Cartwright, D. E.; Enke, K.; Ewing, J. A.; 
Gienapp, H.; Hasselmann, D. E.; Kruseman, P.; Meerbrug, A.; Muller, P.; Olbers, D. J.; Richter, K.; 
Sell, W., and Walden, H., (1973). 
Measurements of wind-wave growth and swell decay during the Joint North Sea Wave Project (JONSWAP). 
Deutsche Hydrographische Zeitschrift A80(12), 95p.

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
