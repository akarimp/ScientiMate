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

mindurationshallow
==================

.. code:: MATLAB

    [tmin, tminhat, Fetchhat, hhat] = mindurationshallow(windvel, Fetch, hmean, CalcMethod, dispout)

Description
-----------

Calculate a minimum required wind duration for wave to be fetch limited in shallow and intermediate water

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
CalcMethod='karimpour';
    | Parametric wave model to be used for calculation 
    | 'spmshallow': Use method by Shore Protection Manual (SPM),
    |     U.S. Army Corps of Engineers (1984) in shallow and intermediate water water
    | 'karimpour': Use method by Karimpour et al. (2017)
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

tmin
    Minimum required wind duration for wind to be fetch limited in (second)
tminhat
    Dimensionless minimum required wind duration for wind to be fetch limited: tminhat=g*tmin/U10
Fetchhat
    Dimensionless wind fetch: Fetchhat=g*Fetch/U10^2
hhat
    | Dimensionless mean water depth along a wind fetch: hhat=g*hmean/U10^2
    | Note, g=9.81: gravitational acceleration
    |     U10: wind velocity

Examples
--------

.. code:: MATLAB

    windvel(:,1)=10.*rand(100,1);
    Fetch(:,1)=10000.*rand(100,1);
    hmean(:,1)=3.*rand(100,1);
    [tmin,tminhat,Fetchhat,hhat]=mindurationshallow(windvel,Fetch,hmean,'karimpour','no');

    windvel=10;
    Fetch(:,1)=[1e3:1000:1e6];
    hmean=3;
    [tmin,tminhat,Fetchhat,hhat]=mindurationshallow(windvel,Fetch,hmean,'karimpour','yes');

References
----------

U.S. Army Corps of Engineers (2015). 
Coastal Engineering Manual. 
Engineer Manual 1110-2-1100, Washington, D.C.: U.S. Army Corps of Engineers.

Department of the Army, Waterways Experiment Station, Corps of Engineers, 
and Coastal Engineering Research Center (1984), 
Shore Protection Manual, Washington, 
D.C., vol. 1, 4th ed., 532 pp.

Karimpour, A., Chen, Q., & Twilley, R. R. (2017). 
Wind Wave Behavior in Fetch and Depth Limited Estuaries. 
Scientific reports, 7, 40654.

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
