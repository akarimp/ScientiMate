.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-10-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

hurricanewavecontourh16
=======================

.. code:: MATLAB

    [Hsgrid, Tpgrid, Hsmax, Tpmax] = hurricanewavecontourh16(xgrid, ygrid, Vgrid, xCenter, yCenter, VtAzmdir, G, fetchCalcMethod, distCalcMethod, dispout)

Description
-----------

| Calculates hurricane wave height and wave period field (contours) on given mesh using method from Hwang (2016) and Hwang & Walsh (2016) 
| For method from Young (1988) use hurricanewavecontoury88
| For method from U.S. Army Corps of Engineers use hurricanewavecontourhcem

Inputs
------

xgrid
    | x (longitude) of points which outputs are calculated at
    | xgrid can be a single point or 1d or 2d array 
ygrid
    | y (latitude) of points which outputs are calculated at
    | ygrid can be a single point or 1d or 2d array 
Vgrid
    Resultant hurricane 1-min averaged wind velocity at 10 m above surface (Vx^2+Vy^2)^0.5 on defined mesh in (m/s)
xCenter
    x (longitude) of hurricane center (track)
yCenter
    y (latitude) of hurricane center (track)
VtAzmdir=0;
    | Hurricane center velocity azimuth (bearing) direction in (Degree)
    | azimuth (bearing) direction which is measured clockwise from the north:
    | 0 (degree): toward North, 90 (degree): toward East, 180 (degree): toward South, 270 (degree): toward West 
G=0.88;
    | Wind gust factor to convert 1-min averaged wind to 10-min averaged wind
    | e.g. Young (2017); Liu et al. (2017)
    | G=U(600s)/U(60s), therefore U(600s)=G*U(60s)
    | G=1/1.1; based on Powell and Houston (1996)
    | G=1/1.08; based on Harper (2013)
    | G=0.88; based on World Meteorological Organization (2015)
    | G=1/1.08 to 1/1.16; based on Liu et al. (2017)
fetchCalcMethod='constant';
    | Effective wind fetch calculation method 
    | 'constant': Use constant coefficients needed for effective wind fetch calculation
    | 'interp': Use interpolateed coefficients needed for effective wind fetch calculation
    | Earth radius coonsidered as mean earth radius=6371000 m
distCalcMethod='gc';
    | Distance calculation method 
    | 'cart': Distances are calculated on cartesian coordinate
    | 'gc': Distances are calculated on Great Circle based on Vincenty formula, Vincenty (1975)
    | Earth radius coonsidered as mean earth radius=6371000 m
dispout='no';
    | Define to display outputs or not
    | 'imagesc': 2 dimensional plot using imagesc or imshow
    | 'pcolor': 2 dimensional plot using pcolor
    | 'contour': 2 dimensional contour plot, number of contour=ncolor
    | 'no': not display 
    | Use dispout='no' if calculation mesh is not 2d array

Outputs
-------

Hsgrid
    Hurricane significant wave height on grid mesh in (m)
Tpgrid
    Hurricane peak wave period on grid mesh in (s)
Hsmax
    Hurricane maximum significant wave height in (m) 
Tpmax
    | Hurricane maximum peak wave period in (s) 
    | Note: Maximum values of wave height and wave period should be limited to fully developed values

Examples
--------

.. code:: MATLAB

    %EXAMPLE 1

    %Creating calculation mesh
    [xgrid,ygrid]=meshgrid(linspace(-98,-68,100),linspace(16,44,100));

    %Longitude of Hurricane Katrine center at max velocity
    longCenter=-88.6;

    %Latitude of Hurricane Katrine center at max velocity
    latCenter=26.3;

    %Hurricane Katrina centeral pressure (Pa) at max velocity
    Pc=90200;

    %Hurricane Katrina translational velocity (m/s) at max velocity
    Vt=5.18467;

    %Hurricane Katrina velocity azimuth (bearing) in (Degree) at max velocity
    VtAzmdir=306.76219;

    %Hurricane Katrina 1-min sustained maximum velocity (m/s) at max velocity
    Vmax=76.5;
    Vmax=Vmax-Vt; %Removing hurricane translation velocity from Vmax
    Vgmax=Vmax/0.8; %Converting surface velocity to gradient velocity

    %34 kt (17.49 m/s) wind radii maximum extent in northeastern quadrant in (m) for Hurricane Katrina at max velocity
    Rknown=370400;
    VRknown=17.49;
    VRknown=VRknown-Vt; %Removing hurricane translation velocity from VRknown
    VgRknown=VRknown/0.8; %Converting surface velocity to gradient velocity

    Pn=101325; %Ambient surface pressure (external pressure) in (Pa)
    Rhoa=1.204; %Air density in (kg/m3)

    %Calculating distance (radius) from hurricane center to each point
    Rgrid=(acos(sin(deg2rad(latCenter)).*sin(deg2rad(ygrid))+cos(deg2rad(latCenter)).*cos(deg2rad(ygrid)).*cos(deg2rad(xgrid)-deg2rad(longCenter)))).*6371000;

    %Generting wind velocity for Hurricane Katrine at max velocity using SLOSH model
    Vggrid=Vgmax.*(2.*32197.*Rgrid)./(32197^2+Rgrid.^2); %Gradient wind velocity
    Vggrid(Rgrid>=423e3)=0; 
    Vgrid=Vggrid.*0.8; %Wind velocity at 10m height

    [Hsgrid,Tpgrid,Hsmax,Tpmax]=hurricanewavecontourh16(xgrid,ygrid,Vgrid,longCenter,latCenter,VtAzmdir,0.88,'constant','gc','contour');

    %EXAMPLE 2

    xgrid=linspace(0,10,100); %(Degree)
    ygrid=ones(1,100).*20; %(Degree)
    longCenter=0; %(Degree)
    latCenter=20; %(Degree)
    Pc=90200; %(Pa)
    Vt=5.18467; %(m/s)
    VtAzmdir=306.76219; %(Degree) 
    Vmax=76.5; %(m/s)
    Vmax=Vmax-Vt;
    Vgmax=Vmax/0.8; %(m/s)
    Rknown=370400; %(m)
    VRknown=17.49; %(m/s)
    VRknown=VRknown-Vt;
    VgRknown=VRknown/0.8; %(m/s)
    Pn=101325; %Ambient surface pressure (external pressure) in (Pa)
    Rhoa=1.204; %Air density in (kg/m3)
    Rgrid=(acos(sin(deg2rad(latCenter)).*sin(deg2rad(ygrid))+cos(deg2rad(latCenter)).*cos(deg2rad(ygrid)).*cos(deg2rad(xgrid)-deg2rad(longCenter)))).*6371000;
    Vggrid=Vgmax.*(2.*32197.*Rgrid)./(32197^2+Rgrid.^2); %Gradient wind velocity
    Vggrid(Rgrid>=423e3)=0; 
    Vgrid=Vggrid.*0.8; %Wind velocity at 10m height

    [Hsgrid,Tpgrid,Hsmax,Tpmax]=hurricanewavecontourh16(xgrid,ygrid,Vgrid,longCenter,latCenter,VtAzmdir,0.88,'constant','gc','no');
    plot(Rgrid,Hsgrid)

References
----------

Data

* www.nhc.noaa.gov/data/
* www.nhc.noaa.gov/data/hurdat/hurdat2-format-nencpac.pdf
* coast.noaa.gov/hurricanes
* www.aoml.noaa.gov/hrd/data_sub/re_anal.html

Harper, B.A. (2013)
Best practice in tropical cyclone wind hazard modelling: In search of data and emptying the skeleton cupboard. 
In Proceedings of the 16th Australasian Wind Engineering Society Workshop, Brisbane, Qld, Australia, 18–19 July 2013

Hwang, P. A. (2016). 
Fetch-and duration-limited nature of surface wave growth inside tropical cyclones: With applications to air–sea exchange and remote sensing. 
Journal of Physical Oceanography, 46(1), 41-56.

Hwang, P. A., & Walsh, E. J. (2016). 
Azimuthal and radial variation of wind-generated surface waves inside tropical cyclones. 
Journal of Physical Oceanography, 46(9), 2605-2621.

Liu, Q., Babanin, A., Fan, Y., Zieger, S., Guan, C., & Moon, I. J. (2017). 
Numerical simulations of ocean surface waves under hurricane conditions: Assessment of existing model performance. 
Ocean Modelling, 118, 73-93.

Powell, M. D., & Houston, S. H. (1996). 
Hurricane Andrew's landfall in South Florida. Part II: Surface wind fields and potential real-time applications. 
Weather and Forecasting, 11(3), 329-349.

World Meteorological Organization. Tropical Cyclone Programme, & Holland, G. J. (2015). 
Global guide to tropical cyclone forecasting. 
Secretariat of the World Meteorological Organization.

Young, I.R. (2017)
A Review of Parametric Descriptions of Tropical Cyclone Wind-Wave Generation.
Atmosphere 2017, 8, 194.

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
