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

hurricanewindh80
================

.. code:: MATLAB

    [Vxgrid, Vygrid, Vgrid, Vgxgrid, Vgygrid, Vggrid, Pgrid, Rgrid, RVmax, thetagrid, thetaVgrid, thetaVgridtan] = hurricanewindh80(xgrid, ygrid, xCenter, yCenter, Pc, Vgmax, Rknown, VgRknown, Rmax, Pn, Rhoa, ...
    VgToVCoeff, inflowdirCalcMethod, Vt, VtAzmdir, backwindCalcMethod, distCalcMethod, flattendata, savedata, dispout)

Description
-----------

Generate hurricane wind and pressure data on given (x,y) points using method from Holland (1980)

Inputs
------

xgrid
    | x (longitude) of points which outputs are calculated at
    | xgrid can be a single point or 1d or 2d array 
ygrid
    | y (latitude) of points which outputs are calculated at
    | ygrid can be a single point or 1d or 2d array 
xCenter
    x (longitude) of hurricane center (track)
yCenter
    y (latitude) of hurricane center (track)
Pc
    Hurricane central surface pressure in (Pa)
Vgmax
    | Maximum hurricane 1-min averaged wind velocity at the gradient level in (m/s)
    | It can be estimated from Vmax=3.44*(1010-Pc*1e-2)^0.644; where Pc is in (Pa)
    | and, Vmax is 1-min averaged wind at a 10-m elevation in (m/s) (Atkinson & Holliday, 1977)
    | Hurricane translation velocity (forward velocity) needs to be removed from Vgmax before applying (e.g. Hu et al., 2012)
    | Hurricane translation velocity (forward velocity) needs to be added after rotational velocity is calculated
Rknown
    | Radius that hurricane wind velocity is known at that radius in (m)
    | Rknown should be larger than radius associated to Vgmax
VgRknown
    Hurricane wind velocity at the gradient level which is known at radius Rknown in (m/s)
Rmax=423e3;
    | Maximum radius of hurricane from hurricane center in (m)
    | Velocity outputs for points with R>Rmax is set to zero
    | Median values of Rmax is 423 km (e.g. Chavas & Emanuel 2010; Lin & Chavas, 2012)
Pn=101325;
    | Ambient surface pressure (external pressure) in (Pa)
    | Standard atmosphere pressure is 101325 (Pa) 
    | Typical value: Pn=101500 (Pa) for the western North Pacific, Pn= 101000 (Pa) for the North Atlantic
    | (Batke et al., 2014)
Rhoa=1.204;
    Air density at the gradient level in (kg/m3)
VgToVCoeff=0.8;
    | Coefficient to convert gradient wind velocity to wind velocity at 10 m above surface as: 
    | V=VgToVCoeff*Vg, if VgToVCoeff=1, then V=Vg
inflowdirCalcMethod='bretschneider';
    | Inflow angle calculation method 
    | 'no': Inflow angle are not calculated, thetaVgrid=thetaVgridtan
    | 'bretschneider': Inflow angle are calculated based on Bretschneider (1972)
    | 'sobey': Inflow angle are calculated based on Sobey et al. (1977)
Vt=0;
    Hurricane central translational velocity in (m/s)
VtAzmdir=0;
    | Hurricane center velocity azimuth (bearing) direction in (Degree)
    | azimuth (bearing) direction which is measured clockwise from the north:
    | 0 (degree): toward North, 90 (degree): toward East, 180 (degree): toward South, 270 (degree): toward West 
backwindCalcMethod='no';
    | Calculation method for adding background wind velocity due to hurricane motion
    | background wind velocity is added to points with R<=Rmax, e.g. Chavas & Emanuel (2010), Lin & Chavas (2012)
    | 'no': background wind velocity is not added to hurricane velocities
    | 'constant': background wind velocity is added as constant value to hurricane velocities
    | 'slosh': Background wind velocity is calculated and added using SLOSH model by Jelesnianski et al. (1992)
    | 'lin': Background wind velocity is calculated and added based on Lin & Chavas (2012)
distCalcMethod='gc';
    | Distance calculation method 
    | 'cart': Distances are calculated on cartesian coordinate
    | 'gc': Distances are calculated on Great Circle based on Vincenty formula, Vincenty (1975)
    | Earth radius coonsidered as mean earth radius=6371000 m
flattendata='no';
    | Define if flat data or not
    | 'no': does not flat the results, outputs are in 3d array
    | 'yes': flat the results, outputs are in 2d array
savedata='no';
    | Define if save data in a file or not in working folder
    | 'no': does not save, 
    | 'yes': save data as ascii 'dat' file, data are flatten regrdless of flattendata value
dispout='no';
    | Define to display outputs or not
    | 'imagesc': 2 dimensional plot using imagesc or imshow
    | 'pcolor': 2 dimensional plot using pcolor
    | 'contour': 2 dimensional contour plot, number of contour=ncolor
    | 'quiver': 2 dimensional vector plot 
    | 'no': not display 
    | Use dispout='no' if calculation mesh is not 2d array
    | if there is more than one time step, only the last one is plotted
    | if flattendata='yes'; then dispout is set as dispout='no';

Outputs
-------

Vxgrid
    | Hurricane 1-min averaged wind velocity at 10 m above surface in x (East) direction on defined mesh in (m/s)
    | Gradient wind velocity converted to wind velocity at 10 m above surface by V=VgToVCoeff*Vg
Vygrid
    | Hurricane 1-min averaged wind velocity at 10 m above surface in y (North) direction on defined mesh in (m/s)
    | Gradient wind velocity converted to wind velocity at 10 m above surface by V=VgToVCoeff*Vg
Vgrid
    | Resultant hurricane 1-min averaged wind velocity at 10 m above surface (Vx^2+Vy^2)^0.5 on defined mesh in (m/s)
    | Gradient wind velocity converted to wind velocity at 10 m above surface by V=VgToVCoeff*Vg
Vgxgrid
    Hurricane 1-min averaged gradient wind velocity at the gradient level in x (East) direction on defined mesh in (m/s)
Vgygrid
    Hurricane 1-min averaged gradient wind velocity at the gradient level in y (North) direction on defined mesh in (m/s)
Vggrid
    Resultant hurricane 1-min averaged gradient wind velocity at the gradient level on defined mesh in (m/s)
Pgrid
    Hurricane surface pressure on defined mesh in (Pa)
Rgrid
    Distance (radius) from hurricane center to each point on the grid
RVmax
    Distance (radius) from hurricane center to a location of maximum hurricane wind velocity (m)
thetagrid
    Angle from hurricane center to each point on the grid in (Degree)
thetaVgrid
    | Inflow angle (trigonometric direction) of hurricane velocity at each grid point in (Degree)
    | Inflow angle: angle between the inwardly spiraling surface wind 
    |     and the circular isobars around the hurricane center (Boose et al., 2004)
thetaVgridtan
    | Angle (trigonometric direction) of hurricane velocity at each grid point in (Degree)
    | thetaVgridtan is tangential angle respect to radius. 
    | Note: Outputs has dimension of [M,N,L] where [M,N] is size of the x-y grid and [L] is number of time steps
    |         If flattendata='yes'; then Outputs has dimension of [M*L,N]
    |     Hurricane translation velocity needs to be added after rotational velocity is calculated 
    |         (e.g. Hu et al., 2012; Lin & Chavas, 2012)
    |     Gradient wind velocity is converted to standard wind height as
    |         wind velocity at 10 m above surface by V=VgToVCoeff*Vg
    |     1-min averaged wind velocity needs to be converted to standard duration such as 
    |         10-min averaged wind by using a gust factor

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
    Vmax=Vmax-Vt; %Removing hurricane translation velocity from Vgmax
    Vgmax=Vmax/0.8; %Converting surface velocity to gradient velocity

    %34 kt (17.49 m/s) wind radii maximum extent in northeastern quadrant in (m) for Hurricane Katrina at max velocity
    Rknown=370400;
    VRknown=17.49;
    VRknown=VRknown-Vt; %Removing hurricane translation velocity from VRknown
    VgRknown=VRknown/0.8; %Converting surface velocity to gradient velocity

    Pn=101325; %Ambient surface pressure (external pressure) in (Pa)
    Rhoa=1.204; %Air density in (kg/m3)

    [Vxgrid,Vygrid,Vgrid,Vgxgrid,Vgygrid,Vggrid,Pgrid,Rgrid,RVmax,thetagrid,thetaVgrid,thetaVgridtan]=hurricanewindh80(xgrid,ygrid,longCenter,latCenter,Pc,Vgmax,Rknown,VgRknown,423e3,Pn,Rhoa,...
        0.8,'bretschneider',Vt,VtAzmdir,'slosh','gc','no','no','quiver');

    %Converting 1-min sustained wind to 10-min averaged wind using gust factor
    %e.g. World Meteorological Organization (2015)
    Vxgrid=Vxgrid*0.88;
    Vygrid=Vygrid*0.88;
    Vgrid=Vgrid*0.88;


    %EXAMPLE 2

    %Creating calculation mesh
    [xgrid,ygrid]=meshgrid(linspace(-98,-68,100),linspace(16,44,100));

    %Longitude of Hurricane Katrine best track
    longtrack=[-75.1;-75.7;-76.2;-76.5;-76.9;-77.7;-78.4;-79.0;-79.6;-80.1;-80.3;-81.3;...
        -82.0;-82.6;-83.3;-84.0;-84.7;-85.3;-85.9;-86.7;-87.7;-88.6;-89.2;-89.6;...
        -89.6;-89.6;-89.6;-89.6;-89.1;-88.6;-88.0;-87.0;-85.3;-82.9];

    %Latitude of Hurricane Katrine best track
    lattrack=[23.1;23.4;23.8;24.5;25.4;26.0;26.1;26.2;26.2;26.0;25.9;25.4;...
        25.1;24.9;24.6;24.4;24.4;24.5;24.8;25.2;25.7;26.3;27.2;28.2;...
        29.3;29.5;30.2;31.1;32.6;34.1;35.6;37.0;38.6;40.1];

    %Hurricane Katrina centeral pressure (Pa)
    Pc=[100800;100700;100700;100600;100300;100000;99700;99400;98800;98400;98300;98700;...
        97900;96800;95900;95000;94200;94800;94100;93000;90900;90200;90500;91300;...
        92000;92300;92800;94800;96100;97800;98500;99000;99400;99600];

    %Hurricane Katrina translational velocity (m/s)
    Vt=[0.00000;3.23091;3.13105;3.86928;4.99513;4.82816;3.27813;2.81998;2.77140;2.53041;...
        1.05928;5.30662;3.60661;2.98269;3.61863;3.43691;3.28168;2.85849;3.20404;4.26279;...
        5.31340;5.18467;5.39195;5.46121;5.66270;1.02958;3.60354;4.63312;8.02540;8.01558;...
        8.12721;8.31580;10.75406;12.28350];
        
    %Hurricane Katrina velocity azimuth (bearing) in (Degree)
    VtAzmdir=[0.00000;298.67291;311.22135;338.70264;338.13626;309.94476;279.18860;280.65053;270.13245;...
    246.10095;240.96690;241.20181;244.79591;249.93382;244.88325;252.71384;270.14459;280.49918;...
    298.94148;299.05364;299.18896;306.76219;329.36839;340.59069;0.00000;0.00000;0.00000;...
        0.00000;15.67775;15.42254;18.00215;29.63266;39.49673;50.29744];

    %Hurricane Katrina 1-min sustained maximum velocity (m/s)
    Vmax=[15.3;15.3;15.3;17.850;20.4;22.950;25.5;28.050;30.6;35.7;35.7;33.150;...
        38.250;43.350;45.9;48.450;51.0;51.0;51.0;63.750;73.950;76.5;71.4;63.750;...
        56.1;56.1;53.550;40.8;25.5;20.4;15.3;15.3;15.3;12.750];

    Vmax=Vmax-Vt; %Removing hurricane translation velocity from Vmax
    Vgmax=Vmax./0.8; %Converting surface velocity to gradient velocity

    %34 kt (17.49 m/s) wind radii maximum extent in northeastern quadrant in (m) for Hurricane Katrina
    RknownRaw=[0;0;0;111120;111120;111120;111120;111120;129640;NaN;129640;138900;...
        138900;138900;166680;240760;240760;259280;259280;296320;333360;370400;370400;370400;...
        NaN;370400;NaN;185200;138900;138900;0;0;0;0];

    %34 kt (17.49 m/s) wind radii maximum extent in northeastern quadrant in (m) for Hurricane Katrina
    Rknown=[0;0;0;111120;111120;111120;111120;111120;129640;129640;129640;138900;...
        138900;138900;166680;240760;240760;259280;259280;296320;333360;370400;370400;370400;...
        370400;370400;277800;185200;138900;138900;0;0;0;0];
    VRknown=ones(34,1).*17.49;
    VRknown=VRknown-Vt; %Removing hurricane translation velocity from VRknown
    VgRknown=VRknown/0.8; %Converting surface velocity to gradient velocity

    Pn=101325; %Ambient surface pressure (external pressure) in (Pa)
    Rhoa=1.204; %Air density in (kg/m3)

    [Vxgrid,Vygrid,Vgrid,Vgxgrid,Vgygrid,Vggrid,Pgrid,Rgrid,RVmax,thetagrid,thetaVgrid,thetaVgridtan]=hurricanewindh80(xgrid,ygrid,longtrack(4:27,1),lattrack(4:27,1),Pc(4:27,1),Vgmax(4:27,1),Rknown(4:27,1),VgRknown(4:27,1),423e3,Pn,Rhoa,...
        0.8,'bretschneider',Vt(4:27,1),VtAzmdir(4:27,1),'slosh','gc','no','no','quiver');

    %Converting 1-min sustained wind to 10-min averaged wind using gust factor
    %e.g. World Meteorological Organization (2015)
    Vxgrid=Vxgrid*0.88;
    Vygrid=Vygrid*0.88;
    Vgrid=Vgrid*0.88;


    EXAMPLE 3

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

    [Vxgrid,Vygrid,Vgrid,Vgxgrid,Vgygrid,Vggrid,Pgrid,Rgrid,RVmax,thetagrid,thetaVgrid,thetaVgridtan]=hurricanewindh80(xgrid,ygrid,longCenter,latCenter,Pc,Vgmax,Rknown,VgRknown,423e3,Pn,Rhoa,...
        0.8,'bretschneider',Vt,VtAzmdir,'slosh','gc','no','no','no');
    plot(Rgrid,Vgrid)

References
----------

Data

* www.nhc.noaa.gov/data/
* www.nhc.noaa.gov/data/hurdat/hurdat2-format-nencpac.pdf
* coast.noaa.gov/hurricanes
* www.aoml.noaa.gov/hrd/data_sub/re_anal.html

Atkinson, G. D., & Holliday, C. R. (1977). 
Tropical cyclone minimum sea level pressure/maximum sustained wind relationship for the western north Pacific. 
Monthly Weather Review, 105(4), 421-427.

Batke, S. P., Jocque, M., & Kelly, D. L. (2014). 
Modelling hurricane exposure and wind speed on a mesoclimate scale: a case study from Cusuco NP, Honduras. 
PloS one, 9(3), e91306.

Boose, E. R., Serrano, M. I., & Foster, D. R. (2004). 
Landscape and regional impacts of hurricanes in Puerto Rico. 
Ecological Monographs, 74(2), 335-352.

Bretschneider, C. L. (1972, January). 
A non-dimensional stationary hurricane wave model. 
In Offshore Technology Conference. Offshore Technology Conference.

Chavas, D. R., & Emanuel, K. A. (2010). 
A QuikSCAT climatology of tropical cyclone size. 
Geophysical Research Letters, 37(18).

Department of the Army, Waterways Experiment Station, Corps of Engineers, 
and Coastal Engineering Research Center (1984), 
Shore Protection Manual, Washington, 
D.C., vol. 1, 4th ed., 532 pp.

Graham and Numm (1959) 
Meteorological Conditions Pertinent to Standard Project Hurricane, Atlantic and Gulf Coasts of United States.
National Hurricane Research Project. U.S. Weather Service, Report no. 33.

Holland, G. J. (1980). 
An analytic model of the wind and pressure profiles in hurricanes. 
Monthly weather review, 108(8), 1212-1218.

Holland, G. (2008). 
A revised hurricane pressure–wind model. 
Monthly Weather Review, 136(9), 3432-3445.

Holland, G. J., Belanger, J. I., & Fritz, A. (2010). 
A revised model for radial profiles of hurricane winds. 
Monthly Weather Review, 138(12), 4393-4401.

Hu, K., Chen, Q., & Kimball, S. K. (2012). 
Consistency in hurricane surface wind forecasting: an improved parametric model. 
Natural hazards, 61(3), 1029-1050.

Jelesnianski, C. P., Chen, J., & Shaffer, W. A. (1992). 
SLOSH: Sea, lake, and overland surges from hurricanes (Vol. 48). 
US Department of Commerce, National Oceanic and Atmospheric Administration, National Weather Service.

Lin, N., & Chavas, D. (2012). 
On hurricane parametric wind and applications in storm surge modeling. 
Journal of Geophysical Research: Atmospheres, 117(D9).

Phadke, A. C., Martino, C. D., Cheung, K. F., & Houston, S. H. (2003). 
Modeling of tropical cyclone winds and waves for emergency management. 
Ocean Engineering, 30(4), 553-578.

Powell, M. D., Vickery, P. J., & Reinhold, T. A. (2003). 
Reduced drag coefficient for high wind speeds in tropical cyclones. 
Nature, 422(6929), 279.

Sobey, R. J., Harper, B. A., & Stark, K. P. (1977). 
Numerical simulation of tropical cyclone storm surge. 
James Cook University of North Queensland, Department of Civil & Systems Engineering.

U.S. Army Corps of Engineers (2015). 
Coastal Engineering Manual. 
Engineer Manual 1110-2-1100, Washington, D.C.: U.S. Army Corps of Engineers.

Valamanesh, V., Myers, A. T., Arwade, S. R., Hajjar, J. F., Hines, E., & Pang, W. (2016). 
Wind-wave prediction equations for probabilistic offshore hurricane hazard analysis. 
Natural Hazards, 83(1), 541-562.

Wei, K., Arwade, S. R., Myers, A. T., Valamanesh, V., & Pang, W. (2017). 
Effect of wind and wave directionality on the structural performance of non‐operational offshore wind turbines supported by jackets during hurricanes. 
Wind Energy, 20(2), 289-303.

World Meteorological Organization. Tropical Cyclone Programme, & Holland, G. J. (2015). 
Global guide to tropical cyclone forecasting. 
Secretariat of the World Meteorological Organization.

Young, I. R., & Vinoth, J. (2013). 
An 'extended fetch' model for the spatial distribution of tropical cyclone wind–waves as observed by altimeter. 
Ocean Engineering, 70, 14-24.

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
