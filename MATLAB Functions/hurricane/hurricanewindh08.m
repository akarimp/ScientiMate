function [Vxgrid, Vygrid, Vgrid, Pgrid, Rgrid, RVmax, thetagrid, thetaVgrid, thetaVgridtan] = hurricanewindh08(xgrid, ygrid, xCenter, yCenter, Pc, Vmax, Rknown, VgRknown, dPcdt, Vt, Rmax, Pn, Rhoa, SST, ...
    VgToVCoeff, inflowdirCalcMethod, VtAzmdir, backwindCalcMethod, distCalcMethod, flattendata, savedata, dispout)
%{
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

hurricanewindh08
================

.. code:: MATLAB

    [Vxgrid, Vygrid, Vgrid, Pgrid, Rgrid, RVmax, thetagrid, thetaVgrid, thetaVgridtan] = hurricanewindh08(xgrid, ygrid, xCenter, yCenter, Pc, Vmax, Rknown, VgRknown, dPcdt, Vt, Rmax, Pn, Rhoa, SST, ...
    VgToVCoeff, inflowdirCalcMethod, VtAzmdir, backwindCalcMethod, distCalcMethod, flattendata, savedata, dispout)

Description
-----------

Generate hurricane wind and pressure data on (x,y) points using method from Holland (2008)

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
Vmax
    | Maximum hurricane 1-min averaged surface wind velocity at 10m height from surface in (m/s)
    | It can be estimated by Vmax=3.44*(1010-Pc*1e-2)^0.644; 
    | where, Vmax is 1-min averaged wind at a 10-m elevation in (m/s) (Atkinson & Holliday, 1977)
    | Hurricane translation velocity (forward velocity) needs to be removed from Vmax before applying (e.g. Hu et al., 2012)
    | Hurricane translation velocity (forward velocity) needs to be added after rotational velocity is calculated
Rknown
    | Radius that hurricane wind velocity is known at that radius in (m)
    | Rknown should be larger than radius associated to Vmax
VgRknown
    Hurricane wind velocity at the gradient level which is known at radius Rknown in (m/s)
dPcdt
    | Hurricane central pressure (Pc) intensity change (dPc/dt) in (hPa/hr)
    | Convert from (Pa) to (hPa): hPa=0.01*Pa;
Vt
    Hurricane translation velocity in (m/s)
Rmax=423e3;
    | Maximum radius of hurricane from hurricane center in (m)
    | Outputs for points with R>Rmax is set to zero
    | Median values of Rmax is 423 km (e.g. Chavas & Emanuel 2010; Lin & Chavas, 2012)
Pn=101325;
    | Ambient surface pressure (external pressure) in (Pa)
    | Standard atmosphere pressure is 101325 (Pa) 
    | Typical value: Pn=101500 (Pa) for the western North Pacific, Pn= 101000 (Pa) for the North Atlantic
    | (Batke et al., 2014)
Rhoa=1.204;
    Air density in (kg/m3)
SST=(28-3*(abs(yCenter)-10)/20)+1;   
    Sea surface temperature in (C), Holland (2008), Ts=SST-1, Ts=28-3*(abs(yCenter)-10)/20
VgToVCoeff=0.8;
    | Coefficient to convert gradient wind velocity to wind velocity at 10 m above surface as: 
    | V=VgToVCoeff*Vg, if VgToVCoeff=1, then V=Vg
inflowdirCalcMethod='bretschneider';
    | Inflow angle calculation method 
    | 'no': Inflow angle are not calculated, thetaVgrid=thetaVgridtan
    | 'bretschneider': Inflow angle are calculated based on Bretschneider (1972)
    | 'sobey': Inflow angle are calculated based on Sobey et al. (1977)
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
    Hurricane 1-min averaged wind velocity at 10 m above surface in x (East) direction on defined mesh in (m/s)
Vygrid
    Hurricane 1-min averaged wind velocity at 10 m above surface in y (North) direction on defined mesh in (m/s)
Vgrid
    Resultant hurricane 1-min averaged wind velocity at 10 m above surface (Vx^2+Vy^2)^0.5 on defined mesh in (m/s)
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
    |               and the circular isobars around the hurricane center (Boose et al., 2004)
thetaVgridtan
    | Angle (trigonometric direction) of hurricane velocity at each grid point in (Degree)
    | thetaVgridtan is tangential angle respect to radius. 
    | Note: Outputs has dimension of [M,N,L] where [M,N] is size of the x-y grid and [L] is number of time steps
    |        If flattendata='yes'; then Outputs has dimension of [M*L,N]
    |    Hurricane translation velocity needs to be added after rotational velocity is calculated 
    |        (e.g. Hu et al., 2012; Lin & Chavas, 2012)
    |    1-min averaged wind velocity needs to be converted to standard duration such as 
    |        10-min averaged wind by using a gust factor

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
    dPcdt=-1.16667; %Hurricane central pressure (Pc) intensity change in (hPa/hr)
    SST=27.77; %Sea surface temperature in (C)

    [Vxgrid,Vygrid,Vgrid,Pgrid,Rgrid,RVmax,thetagrid,thetaVgrid,thetaVgridtan]=hurricanewindh08(xgrid,ygrid,longCenter,latCenter,Pc,Vmax,Rknown,VgRknown,dPcdt,Vt,423e3,Pn,Rhoa,SST,...
        0.8,'bretschneider',VtAzmdir,'slosh','gc','no','no','quiver');

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

    %Hurricane central pressure (Pc) intensity change in (hPa/hr)
    dPcdt=[0.00000;-0.16667;0.00000;-0.16667;-0.50000;-0.50000;-0.50000;-0.50000;-1.00000;-0.66667;-0.16667;...
        0.66667;-1.33333;-1.83333;-1.50000;-1.50000;-1.33333;1.00000;-1.16667;-1.83333;-3.50000;-1.16667;...
        0.50000;1.33333;1.16667;0.50000;0.83333;3.33333;2.16667;2.83333;1.16667;0.83333;0.66667;...
        0.33333];

    Pn=101325; %Ambient surface pressure (external pressure) in (Pa)
    Rhoa=1.204; %Air density in (kg/m3)
    SST=27.77; %Sea surface temperature in (C)

    [Vxgrid,Vygrid,Vgrid,Pgrid,Rgrid,RVmax,thetagrid,thetaVgrid,thetaVgridtan]=hurricanewindh08(xgrid,ygrid,longtrack(4:27,1),lattrack(4:27,1),Pc(4:27,1),Vmax(4:27,1),Rknown(4:27,1),VgRknown(4:27,1),dPcdt(4:27,1),Vt(4:27,1),423e3,Pn,Rhoa,SST,...
        0.8,'bretschneider',VtAzmdir(4:27,1),'slosh','gc','no','no','quiver');

    %Converting 1-min sustained wind to 10-min averaged wind using gust factor
    %e.g. World Meteorological Organization (2015)
    Vxgrid=Vxgrid*0.88;
    Vygrid=Vygrid*0.88;
    Vgrid=Vgrid*0.88;


    %EXAMPLE 3

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

    dPcdt=3; %Hurricane central pressure (Pc) intensity change in (hPa/hr)
    SST=27.77; %Sea surface temperature in (C)

    [Vxgrid,Vygrid,Vgrid,Pgrid,Rgrid,RVmax,thetagrid,thetaVgrid,thetaVgridtan]=hurricanewindh08(xgrid,ygrid,longCenter,latCenter,Pc,Vmax,Rknown,VgRknown,dPcdt,Vt,423e3,Pn,Rhoa,SST,...
        0.8,'bretschneider',VtAzmdir,'slosh','gc','no','no','no');
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
%}

%--------------------------------------------------------------------------
%CODE
%--------------------------------------------------------------------------
%Assign default values

switch nargin
    case 10
        Rmax=423e3; Pn=101325; Rhoa=1.204; SST=(28-3.*(abs(yCenter)-10)./20)+1; VgToVCoeff=0.8; inflowdirCalcMethod='bretschneider'; VtAzmdir=0; backwindCalcMethod='no'; distCalcMethod='gc'; flattendata='no'; savedata='no'; dispout='no';
    case 11
        Pn=101325; Rhoa=1.204; SST=(28-3.*(abs(yCenter)-10)./20)+1; VgToVCoeff=0.8; inflowdirCalcMethod='bretschneider'; VtAzmdir=0; backwindCalcMethod='no'; distCalcMethod='gc'; flattendata='no'; savedata='no'; dispout='no';
    case 12
        Rhoa=1.204; SST=(28-3.*(abs(yCenter)-10)./20)+1; VgToVCoeff=0.8; inflowdirCalcMethod='bretschneider'; VtAzmdir=0; backwindCalcMethod='no'; distCalcMethod='gc'; flattendata='no'; savedata='no'; dispout='no';
    case 13
        SST=(28-3.*(abs(yCenter)-10)./20)+1; VgToVCoeff=0.8; inflowdirCalcMethod='bretschneider'; VtAzmdir=0; backwindCalcMethod='no'; distCalcMethod='gc'; flattendata='no'; savedata='no'; dispout='no';
    case 14
        VgToVCoeff=0.8; inflowdirCalcMethod='bretschneider'; VtAzmdir=0; backwindCalcMethod='no'; distCalcMethod='gc'; flattendata='no'; savedata='no'; dispout='no';
    case 15
        inflowdirCalcMethod='bretschneider'; VtAzmdir=0; backwindCalcMethod='no'; distCalcMethod='gc'; flattendata='no'; savedata='no'; dispout='no';
    case 16
        VtAzmdir=0; backwindCalcMethod='no'; distCalcMethod='gc'; flattendata='no'; savedata='no'; dispout='no';
    case 17
        backwindCalcMethod='no'; distCalcMethod='gc'; flattendata='no'; savedata='no'; dispout='no';
    case 18
        distCalcMethod='gc'; flattendata='no'; savedata='no'; dispout='no';
    case 19
        flattendata='no'; savedata='no'; dispout='no';
    case 20
        savedata='no'; dispout='no';
    case 21
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(xCenter)==1
    xCenter=xCenter';
end

if isrow(yCenter)==1
    yCenter=yCenter';
end

if isrow(Pc)==1
    Pc=Pc';
end

if isrow(Vmax)==1
    Vmax=Vmax';
end

if isrow(Rknown)==1
    Rknown=Rknown';
end

if isrow(VgRknown)==1
    VgRknown=VgRknown';
end

if isrow(dPcdt)==1
    dPcdt=dPcdt';
end

if isrow(Pn)==1
    Pn=Pn';
end

if isrow(SST)==1
    SST=SST';
end

if isrow(Vt)==1
    Vt=Vt';
end

if isrow(VtAzmdir)==1
    VtAzmdir=VtAzmdir';
end

%--------------------------------------------------------------------------
%Pre-assigning array

[Mgrid,Ngrid]=size(xgrid);
if length(xCenter(:,1))>1
    Rgrid=zeros(Mgrid,Ngrid,length(xCenter(:,1))); %Pre-assigning array
    thetagrid=zeros(Mgrid,Ngrid,length(xCenter(:,1))); %Pre-assigning array
    Pgrid=zeros(Mgrid,Ngrid,length(xCenter(:,1))); %Pre-assigning array
    Vgrid=zeros(Mgrid,Ngrid,length(xCenter(:,1))); %Pre-assigning array
    thetaVgridtan=zeros(Mgrid,Ngrid,length(xCenter(:,1))); %Pre-assigning array
    thetaVgrid=zeros(Mgrid,Ngrid,length(xCenter(:,1))); %Pre-assigning array
    Beta=zeros(Mgrid,Ngrid,length(xCenter(:,1))); %Pre-assigning array
else
    Rgrid=zeros(Mgrid,Ngrid); %Pre-assigning array
    thetagrid=zeros(Mgrid,Ngrid); %Pre-assigning array
    Pgrid=zeros(Mgrid,Ngrid); %Pre-assigning array
    Vgrid=zeros(Mgrid,Ngrid); %Pre-assigning array
    thetaVgridtan=zeros(Mgrid,Ngrid); %Pre-assigning array
    thetaVgrid=zeros(Mgrid,Ngrid); %Pre-assigning array
    Beta=zeros(Mgrid,Ngrid); %Pre-assigning array
end

RVmax=zeros(length(xCenter(:,1)),1); %Pre-assigning array
RVmaxs=zeros(length(xCenter(:,1)),1); %Pre-assigning array
VtR=zeros(Mgrid,Ngrid); %Pre-assigning array

%--------------------------------------------------------------------------
%Calculating distance (radius) from hurricane center to each point

for i=1:length(xCenter(:,1))

    %Calculating distance using cartesian formula
    if strcmp(distCalcMethod,'cart')==1

        distxy=sqrt((xgrid-xCenter(i,1)).^2+(ygrid-yCenter(i,1)).^2); %Calculating distance from (x,y) to (x(1),y(1))

        %Calculating angle of the line between start and end points

        thetarad=atan2(ygrid-yCenter(i,1),xgrid-xCenter(i,1)); %Angle in radian
        theta=rad2deg(thetarad); %Angle in degree

        %Add 360 to all numbers to have them all positive
        %Use mod(360) to take care of the ones larger than 360
        theta=mod((theta+360),360); 

        Rgrid(:,:,i)=distxy; %Distance from hurricane center to each point in (m)
        thetagrid(:,:,i)=theta; %Angle from hurricane center to each point in (degree)

    %Calculating distance using Vincenty formula
    elseif strcmp(distCalcMethod,'gc')==1

        %Converting to radian
        lat1rad=deg2rad(yCenter(i,1));
        lon1rad=deg2rad(xCenter(i,1));
        lat2rad=deg2rad(ygrid);
        lon2rad=deg2rad(xgrid);

        deltalatrad21=lat2rad-lat1rad;
        deltalonrad21=lon2rad-lon1rad;

        REarth=6371000; %Earth radius in (m), mean earth radius=6371000 m
        deltasigma=atan2(sqrt((cos(lat2rad).*sin(deltalonrad21)).^2+(cos(lat1rad).*sin(lat2rad)-sin(lat1rad).*cos(lat2rad).*cos(deltalonrad21)).^2),sin(lat1rad).*sin(lat2rad)+cos(lat1rad).*cos(lat2rad).*cos(deltalonrad21)); %Central angle
        arclen=REarth.*deltasigma; %Total distance of the line 
        distxy=arclen;

        %Calculating azimuth (bearing) between start and end of the line

        azimuthrad=atan2(sin(deltalonrad21).*cos(lat2rad),cos(lat1rad).*sin(lat2rad)-sin(lat1rad).*cos(lat2rad).*cos(deltalonrad21)); %Azimuth (bearing) in radian
        azimuth=rad2deg(azimuthrad); %Azimuth (bearing) in degree

        %Add 360 to all numbers to have them all positive
        %Use mod(360) to take care of the ones larger than 360
        azimuth=mod((azimuth+360),360); 

        Rgrid(:,:,i)=distxy; %Distance from hurricane center to each point in (m)
        thetagrid(:,:,i)=azimuth; %Angle from hurricane center to each point in (degree)

    end

end

%--------------------------------------------------------------------------
%Calculating hurricane wind and pressure using Holland (1980)

%Converting velocity 10 m above surface to gradient (mean boundary-layer) wind velocity
%e.g. Graham & Numm (1959), Young & Vinoth (2013), Phadke et al. (2003), Powell et al. (2003), Valamanesh et al. (2016), Wei et al. (2017)
%Shore Protection Manual (SPM), U.S. Army Corps of Engineers (1984)
%Coastal Engineering Manual (CEM), U.S. Army Corps of Engineers (2015)
Vgmax=Vmax./VgToVCoeff;

deltaP=(Pn-Pc); %Pressure drop from a defined external pressure to the cyclone center
%deltaP=deltaP.*0.01; %Convert from (Pa) to (hPa)

B=Vgmax.^2.*Rhoa.*exp(1)./(deltaP); %Hurricane shape parameter
B(B>2.5)=2.5; %B is typically between 1 and 2.5
B(B<1)=1; %B is typically between 1 and 2.5

%w=7.2921150e-5; %Rotational frequency of the earth (Radian/s)
%f=2.*w.*sin(deg2rad(yCenter)); %Coriolis parameter (Coriolis frequency)

for i=1:length(xCenter(:,1))

    %Initial guess for RVmax using Rankine vortex model
    %Rankine vortex model for R>=RVmax: V=Vmax*(RVmax/R)^x, where 0.4<x<0.6 (e,g. Phadke et al., 2003; Holland et al., 2010)
    %Rankine vortex model for R<RVmax: V=Vmax*(R/RVmax)^x, where 0.4<x<0.6 (e,g. Phadke et al., 2003; Holland et al., 2010)
    RVmaxini=Rknown(i,1)*(VgRknown(i,1)/Vmax(i,1))^2; %Initial guess

    %Calculating radius of maximum wind from known velocity and its radius
    fun = @(x)(VgRknown(i,1)-(B(i,1)/Rhoa*(x/Rknown(i,1))^B(i,1)*deltaP(i,1)*exp(-(x/Rknown(i,1))^B(i,1)))^0.5); %Hurricane radius for Vmax without coriolis parameter
    RVmax(i,1)=fzero(fun,RVmaxini); %Radius of maximum wind

    %Surface pressure at radius R
    Ps=Pc(i,1)+deltaP(i,1).*exp(-(RVmax(i,1)./Rgrid(:,:,i)).^B(i,1)); %Surface pressure at radius R from hurricane center
    
    Pgrid(:,:,i)=Ps; %Assigning hurricane pressure at each track point

end

%Surface pressure at radius RVmax
PsRVmax=Pc+deltaP.*exp(-1); %Surface pressure at radius R equal to RVmax from hurricane center

%--------------------------------------------------------------------------
%Calculating hurricane wind and pressure using Holland (2008)

deltaPshPa=deltaP.*0.01; %Convert from (Pa) to (hPa)

Ts=SST-1; %Ts is surface air temperature in (C) and SST is sea surface temperature in (C)
%Ts=28-3.*(yCenter-10)./20; %Surface air temperature in (C)
qs=0.9.*3.802./PsRVmax.*exp(17.67.*Ts./(243.5+Ts)); %Vapor pressure at an assumed relative humidity of 90%
Tvs=(Ts+273.15).*(1+0.81.*qs); %Virtual surface temperature in (K)
Rhoas=(PsRVmax./(286.9.*Tvs)); %Air density at surface (kg/m^3)

xs=0.6.*(1-deltaPshPa./215);
bs=-4.4e-5.*deltaPshPa.^2+0.01.*deltaPshPa+0.03.*dPcdt-0.014.*abs(yCenter)+0.15.*Vt.^xs+1; %Hurricane shape parameter at surface
Vmaxs=(bs.*deltaP./(Rhoas.*exp(1))).^0.5; %Hurricane maximum wind velocity at surface level (10m above surface)

for i=1:length(xCenter(:,1))

    %Initial guess for RVmax using Rankine vortex model
    %Rankine vortex model for R>=RVmax: V=Vmax*(RVmax/R)^x, where 0.4<x<0.6 (e,g. Phadke et al., 2003; Holland et al., 2010)
    %Rankine vortex model for R<RVmax: V=Vmax*(R/RVmax)^x, where 0.4<x<0.6 (e,g. Phadke et al., 2003; Holland et al., 2010)
    %RVmaxini=Rknown(i,1)*(VRknown(i,1)/Vmaxs(i,1))^2 %Initial guess

    %Calculating radius of maximum wind from known velocity and its radius
    %fun = @(x)(VRknown(i,1)-(bs(i,1)/Rhoas(i,1)*(x/Rknown(i,1))^bs(i,1)*deltaP(i,1)*exp(-(x/Rknown(i,1))^bs(i,1)))^xs(i,1)); %Hurricane radius for Vmax without coriolis parameter
    %RVmaxs(i,1)=fzero(fun,RVmaxini); %Radius of maximum wind
    RVmaxs(i,1)=RVmax(i,1);

    %Surface wind velocity at 10 m height at radius R
    %VR=(bs(i,1)./Rhoas(i,1).*(RVmaxs(i,1)./Rgrid(:,:,i)).^bs(i,1).*deltaP(i,1).*exp(-(RVmaxs(i,1)./Rgrid(:,:,i)).^bs(i,1))).^xs(i,1); %Gradient wind velocity at radius R without oriolis parameter
    VR=Vmaxs(i,1).*((RVmaxs(i,1)./Rgrid(:,:,i)).^bs(i,1).*exp(1-(RVmaxs(i,1)./Rgrid(:,:,i)).^bs(i,1))).^xs(i,1); %Gradient wind velocity at radius R without oriolis parameter

    VR(Rgrid(:,:,i)>=Rmax)=0; %e.g. Chavas & Emanuel (2010), Lin & Chavas (2012)
    VR(VR<0)=0;

    Vgrid(:,:,i)=VR; %Assigning hurricane surface velocity at each track point

end

%--------------------------------------------------------------------------
%Calculating hurricane velocity tangential and inflow angle and inflow velocity in x (East) and y (North) directions
%Tangential velocity has a right angle respect to radius
%thetaVgridtan is tangential angle respect to radius 
%thetaVgrid is inflow angle calculated based on Bretschneider (1972)
%Inflow angle: angle between the inwardly spiraling surface wind 
%and the circular isobars around the hurricane center (Boose et al., 2004)

for i=1:length(xCenter(:,1))

    %Calculating inflow angle based on Bretschneider (1972), e.g Phadke et al.,(2003), 
    if strcmp(inflowdirCalcMethod,'bretschneider')==1
        %Jeong, C. K., Panchang, V., & Demirbilek, Z. (2012). 
        %Parametric adjustments to the rankine vortex wind model for Gulf of Mexico hurricanes. 
        %Journal of Offshore Mechanics and Arctic Engineering, 134(4), 041102.
        
        %Beta is an inflow angle inward from the tangential flow (Martino et al., 2001)
        Beta(Rgrid<RVmax(i,1))=10+(1+Rgrid(Rgrid<RVmax(i,1))./RVmax(i,1));
        Beta(Rgrid>=RVmax(i,1) & Rgrid<1.2*RVmax(i,1))=20+25.*(Rgrid(Rgrid>=RVmax(i,1) & Rgrid<1.2*RVmax(i,1))./RVmax(i,1)-1);
        Beta(Rgrid>=1.2*RVmax(i,1))=25;

    %Calculating inflow angle based on Sobey et al. (1977), e.g. Martino et al. (2001)
    elseif strcmp(inflowdirCalcMethod,'sobey')==1
        %Sobey, R. J., Harper, B. A., & Stark, K. P. (1977). 
        %Numerical simulation of tropical cyclone storm surge. 
        %James Cook University of North Queensland, Department of Civil & Systems Engineering.

        %Martino, C. D., Cheung, K. F., Phadke, A. C., & Houston, S. H. (2001, January). 
        %Modeling of hurricane waves in Hawaiian waters. 
        %In The Eleventh International Offshore and Polar Engineering Conference. International Society of Offshore and Polar Engineers.    
        
        %Beta is an inflow angle inward from the tangential flow (Martino et al., 2001)
        Beta(Rgrid<RVmax(i,1))=10.*(Rgrid(Rgrid<RVmax(i,1))./RVmax(i,1));
        Beta(Rgrid>=RVmax(i,1) & Rgrid<1.2*RVmax(i,1))=10+75.*(Rgrid(Rgrid>=RVmax(i,1) & Rgrid<1.2*RVmax(i,1))./RVmax(i,1)-1);
        Beta(Rgrid>=1.2*RVmax(i,1))=25;

    else
        Beta=0;
    
    end

end

%Calculating velocity vector direction using cartesian formula
if strcmp(distCalcMethod,'cart')==1

    %Hurricane is in northern hemisphere
    if mean(yCenter(:,1))>=0
        %In northern hemisphere, velocity vector rotate counter clockwise with right angle respect to radius
        thetaVgridtan=thetagrid+90; %Angle of velocity vector in degree
        thetaVgrid=thetaVgridtan+Beta; %Inflow angle of velocity vector in degree

    %Hurricane is in southern hemisphere
    elseif mean(yCenter(:,1))<0
        %In southern hemisphere, velocity vector rotate clockwise with right angle respect to radius
        thetaVgridtan=thetagrid-90; %Angle of velocity vector in degree
        thetaVgrid=thetaVgridtan-Beta; %Inflow angle of velocity vector in degree

    end

    %Add 360 to all numbers to have them all positive
    %Use mod(360) to take care of the ones larger than 360
    thetaVgridtan=mod((thetaVgridtan+360),360); 
    thetaVgrid=mod((thetaVgrid+360),360); 

%Calculating velocity vector direction using great circle formula
elseif strcmp(distCalcMethod,'gc')==1

    %Converting azimuth (bearing) to trigonometric direction
    thetagrid=-thetagrid+90;

    %Hurricane is in northern hemisphere
    if mean(yCenter(:,1))>=0
        %In northern hemisphere, velocity vector rotate counter clockwise with right angle respect to radius
        thetaVgridtan=thetagrid+90; %Angle of velocity vector in degree
        thetaVgrid=thetaVgridtan+Beta; %Inflow angle of velocity vector in degree

    %Hurricane is in southern hemisphere
    elseif mean(yCenter(:,1))<0
        %In southern hemisphere, velocity vector rotate clockwise with right angle respect to radius
        thetaVgridtan=thetagrid-90; %Angle of velocity vector in degree
        thetaVgrid=thetaVgridtan-Beta; %Inflow angle of velocity vector in degree

    end

    %Add 360 to all numbers to have them all positive
    %Use mod(360) to take care of the ones larger than 360
    thetaVgridtan=mod((thetaVgridtan+360),360); 
    thetaVgrid=mod((thetaVgrid+360),360); 

end

Vxgrid=Vgrid.*cos(deg2rad(thetaVgrid)); %Hurricane surafec velocity in x (East) direction
Vygrid=Vgrid.*sin(deg2rad(thetaVgrid)); %Hurricane surface velocity in y (North) direction

%--------------------------------------------------------------------------
%Adding background wind velocity due to hurricane motion to hurricane wind velocity

%Converting azimuth (bearing) to trigonometric direction
VtTrigdir=-VtAzmdir+90;

%Add 360 to all numbers to have them all positive
%Use mod(360) to take care of the ones larger than 360
VtTrigdir=mod((VtTrigdir+360),360); 

%Adding background wind velocity as constant value to all hurricane velocities
if strcmp(backwindCalcMethod,'constant')==1

    for i=1:length(xCenter(:,1))
        VtR(:)=Vt(i,1); %Translational wind velocity at radius R
        VtR(Rgrid(:,:,i)>=Rmax)=0; %e.g. Chavas & Emanuel (2010), Lin & Chavas (2012)

        Vxgrid(:,:,i)=Vxgrid(:,:,i)+VtR.*cos(deg2rad(VtTrigdir(i,1))); %Hurricane velocity in x (East) direction
        Vygrid(:,:,i)=Vygrid(:,:,i)+VtR.*sin(deg2rad(VtTrigdir(i,1))); %Hurricane velocity in y (North) direction
    end

    Vgrid=sqrt(Vxgrid.^2+Vygrid.^2); %Assigning hurricane velocity at each track point

%Adding background wind velocity using SLOSH model by Jelesnianski et al. (1992)
elseif strcmp(backwindCalcMethod,'slosh')==1

    for i=1:length(xCenter(:,1))
        VtR=Vt(i,1).*(RVmax(i,1).*Rgrid(:,:,i))./((RVmax(i,1))^2+(Rgrid(:,:,i)).^2); %Translational wind velocity at radius R
        VtR(VtR<0)=0;
        VtR(Rgrid(:,:,i)>=Rmax)=0; %e.g. Chavas & Emanuel (2010), Lin & Chavas (2012)

        Vxgrid(:,:,i)=Vxgrid(:,:,i)+VtR.*cos(deg2rad(VtTrigdir(i,1))); %Hurricane velocity in x (East) direction
        Vygrid(:,:,i)=Vygrid(:,:,i)+VtR.*sin(deg2rad(VtTrigdir(i,1))); %Hurricane velocity in y (North) direction
    end

    Vgrid=sqrt(Vxgrid.^2+Vygrid.^2); %Assigning hurricane velocity at each track point

%Adding background wind velocity using Lin & Chavas (2012)
elseif strcmp(backwindCalcMethod,'lin')==1

    for i=1:length(xCenter(:,1))
        %Hurricane is in northern hemisphere
        if mean(yCenter(:,1))>=0
            %In northern hemisphere, velocity vector rotate counter clockwise with right angle respect to radius

            VtR(:)=0.55*Vt(i,1); %Translational wind velocity at radius R
            VtR(Rgrid(:,:,i)>=Rmax)=0; %e.g. Chavas & Emanuel (2010), Lin & Chavas (2012)

            Vxgrid(:,:,i)=Vxgrid(:,:,i)+VtR.*cos(deg2rad(VtTrigdir(i,1)+20)); %Hurricane velocity in x (East) direction
            Vygrid(:,:,i)=Vygrid(:,:,i)+VtR.*sin(deg2rad(VtTrigdir(i,1)+20)); %Hurricane velocity in y (North) direction

        %Hurricane is in southern hemisphere
        elseif mean(yCenter(:,1))<0
            %In southern hemisphere, velocity vector rotate clockwise with right angle respect to radius

            VtR(:)=0.55*Vt(i,1); %Translational wind velocity at radius R
            VtR(Rgrid(:,:,i)>=Rmax)=0; %e.g. Chavas & Emanuel (2010), Lin & Chavas (2012)

            Vxgrid(:,:,i)=Vxgrid(:,:,i)+VtR.*cos(deg2rad(VtTrigdir(i,1)-20)); %Hurricane velocity in x (East) direction
            Vygrid(:,:,i)=Vygrid(:,:,i)+VtR.*sin(deg2rad(VtTrigdir(i,1)-20)); %Hurricane velocity in y (North) direction

        end
    end

    Vgrid=sqrt(Vxgrid.^2+Vygrid.^2); %Assigning hurricane velocity at each track point

end

%--------------------------------------------------------------------------
%Flating results

if strcmp(flattendata,'yes')==1

    [M,N,L]=size(Vgrid);
    Vxgrid=reshape(Vxgrid,M*L,N);
    Vygrid=reshape(Vygrid,M*L,N);
    Vgrid=reshape(Vgrid,M*L,N);
    Pgrid=reshape(Pgrid,M*L,N);

    dispout='no'; %Results of flatten array cannot be plotted

end

%--------------------------------------------------------------------------
%Saving data

if strcmp(savedata,'yes')==1

    if strcmp(flattendata,'yes')==1

        dlmwrite('Vx.dat',Vxgrid,'delimiter',' ')
        dlmwrite('Vy.dat',Vygrid,'delimiter',' ')
        dlmwrite('V.dat',Vgrid,'delimiter',' ')
        dlmwrite('P.dat',Pgrid,'delimiter',' ')

    %Flating results if they were not flatten before
    elseif strcmp(flattendata,'no')==1

        [M,N,L]=size(Vgrid);
        Vxgridrshp=reshape(Vxgrid,M*L,N);
        Vygridrshp=reshape(Vygrid,M*L,N);
        Vgridrshp=reshape(Vgrid,M*L,N);
        Pgridrshp=reshape(Pgrid,M*L,N);
 
        dlmwrite('Vx.dat',Vxgridrshp,'delimiter',' ')
        dlmwrite('Vy.dat',Vygridrshp,'delimiter',' ')
        dlmwrite('V.dat',Vgridrshp,'delimiter',' ')
        dlmwrite('P.dat',Pgridrshp,'delimiter',' ')

    end

end

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'imagesc')==1 %2 dimensional plot using imagesc or imshow
    
    subplot(1,2,1)
    imagesc(xgrid,ygrid,Vgrid(:,:,end))
    set(gca,'YDir','normal'); %Correct y axis direction

    subplot(1,2,2)
    imagesc(xgrid,ygrid,Pgrid(:,:,end))
    set(gca,'YDir','normal'); %Correct y axis direction

elseif strcmp(dispout,'pcolor')==1 %2 dimensional plot using pcolor

    subplot(1,2,1)
    pcol=pcolor(xgrid,ygrid,Vgrid(:,:,end));
    set(pcol,'edgecolor','none')

    subplot(1,2,2)
    pcol=pcolor(xgrid,ygrid,Pgrid(:,:,end));
    set(pcol,'edgecolor','none')

elseif strcmp(dispout,'contour')==1 %2 dimensional contour plot

    subplot(1,2,1)
    [Cont,hCont]=contourf(xgrid,ygrid,Vgrid(:,:,end));
    set(hCont,'LineColor','none');

    subplot(1,2,2)
    [Cont,hCont]=contourf(xgrid,ygrid,Pgrid(:,:,end));
    set(hCont,'LineColor','none');

elseif strcmp(dispout,'quiver')==1 %Surface plot
    
    subplot(1,2,1)
    [Cont,hCont]=contour(xgrid,ygrid,Vgrid(:,:,end));
    hold on
    qvr1=quiver(xgrid,ygrid,Vxgrid(:,:,end),Vygrid(:,:,end));
    %set(qvr1,'FaceColor','interp','EdgeColor','none')

    subplot(1,2,2)
    imagesc(xgrid,ygrid,Pgrid(:,:,end))
    set(gca,'YDir','normal'); %Correct y axis direction

end

%Setting plot properties
if (strcmp(dispout,'imagesc')==1 | strcmp(dispout,'pcolor')==1 | strcmp(dispout,'contour')==1 | strcmp(dispout,'quiver')==1)

    subplot(1,2,1)

    %Setting axis limits
    %xlim([min(xgrid(:)),max(xgrid(:))])
    %ylim([min(ygrid(:)),max(ygrid(:))])

    %Plotting colorbar
    colorbar;

    %Setting label
    xlabel('x')
    ylabel('y')
    title('Velocity (m/s)')
    
    subplot(1,2,2)

    %Setting axis limits
    %xlim([min(xgrid(:)),max(xgrid(:))])
    %ylim([min(ygrid(:)),max(ygrid(:))])

    %Plotting colorbar
    colorbar;

    %Setting label
    xlabel('x')
    ylabel('y')
    title('Pressure (Pa)')
    
end

%--------------------------------------------------------------------------
