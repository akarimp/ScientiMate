.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-11-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

hurricanebackgroundwind
=======================

.. code:: MATLAB

    [Vxgridbackwind, Vygridbackwind, Vgridbackwind, Rgrid] = hurricanebackgroundwind(xgrid, ygrid, Vxgrid, Vygrid, xCenter, yCenter, Vt, VtAzmdir, RVmax, backwindCalcMethod, Rmax, distCalcMethod, dispout)

Description
-----------

Calculate and add background wind velocity due to hurricane front motion to hurricane rotational wind velocity

Inputs
------

xgrid
    | x (longitude) of points which outputs are calculated at as a [M*N] array 
    | xgrid can be a single point or 1d or 2d array 
ygrid
    | y (latitude) of points which outputs are calculated at as a [M*N] array 
    | ygrid can be a single point or 1d or 2d array
Vxgrid
    | Hurricane wind velocity in x (East) direction as a [M*N*L] array in (m/s)
    | L is a number of time steps
    | If only background wind is required, set Vxgrid=0 and Vygrid=0
Vygrid
    | Hurricane wind velocity in y (North) direction as a [M*N*L] array in (m/s)
    | L is a number of time steps
    | If only background wind is required, set Vxgrid=0 and Vygrid=0
xCenter
    x (longitude) of hurricane center (track) as a [L] array
yCenter
    y (latitude) of hurricane center (track) as a [L] array
Vt
    Hurricane central translational velocity as a [L] array in (m/s)
VtAzmdir
    | Hurricane center velocity azimuth (bearing) direction as a [L] array in (Degree)
    | azimuth (bearing) direction which is measured clockwise from the north:
    | 0 (degree): toward North, 90 (degree): toward East, 180 (degree): toward South, 270 (degree): toward West 
RVmax
    Distance (radius) from hurricane center to a location of maximum hurricane wind velocity as a [L] array in (m)
backwindCalcMethod='slosh';
    | Calculation method for adding background wind velocity due to hurricane motion
    | background wind velocity is added to points with R<=Rmax, e.g. Chavas & Emanuel (2010), Lin & Chavas (2012)
    | 'no': background wind velocity is not added to hurricane velocities
    | 'constant': background wind velocity is added as constant value to hurricane velocities
    | 'slosh': Background wind velocity is calculated and added using SLOSH model by Jelesnianski et al. (1992)
    | 'lin': Background wind velocity is calculated and added based on Lin & Chavas (2012)
Rmax=423e3;
    | Maximum radius of hurricane from hurricane center in (m)
    | Outputs for points with R>Rmax is set to zero
    | Median values of Rmax is 423 km (e.g. Chavas & Emanuel 2010; Lin & Chavas, 2012)
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
    | 'quiver': 2 dimensional vector plot 
    | 'no': not display 
    | Use dispout='no' if calculation mesh is not 2d array
    | if there is more than one time step, only the last one is plotted
    | if flattendata='yes'; then dispout is set as dispout='no';

Outputs
-------

Vxgridbackwind
    Hurricane wind velocity with background wind in x (East) direction on defined mesh in (m/s)
Vygridbackwind
    Hurricane wind velocity with background wind in y (North) direction on defined mesh in (m/s)
Vgridbackwind
    | Resultant hurricane wind velocity (Vx^2+Vy^2)^0.5 with background wind on defined mesh in (m/s)
    | Gradient wind velocity converted to wind velocity at 10 m above surface by V=VgToVCoeff*Vg
Rgrid
    | Distance (radius) from hurricane center to each point on the grid
    | Note: Outputs has dimension of [M,N,L] where [M,N] is size of the x-y grid and [L] is number of time steps

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

    %Hurricane Katrina translational velocity (m/s) at max velocity
    Vt=5.18467;

    %Hurricane Katrina velocity azimuth (bearing) in (Degree) at max velocity
    VtAzmdir=306.76219;

    %Hurricane Katrina 1-min sustained maximum velocity (m/s) at max velocity
    Vmax=76.5;
    Vmax=Vmax-Vt; %Removing hurricane translation velocity from Vmax
    Vgmax=Vmax/0.8; %Converting surface velocity to gradient velocity

    %Calculating distance and angle using spherical law of cosines
    Rgrid=(acos(sin(deg2rad(latCenter)).*sin(deg2rad(ygrid))+cos(deg2rad(latCenter)).*cos(deg2rad(ygrid)).*cos(deg2rad(xgrid)-deg2rad(longCenter)))).*6371000; %Radius
    thetagrid=atan2(sin(deg2rad(xgrid)-deg2rad(longCenter)).*cos(deg2rad(ygrid)),cos(deg2rad(latCenter)).*sin(deg2rad(ygrid))-sin(deg2rad(latCenter)).*cos(deg2rad(ygrid)).*cos(deg2rad(xgrid)-deg2rad(longCenter))); %Azimuth in radian
    thetagrid=-thetagrid+pi/2; %Converting azimuth to trigonometric direction
    thetagrid=thetagrid+pi/2; %Angle of velocity vector in degree (right angle respect to radius)

    %Calculating hurricane velocity at each radius using SLOSH model
    RVmax=32197; %Radius from hurricane center to a location of maximum hurricane wind
    Vgrid=Vgmax.*(2.*RVmax.*Rgrid)./((RVmax)^2+(Rgrid).^2); %Hurricane wind velocity at radius R
    Vxgrid=Vgrid.*cos(thetagrid); %Hurricane velocity in x (East) direction
    Vygrid=Vgrid.*sin(thetagrid); %Hurricane velocity in y (North) direction

    [Vxgridbackwind,Vygridbackwind,Vgridbackwind,Rgrid]=hurricanebackgroundwind(xgrid,ygrid,Vxgrid,Vygrid,longCenter,latCenter,Vt,VtAzmdir,RVmax,'slosh',423e3,'gc','quiver');


    %EXAMPLE 2

    %Creating calculation mesh
    [xgrid,ygrid]=meshgrid(linspace(-98,-68,100),linspace(16,44,100));

    %Longitude of Hurricane Katrine center at max velocity
    longCenter=-88.6;

    %Latitude of Hurricane Katrine center at max velocity
    latCenter=26.3;

    %Hurricane Katrina translational velocity (m/s) at max velocity
    Vt=5.18467;

    %Hurricane Katrina velocity azimuth (bearing) in (Degree) at max velocity
    VtAzmdir=306.76219;

    RVmax=32197; %Radius from hurricane center to a location of maximum hurricane wind

    Vxgrid=0; %Hurricane velocity in x (East) direction
    Vygrid=0; %Hurricane velocity in y (North) direction

    [Vxgridbackwind,Vygridbackwind,Vgridbackwind,Rgrid]=hurricanebackgroundwind(xgrid,ygrid,Vxgrid,Vygrid,longCenter,latCenter,Vt,VtAzmdir,RVmax,'slosh',423e3,'gc','quiver');


References
----------

Data

* www.nhc.noaa.gov/data/
* www.nhc.noaa.gov/data/hurdat/hurdat2-format-nencpac.pdf
* coast.noaa.gov/hurricanes
* www.aoml.noaa.gov/hrd/data_sub/re_anal.html

Chavas, D. R., & Emanuel, K. A. (2010). 
A QuikSCAT climatology of tropical cyclone size. 
Geophysical Research Letters, 37(18).

Jelesnianski, C. P., Chen, J., & Shaffer, W. A. (1992). 
SLOSH: Sea, lake, and overland surges from hurricanes (Vol. 48). 
US Department of Commerce, National Oceanic and Atmospheric Administration, National Weather Service.

Lin, N., & Chavas, D. (2012). 
On hurricane parametric wind and applications in storm surge modeling. 
Journal of Geophysical Research: Atmospheres, 117(D9).

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
