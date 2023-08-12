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

hurricanewindinflowangle
========================

.. code:: MATLAB

    [Vxgrid, Vygrid, thetaVgrid, thetaVgridtan, thetagrid, Rgrid] = hurricanewindinflowangle(xgrid, ygrid, Vgrid, xCenter, yCenter, RVmax, inflowdirCalcMethod, distCalcMethod, dispout)

Description
-----------

Calculate hurricane velocity tangential and inflow angle and inflow velocity in x (East) and y (North) directions

Inputs
------

xgrid
    | x (longitude) of points which outputs are calculated at as a [M*N] array 
    | xgrid can be a single point or 1d or 2d array 
ygrid
    | y (latitude) of points which outputs are calculated at as a [M*N] array 
    | ygrid can be a single point or 1d or 2d array
Vgrid
    | Resultant hurricane wind velocity (Vx^2+Vy^2)^0.5 on (xgrid,ygrid) as a [M*N*L] array in (m/s)
    | L is a number of time steps
    | If only angle values are required, then set Vgrid equal to an arbitary constant such as Vgrid=1
    | For demonstaration purpose, set Vgrid equal to an arbitary constant such as Vgrid=1
xCenter
    x (longitude) of hurricane center (track) as a [L] array
yCenter
    y (latitude) of hurricane center (track) as a [L] array
RVmax
    Distance (radius) from hurricane center to a location of maximum hurricane wind velocity as a [L] array in (m)
inflowdirCalcMethod='bretschneider';
    | Inflow angle calculation method 
    | 'no': Inflow angle are not calculated, thetaVgrid=thetaVgridtan
    | 'bretschneider': Inflow angle are calculated based on Bretschneider (1972)
    | 'sobey': Inflow angle are calculated based on Sobey et al. (1977)
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

Vxgrid
    Hurricane wind velocity after applying inflow angle in x (East) direction on defined mesh in (m/s)
Vygrid
    Hurricane wind velocity after applying inflow angle in y (North) direction on defined mesh in (m/s)
thetaVgrid
    | Inflow angle (trigonometric direction) of hurricane velocity at each grid point in (Degree)
    | Inflow angle: angle between the inwardly spiraling surface wind 
    |     and the circular isobars around the hurricane center (Boose et al., 2004)
thetaVgridtan
    | Angle (trigonometric direction) of hurricane velocity at each grid point in (Degree)
    | thetaVgridtan is tangential angle respect to radius. 
    | Note: Outputs has dimension of [M,N,L] where [M,N] is size of the x-y grid and [L] is number of time steps
thetagrid
    Angle from hurricane center to each point on the grid in (Degree)
Rgrid
    Distance (radius) from hurricane center to each point on the grid

Examples
--------

.. code:: MATLAB

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

    %Calculating distance using spherical law of cosines
    Rgrid=(acos(sin(deg2rad(latCenter)).*sin(deg2rad(ygrid))+cos(deg2rad(latCenter)).*cos(deg2rad(ygrid)).*cos(deg2rad(xgrid)-deg2rad(longCenter)))).*6371000; %Radius

    %Calculating hurricane velocity at each radius using SLOSH model
    RVmax=32197; %Radius from hurricane center to a location of maximum hurricane wind
    Vgrid=Vgmax.*(2.*RVmax.*Rgrid)./((RVmax)^2+(Rgrid).^2); %Hurricane wind velocity at radius R

    [Vxgrid,Vygrid,thetaVgrid,thetaVgridtan,thetagrid,Rgrid]=hurricanewindinflowangle(xgrid,ygrid,Vgrid,longCenter,latCenter,RVmax,'bretschneider','gc','quiver');

References
----------

Data

* www.nhc.noaa.gov/data/
* www.nhc.noaa.gov/data/hurdat/hurdat2-format-nencpac.pdf
* coast.noaa.gov/hurricanes
* www.aoml.noaa.gov/hrd/data_sub/re_anal.html

Boose, E. R., Serrano, M. I., & Foster, D. R. (2004). 
Landscape and regional impacts of hurricanes in Puerto Rico. 
Ecological Monographs, 74(2), 335-352.

Bretschneider, C. L. (1972, January). 
A non-dimensional stationary hurricane wave model. 
In Offshore Technology Conference. Offshore Technology Conference.

Phadke, A. C., Martino, C. D., Cheung, K. F., & Houston, S. H. (2003). 
Modeling of tropical cyclone winds and waves for emergency management. 
Ocean Engineering, 30(4), 553-578.

Sobey, R. J., Harper, B. A., & Stark, K. P. (1977). 
Numerical simulation of tropical cyclone storm surge. 
James Cook University of North Queensland, Department of Civil & Systems Engineering.

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
