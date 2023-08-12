.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 202-12-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

scientimate.windfetch
=====================

.. code:: python

    wind_fetch, z_mean, z_point, bed_slope_mean, x_fetch, y_fetch, z_fetch, dist_along_fetch = scientimate.windfetch(xgrid, ygrid, zgrid, x_point, y_point, winddir=0, waterlevel=0, shorelinelevel=0, n_midpoints=100, distCalcMethod='gc', interpMethod='nearest', dispout='no')

Description
-----------

Calculate a wind fecth and z (elevation) profile along a path over water for a given 2d x-y domain (map, image, ...)

Inputs
------

xgrid
    x (longitude) data as a [M*N] array
ygrid
    y (latitude) data as a [M*N] array
zgrid
    | z (elevation) data as a [M*N] array
    | z>=0 is considered land, z<0 is considered water
x_point
    | x (longitude) of the points to find fetches toward those points
    | If x data are longitude, they should be in Degree
y_point
    | y (latitude) of the points to find fetches toward those points
    | If y data are latitude, they should be in Degree
winddir=0
    | Meteorological wind direction in (Degree)
    | It represents a direction wind comes from and is measured counter-clockwise from the North
    | 0 (degree): from North, 90 (degree): from East, 180 (degree): from South, 270 (degree): from West
waterlevel=0
    | Water surface level (Should have the same datum as zgrid)
    | Size of waterlevel should be either 1 or should be equal to size of winddir 
    | If water level is zero (waterlevel=0), it means that the water level is at a level of surface elevation (z=0) in bathymetry data  
    | If water surface level is positive (waterlevel>0), it means that the water level is above zero elevation (z=0) in bathymetry data  
    |     ,therefore water level will be subtracted from bathymetry elevation 
    | If water surface level is negative (waterlevel<0), it means that the water level is below zero elevation (z=0) in bathymetry data  
    |     ,therefore water level will be added to bathymetry elevation 
shorelinelevel=0
    | Shoreline elevation, 
    | z>=shorelinelevel is considered land, z<shorelinelevel is considered water
n_midpoints=100
    | Number of middle points generated between first point and last point
    | if n_midpoints=1 then 1 point between first point and last point is generated
distCalcMethod='gc'
    | Distance calculation method 
    | 'cart': Distances are calculated on cartesian coordinate
    | 'gc': Distances are calculated on Great Circle based on Vincenty formula, Vincenty (1975)
    | Earth radius considered as mean earth radius=6371000 m
interpMethod='nearest'
    | Interpolation method 
    | 'linear': Use default or 'linear' method to interpolate
    | 'nearest': Use nearest neighbor method to interpolate
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

wind_fetch
    | Total wind fecth distance along a path from shoreline to point (x_point,y_point) with angle of winddir
    | If input data are latitude and longitude in Degree, wind_fetch is in m
z_mean
    Mean of z (elevation) along a wind fetch
z_point
    Value of z (elevation) at (x_point,y_point) location
bed_slope_mean
    Mean of bed slope along a wind fetch
x_fetch
    x (longitude) along a path at given points (x_fetch,y_fetch)
y_fetch
    x (latitude) along a path at given points (x_fetch,y_fetch)
z_fetch
    z (elevation, ...) data along a path at given points (x_fetch,y_fetch)
dist_along_fetch
    | Distance at each points of (x_fetch,y_fetch) from the domain boundary
    | If input data are latitude and longitude in Degree, dist_along_fetch is in m

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np
    import os

    xgrid,ygrid=np.meshgrid(np.linspace(-5,5,100),np.linspace(-5,5,100))
    R=np.sqrt(xgrid**2+ygrid**2)
    zgrid=np.sin(R)/R
    x_point=-4
    y_point=-2
    winddir=90
    waterlevel=0
    shorelinelevel=0
    n_midpoints=100
    wind_fetch,z_mean,z_point,bed_slope_mean,x_fetch,y_fetch,z_fetch,dist_along_fetch=sm.windfetch(xgrid,ygrid,zgrid,x_point,y_point,winddir,waterlevel,shorelinelevel,n_midpoints,'cart','nearest','yes')

    xgrid,ygrid=np.meshgrid(np.linspace(-5,5,100),np.linspace(-5,5,100))
    zgrid=ygrid*np.sin(xgrid)-xgrid*np.cos(ygrid)-3
    x_point=0
    y_point=0
    winddir=np.arange(0,360+30,30)
    waterlevel=np.ones(np.size(winddir))*0.1
    shorelinelevel=0
    n_midpoints=100
    wind_fetch,z_mean,z_point,bed_slope_mean,x_fetch,y_fetch,z_fetch,dist_along_fetch=sm.windfetch(xgrid,ygrid,zgrid,x_point,y_point,winddir,waterlevel,shorelinelevel,n_midpoints,'cart','nearest','yes')

    #Download Persina Gulf and Gulf of Oman with coordinate exteneds xmin=47, xmax=63, ymin=19, ymax=31 from https://maps.ngdc.noaa.gov/viewers/grid-extract/index.html
    xyzfilename='xyz.csv'; #e.g. xyzfilename='PersianGulf.csv'
    xyzfilelocation=os.getcwd(); #e.g. xyzfilelocation=r'C:/datafolder'
    x,y,z=sm.readxyz(xyzfilename,xyzfilelocation,1,'all');
    xgrid,ygrid,zgrid=sm.interpxyz2grid(x,y,z,100,'points',np.nanmin(x),np.nanmax(x),np.nanmin(y),np.nanmax(y),np.nanmin(z),np.nanmax(z),'all','nearest','no');
    x_point=58.0 #Or x_point=52.0
    y_point=24.5 #Or y_point=26.0
    winddir=[0,40,135,280]
    waterlevel=0
    shorelinelevel=0
    n_midpoints=100
    wind_fetch,z_mean,z_point,bed_slope_mean,x_fetch,y_fetch,z_fetch,dist_along_fetch=sm.windfetch(xgrid,ygrid,zgrid,x_point,y_point,winddir,waterlevel,shorelinelevel,n_midpoints,'gc','nearest','yes');

References
----------

Vincenty, T. (1975). 
Direct and inverse solutions of geodesics on the ellipsoid with application of nested equations. 
Survey review, 23(176), 88-93.

* http://www.movable-type.co.uk/scripts/latlong.html
* http://edwilliams.org/gccalc.htm
* http://edwilliams.org/avform.htm
* https://www.nhc.noaa.gov/gccalc.shtml

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
