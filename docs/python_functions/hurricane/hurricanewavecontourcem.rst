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

scientimate.hurricanewavecontourcem
===================================

.. code:: python

    Hgrid = scientimate.hurricanewavecontourcem(xgrid, ygrid, xCenter, yCenter, Hmax, RVmax, VtAzmdir=0, CalcMethod='cem', distCalcMethod='gc', dispout='no')

Description
-----------

| Calculate hurricane wave height field (contours) on given mesh using 
| method from Shore Protection Manual (SPM), U.S. Army Corps of Engineers (1984), and Young (1988)
| For method from Young (1988) use hurricanewavecontoury88
| For method from Hwang (2016) and Hwang & Walsh (2016) use hurricanewavecontourh16

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
Hmax
    Hurricane maximum wave height
RVmax
    Distance (radius) from hurricane center to a location of maximum hurricane wind velocity (m)
VtAzmdir=0
    | Hurricane center velocity azimuth (bearing) direction in (Degree)
    | azimuth (bearing) direction which is measured clockwise from the north:
    | 0 (degree): toward North, 90 (degree): toward East, 180 (degree): toward South, 270 (degree): toward West 
CalcMethod='cem'
    | Hurricane wave contour calculation method 
    | 'spm': Use method by Shore Protection Manual (SPM),
    |     U.S. Army Corps of Engineers (1984) in deep water
    | 'cem': Use method by Coastal Engineering Manual (CEM),
    |     U.S. Army Corps of Engineers (2015)
    |     For hurricane with maximum velocity of Vmax>40 (m/s) and translational velocity of Vt<12 (m/s)
    | 'cemnext': Use method by Coastal Engineering Manual (CEM),
    |     U.S. Army Corps of Engineers (Next Release)
    |     For hurricane with maximum velocity of Vmax=40 (m/s) and translational velocity of Vt=7.5 (m/s)
    | 'youngvt2': Use method by Young (1988)
    |     For hurricane with maximum velocity of Vmax=40 (m/s) and translational velocity of Vt=2.5 (m/s)
    | 'youngvt5': Use method by Young (1988)
    |     For hurricane with maximum velocity of Vmax=40 (m/s) and translational velocity of Vt=5 (m/s)
    | 'youngvt7': Use method by Young (1988)
    |     For hurricane with maximum velocity of Vmax=40 (m/s) and translational velocity of Vt=7.5 (m/s)
    | 'youngvt10': Use method by Young (1988)
    |     For hurricane with maximum velocity of Vmax=40 (m/s) and translational velocity of Vt=10 (m/s)
distCalcMethod='gc'
    | Distance calculation method 
    | 'cart': Distances are calculated on cartesian coordinate
    | 'gc': Distances are calculated on Great Circle based on Vincenty formula, Vincenty (1975)
    | Earth radius coonsidered as mean earth radius=6371000 m
dispout='no'
    | Define to display outputs or not
    | 'imagesc': 2 dimensional plot using imagesc or imshow
    | 'pcolor': 2 dimensional plot using pcolor
    | 'contour': 2 dimensional contour plot, number of contour=ncolor
    | 'no': not display 
    | Use dispout='no' if calculation mesh is not 2d array
    | if there is more than one time step, only the last one is plotted

Outputs
-------

Hgrid
    Hurricane wave height on grid mesh

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np


    #EXAMPLE 1

    #Creating calculation mesh
    xgrid,ygrid=np.meshgrid(np.linspace(-98,-68,100),np.linspace(16,44,100))

    #Longitude of Hurricane Katrine center at max velocity
    longCenter=-88.6

    #Latitude of Hurricane Katrine center at max velocity
    latCenter=26.3

    #Hurricane Katrina maximum significant wave height (m) at max velocity
    Hmax=24.9821

    #Hurricane Katrina radius from hurricane center to a location of maximum hurricane wind velocity (m) at max velocity
    RVmax=6.2750e+004

    #Hurricane Katrina velocity azimuth (bearing) in (Degree) at max velocity
    VtAzmdir=306.76219

    Hgrid=sm.hurricanewavecontourcem(xgrid,ygrid,longCenter,latCenter,Hmax,RVmax,VtAzmdir,'cem','gc','contour')


    #EXAMPLE 2

    #Creating calculation mesh
    xgrid,ygrid=np.meshgrid(np.linspace(-98,-68,100),np.linspace(16,44,100))

    #Longitude of Hurricane Katrine best track
    longtrack=[-75.1,-75.7,-76.2,-76.5,-76.9,-77.7,-78.4,-79.0,-79.6,-80.1,-80.3,-81.3,\
        -82.0,-82.6,-83.3,-84.0,-84.7,-85.3,-85.9,-86.7,-87.7,-88.6,-89.2,-89.6,\
        -89.6,-89.6,-89.6,-89.6,-89.1,-88.6,-88.0,-87.0,-85.3,-82.9]

    #Latitude of Hurricane Katrine best track
    lattrack=[23.1,23.4,23.8,24.5,25.4,26.0,26.1,26.2,26.2,26.0,25.9,25.4,\
        25.1,24.9,24.6,24.4,24.4,24.5,24.8,25.2,25.7,26.3,27.2,28.2,\
        29.3,29.5,30.2,31.1,32.6,34.1,35.6,37.0,38.6,40.1]

    #Hurricane Katrina maximum significant wave height
    Hmax=[0,0,0,4.3788,4.9295,5.5527,6.2110,6.8516,7.5428,9.1513,8.5021,8.6332,10.1511,11.3434,\
        12.3171,13.5606,14.1226,14.4931,14.1972,19.9683,24.0121,24.9821,23.0419,19.9342,16.5366,\
        14.5246,14.8050,0,0,0,0,0,0,0]

    #Hurricane Katrina radius from hurricane center to a location of maximum hurricane wind velocity (m)
    RVmax=[0,0,0,8.0290e+004,5.6029e+004,4.2063e+004,3.6769e+004,3.3849e+004,3.1352e+004,3.3405e+004,3.3773e+004,\
        3.2657e+004,3.1122e+004,2.7037e+004,2.6512e+004,3.3476e+004,3.0881e+004,4.0266e+004,3.2433e+004,\
        5.1747e+004,5.7297e+004,6.2750e+004,5.3376e+004,4.3074e+004,3.1790e+004,4.3114e+004,2.7800e+004,\
        0,0,0,0,0,0,0]

    #Hurricane Katrina velocity azimuth (bearing) in (Degree)
    VtAzmdir=[0.00000,298.67291,311.22135,338.70264,338.13626,309.94476,279.18860,280.65053,270.13245,\
    246.10095,240.96690,241.20181,244.79591,249.93382,244.88325,252.71384,270.14459,280.49918,\
    298.94148,299.05364,299.18896,306.76219,329.36839,340.59069,0.00000,0.00000,0.00,\
        0.00000,15.67775,15.42254,18.00215,29.63266,39.49673,50.29744]

    Hgrid=sm.hurricanewavecontourcem(xgrid,ygrid,longtrack[3:27],lattrack[3:27],Hmax[3:27],RVmax[3:27],VtAzmdir[3:27],'cem','gc','contour')


References
----------

Data

* www.nhc.noaa.gov/data/
* www.nhc.noaa.gov/data/hurdat/hurdat2-format-nencpac.pdf
* coast.noaa.gov/hurricanes
* www.aoml.noaa.gov/hrd/data_sub/re_anal.html

Department of the Army, Waterways Experiment Station, Corps of Engineers, 
and Coastal Engineering Research Center (1984), 
Shore Protection Manual, Washington, 
D.C., vol. 1, 4th ed., 532 pp.

U.S. Army Corps of Engineers (2015). 
Coastal Engineering Manual. 
Engineer Manual 1110-2-1100, Washington, D.C.: U.S. Army Corps of Engineers.

Young, I. R. (1988). 
Parametric hurricane wave prediction model. 
Journal of Waterway, Port, Coastal, and Ocean Engineering, 114(5), 637-652.

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
