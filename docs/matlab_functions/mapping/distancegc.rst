.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-07-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

distancegc
==========

.. code:: MATLAB

    [arclen, azimuthdir, metedir] = distancegc(lat1, lon1, lat2, lon2, CalcMethod, R, dispout)

Description
-----------

Calculate distance and azimuth (bearing) between (Latitude,Longitude) points using Great Circle

Inputs
------

lat1
    Latitude (y) of start point (first point) in (Degree)
lon1
    Longitude (x) of start point (first point) in (Degree)
lat2
    Latitude (y) of end point (last point) in (Degree)
lon2
    Longitude (x) of end point (last point) in (Degree)
CalcMethod='haversine';
    | Distance calculation method 
    | 'cos': Spherical law of cosines
    | 'haversine': Haversine formula
    | 'vincenty': Vincenty formula, Vincenty (1975)
R=6371000;
    Earth radius in (m), mean earth radius=6371000 m
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

arclen
    Total distance from start point to end point in (m)
azimuthdir
    | Azimuth (bearing or compass direction) from start point to end point in (Degree)
    | 0 (degree): toward North, 90 (degree): toward East, 180 (degree): toward South, 270 (degree): toward West 
metedir
    | Meteorological direction from start point to end point in (Degree)
    | 0 (degree): from North, 90 (degree): from East, 180 (degree): from South, 270 (degree): from West 

Examples
--------

.. code:: MATLAB

    lat1=29.5; %First point 
    lon1=-89.4; %First point 
    lat2=29.7; %last point
    lon2=-89.4; %last point
    [arclen,azimuthdir,metedir]=distancegc(lat1,lon1,lat2,lon2);

    lat1=[29.5;29]; %First point 
    lon1=[-89.4;-89]; %First point 
    lat2=[29.7;30]; %Last point
    lon2=[-89.4;-90]; %Last point
    [arclen,azimuthdir,metedir]=distancegc(lat1,lon1,lat2,lon2,'haversine',6371000,'yes');

References
----------

Vincenty, T. (1975). 
Direct and inverse solutions of geodesics on the ellipsoid with application of nested equations. 
Survey review, 23(176), 88-93.

| http://www.movable-type.co.uk/scripts/latlong.html
| https://en.wikipedia.org/wiki/Great-circle_distance
| https://en.wikipedia.org/wiki/Great-circle_navigation
| http://www.codeguru.com/cpp/cpp/algorithms/article.php/c5115/Geographic-Distance-and-Azimuth-Calculations.htm

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
