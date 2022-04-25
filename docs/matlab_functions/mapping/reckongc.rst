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

reckongc
========

.. code:: MATLAB

    [lat2, lon2] = reckongc(lat1, lon1, arclen, azimuthdir, R, dispout)

Description
-----------

Calculate end point (Latitude,Longitude) from start point (Latitude,Longitude) and distance and azimuth (bearing) using Great Circle

Inputs
------

lat1
    Latitude (y) of start point (first point) in (Degree)
lon1
    Longitude (x) of start point (first point) in (Degree)
arclen
    Total distance from start point to end point in (m)
azimuthdir
    | Azimuth (bearing or compass direction) from start point to end point in (Degree)
    | 0 (degree): toward North, 90 (degree): toward East, 180 (degree): toward South, 270 (degree): toward West 
R=6371000;
    Earth radius in (m), mean earth radius=6371000 m
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

lat2
    Latitude (y) of end point (last point) in (Degree)
lon2
    Longitude (x) of end point (last point) in (Degree)

Examples
--------

.. code:: MATLAB

    lat1=29.5; %First point 
    lon1=-89.4; %First point 
    arclen=22239; %Arc length
    azimuthdir=0; %Azimuth
    [lat2,lon2]=reckongc(lat1,lon1,arclen,azimuthdir);

    lat1=[29.5;29]; %First point 
    lon1=[-89.4;-89]; %First point 
    arclen=[22239;147410]; %Arc length
    azimuthdir=[0;319.21]; %Azimuth
    [lat2,lon2]=reckongc(lat1,lon1,arclen,azimuthdir,6371000,'yes');

References
----------

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
