.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-08-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

scientimate.intersectgc
=======================

.. code:: python

    latint1, lonint1, latint2, lonint2, latintn, lonintn, latints, lonints, latinte, loninte, latintw, lonintw = scientimate.intersectgc(lat1, lon1, lat2, lon2, lat3, lon3, lat4, lon4, dispout='no')

Description
-----------

Find intersection point between two line segments (line edges) on Great Circle

Inputs
------

lat1
                            :Latitude (y) of start point of first segment as (lon1,lat1) in (Degree)
lon1
                            :Longitude (x) of start point of first segment as (lon1,lat1) in (Degree)
lat2
                            :Latitude (y) of end point of first segment as (lon2,lat2) in (Degree)
lon2
                            :Longitude (x) of end point of first segment as (lon2,lat2) in (Degree)
                            :First segment: p1(lon1,lat1) to p2(lon2,lat2)
lat3
                            :Latitude (y) of start point of second segment as (lon3,lat3) in (Degree)
lon3
                            :Longitude (x) of start point of second segment as (lon3,lat3) in (Degree)
lat4
                            :Latitude (y) of end point of second segment as (lon4,lat4) in (Degree)
lon4
                            :Longitude (x) of end point of second segment as (lon4,lat4) in (Degree)
                            :Second segment: p3(lon3,lat3) to p4(lon4,lat4)
dispout='no'
                            :Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

latint1
                            :Latitude of first intersection point between two segments in (Degree)
lonint1
                            :Longitude of first intersection point between two segments in (Degree)
latint2
                            :Latitude of second intersection point between two segments in (Degree)
lonint2
                            :Longitude of second intersection point between two segments in (Degree)
                            :Two path on sphere have two intesections on opposite side of sphere 
                            :(latint1,lonint1) and (latint2,lonint2) are antipodes of each other
latintn
                            :Latitude of intersection point between two segments in Northern Hemisphere in (Degree)
lonintn
                            :Longitude of intersection point between two segments in Northern Hemisphere in (Degree)
latints
                            :Latitude of intersection point between two segments in Southern Hemisphere in (Degree)
lonints
                            :Longitude of intersection point between two segments in Southern Hemisphere in (Degree)
latinte
                            :Latitude of intersection point between two segments in Eastern Hemisphere in (Degree)
loninte
                            :Longitude of intersection point between two segments in Eastern Hemisphere in (Degree)
latintw
                            :Latitude of intersection point between two segments in Western Hemisphere in (Degree)
lonintw
                            :Longitude of intersection point between two segments in Western Hemisphere in (Degree)

Examples
--------

.. code:: python

    import scientimate as sm

    #Segment 1:
    lat1=-90
    lon1=0
    lat2=90 
    lon2=0
    #Segment 2:
    lat3=0
    lon3=-90
    lat4=0
    lon4=90
    latint1,lonint1,latint2,lonint2,latintn,lonintn,latints,lonints,latinte,loninte,latintw,lonintw=sm.intersectgc(lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4,'yes')

    #Segment 1:
    lat1=29.079710
    lon1=-90.457758
    lat2=29.344989
    lon2=-90.476782
    #Segment 2:
    lat3=29.190111
    lon3=-90.598023
    lat4=29.206355
    lon4=-90.337782
    latint1,lonint1,latint2,lonint2,latintn,lonintn,latints,lonints,latinte,loninte,latintw,lonintw=sm.intersectgc(lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4,'yes')
    #latint1 = -29.198
    #lonint1 =  89.534
    #latint2 =  29.198
    #lonint2 = -90.466

    #Segment 1:
    lat1=47.94713
    lon1=-131.211073
    lat2=24.207076
    lon2=-83.815088
    #Segment 2:
    lat3=28.257645
    lon3=-95.964404
    lat4=28.343359
    lon4=-42.233815
    latint1,lonint1,latint2,lonint2,latintn,lonintn,latints,lonints,latinte,loninte,latintw,lonintw=sm.intersectgc(lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4,'yes')
    #latint1 =  29.396
    #lonint1 = -89.930
    #latint2 = -29.396
    #lonint2 =  90.070

    #Segment 1:
    lat1=[29.079710,47.94713]
    lon1=[-90.457758,-131.211073]
    lat2=[29.344989,24.207076]
    lon2=[-90.476782,-83.815088]
    #Segment 2:
    lat3=[29.190111,28.257645]
    lon3=[-90.598023,-95.964404]
    lat4=[29.206355,28.343359]
    lon4=[-90.337782,-42.233815]
    latint1,lonint1,latint2,lonint2,latintn,lonintn,latints,lonints,latinte,loninte,latintw,lonintw=sm.intersectgc(lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4,'yes')

References
----------

| https://en.wikipedia.org/wiki/Geographic_coordinate_conversion
| https://en.wikipedia.org/wiki/Geographic_coordinate_system
| https://en.wikipedia.org/wiki/Earth_radius
| https://stackoverflow.com/questions/29465468/python-intersection-point-of-two-great-circles-lat-long
| https://stackoverflow.com/questions/2954337/great-circle-rhumb-line-intersection
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
