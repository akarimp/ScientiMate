function [arclen, azimuthdir, metedir] = distancegc(lat1, lon1, lat2, lon2, CalcMethod, R, dispout)
%{
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
%}

%--------------------------------------------------------------------------
%CODE
%--------------------------------------------------------------------------
%Assign default values

switch nargin
    case 4
        CalcMethod='haversine'; R=6371000; dispout='no';
    case 5
        R=6371000; dispout='no';
    case 6
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(lat1)==1
    lat1=lat1';
end

if isrow(lon1)==1
    lon1=lon1';
end

if isrow(lat2)==1
    lat2=lat2';
end

if isrow(lon2)==1
    lon2=lon2';
end

%--------------------------------------------------------------------------
%Calculating distace between start and end of the line (total length)

%Converting to radian
lat1rad=deg2rad(lat1);
lon1rad=deg2rad(lon1);
lat2rad=deg2rad(lat2);
lon2rad=deg2rad(lon2);

deltalatrad21=lat2rad-lat1rad;
deltalonrad21=lon2rad-lon1rad;

%Calculating distance using spherical law of cosines
if strcmp(CalcMethod,'cos')==1

    deltasigma=acos(sin(lat1rad).*sin(lat2rad)+cos(lat1rad).*cos(lat2rad).*cos(deltalonrad21)); %Central angle
    arclen=R.*deltasigma; %Total distance of the line 

%Calculating distance using Haversine formula
elseif strcmp(CalcMethod,'haversine')==1

    a=(sin(deltalatrad21./2)).^2+cos(lat1rad).*cos(lat2rad).*(sin(deltalonrad21./2)).^2; %Angular distance in radians
    %c=2.*atan2(sqrt(a),sqrt(1-a)); %Square of half the chord length between the points (central angle)
    deltasigma=2.*atan2(sqrt(a),sqrt(1-a)); %Central angle
    arclen=R.*deltasigma; %Total distance of the line 

%Calculating distance using Vincenty formula
elseif strcmp(CalcMethod,'vincenty')==1

    deltasigma=atan2(sqrt((cos(lat2rad).*sin(deltalonrad21)).^2+(cos(lat1rad).*sin(lat2rad)-sin(lat1rad).*cos(lat2rad).*cos(deltalonrad21)).^2),sin(lat1rad).*sin(lat2rad)+cos(lat1rad).*cos(lat2rad).*cos(deltalonrad21)); %Central angle
    arclen=R.*deltasigma; %Total distance of the line 

end

%--------------------------------------------------------------------------
%Calculating azimuth (bearing) from start point to end point

azimuthdirrad=atan2(sin(deltalonrad21).*cos(lat2rad),cos(lat1rad).*sin(lat2rad)-sin(lat1rad).*cos(lat2rad).*cos(deltalonrad21)); %Azimuth (bearing) in radian
azimuthdir=rad2deg(azimuthdirrad); %Azimuth (bearing) in degree

%Add 360 to all numbers to have them all positive
%Use mod(360) to take care of the ones larger than 360
azimuthdir=mod((azimuthdir+360),360); 

%--------------------------------------------------------------------------
%Calculating meteorological direction
metedir=azimuthdir+180; %Convert azimuth (bearing or compass) angle to meteorological direction

%Add 360 to all numbers to have them all positive
%Use mod(360) to take care of the ones larger than 360
metedir=mod((metedir+360),360); 

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1

    scatter(lon1,lat1,'filled')
    hold on
    scatter(lon2,lat2,'filled')

    xlabel('Longitude (Degree)')
    ylabel('Latitude (Degree)')
    legend('Start Point','End Point')

end

%--------------------------------------------------------------------------
