function [latint1, lonint1, latint2, lonint2, latintn, lonintn, latints, lonints, latinte, loninte, latintw, lonintw] = intersectgc(lat1, lon1, lat2, lon2, lat3, lon3, lat4, lon4, dispout)
%{
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

intersectgc
===========

.. code:: MATLAB

    [latint1, lonint1, latint2, lonint2, latintn, lonintn, latints, lonints, latinte, loninte, latintw, lonintw] = intersectgc(lat1, lon1, lat2, lon2, lat3, lon3, lat4, lon4, dispout)

Description
-----------

Find intersection point between two line segments (line edges) on Great Circle

Inputs
------

lat1
    Latitude (y) of start point of first segment as (lon1,lat1) in (Degree)
lon1
    Longitude (x) of start point of first segment as (lon1,lat1) in (Degree)
lat2
    Latitude (y) of end point of first segment as (lon2,lat2) in (Degree)
lon2
    | Longitude (x) of end point of first segment as (lon2,lat2) in (Degree)
    | First segment: p1(lon1,lat1) to p2(lon2,lat2)
lat3
    Latitude (y) of start point of second segment as (lon3,lat3) in (Degree)
lon3
    Longitude (x) of start point of second segment as (lon3,lat3) in (Degree)
lat4
    Latitude (y) of end point of second segment as (lon4,lat4) in (Degree)
lon4
    | Longitude (x) of end point of second segment as (lon4,lat4) in (Degree)
    | Second segment: p3(lon3,lat3) to p4(lon4,lat4)
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

latint1
    Latitude of first intersection point between two segments in (Degree)
lonint1
    Longitude of first intersection point between two segments in (Degree)
latint2
    Latitude of second intersection point between two segments in (Degree)
lonint2
    | Longitude of second intersection point between two segments in (Degree)
    | Two path on sphere have two intesections on opposite side of sphere 
    | (latint1,lonint1) and (latint2,lonint2) are antipodes of each other
latintn
    Latitude of intersection point between two segments in Northern Hemisphere in (Degree)
lonintn
    Longitude of intersection point between two segments in Northern Hemisphere in (Degree)
latints
    Latitude of intersection point between two segments in Southern Hemisphere in (Degree)
lonints
    Longitude of intersection point between two segments in Southern Hemisphere in (Degree)
latinte
    Latitude of intersection point between two segments in Eastern Hemisphere in (Degree)
loninte
    Longitude of intersection point between two segments in Eastern Hemisphere in (Degree)
latintw
    Latitude of intersection point between two segments in Western Hemisphere in (Degree)
lonintw
    Longitude of intersection point between two segments in Western Hemisphere in (Degree)

Examples
--------

.. code:: MATLAB

    %Segment 1:
    lat1=-90;
    lon1=0;
    lat2=90; 
    lon2=0;
    %Segment 2:
    lat3=0;
    lon3=-90;
    lat4=0;
    lon4=90;
    [latint1,lonint1,latint2,lonint2,latintn,lonintn,latints,lonints,latinte,loninte,latintw,lonintw]=intersectgc(lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4,'yes');

    %Segment 1:
    lat1=29.079710;
    lon1=-90.457758;
    lat2=29.344989;
    lon2=-90.476782;
    %Segment 2:
    lat3=29.190111;
    lon3=-90.598023;
    lat4=29.206355;
    lon4=-90.337782;
    [latint1,lonint1,latint2,lonint2,latintn,lonintn,latints,lonints,latinte,loninte,latintw,lonintw]=intersectgc(lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4,'yes');
    %latint1 = -29.198
    %lonint1 =  89.534
    %latint2 =  29.198
    %lonint2 = -90.466

    %Segment 1:
    lat1=47.94713;
    lon1=-131.211073;
    lat2=24.207076;
    lon2=-83.815088;
    %Segment 2:
    lat3=28.257645;
    lon3=-95.964404;
    lat4=28.343359;
    lon4=-42.233815;
    [latint1,lonint1,latint2,lonint2,latintn,lonintn,latints,lonints,latinte,loninte,latintw,lonintw]=intersectgc(lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4,'yes');
    %latint1 =  29.396
    %lonint1 = -89.930
    %latint2 = -29.396
    %lonint2 =  90.070

    %Segment 1:
    lat1=[29.079710;47.94713];
    lon1=[-90.457758;-131.211073];
    lat2=[29.344989;24.207076];
    lon2=[-90.476782;-83.815088];
    %Segment 2:
    lat3=[29.190111;28.257645];
    lon3=[-90.598023;-95.964404];
    lat4=[29.206355;28.343359];
    lon4=[-90.337782;-42.233815];
    [latint1,lonint1,latint2,lonint2,latintn,lonintn,latints,lonints,latinte,loninte,latintw,lonintw]=intersectgc(lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4,'yes');

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
%}

%--------------------------------------------------------------------------
%CODE
%--------------------------------------------------------------------------
%Assign default values

switch nargin
    case 8
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

if isrow(lat3)==1
    lat3=lat3';
end

if isrow(lon3)==1
    lon3=lon3';
end

if isrow(lat4)==1
    lat4=lat4';
end

if isrow(lon4)==1
    lon4=lon4';
end

%--------------------------------------------------------------------------
%Find location of intersection

%Converting to radian
lat1rad=deg2rad(lat1);
lon1rad=deg2rad(lon1);
lat2rad=deg2rad(lat2);
lon2rad=deg2rad(lon2);
lat3rad=deg2rad(lat3);
lon3rad=deg2rad(lon3);
lat4rad=deg2rad(lat4);
lon4rad=deg2rad(lon4);

%Earth radius
a=63781370; %Equatorial radius (m)
b=63567523; %Polar radius (m)

%Distance from the surface to the z axis along the ellipsoid normal (transverse radius of curvature)
N1(:,1)=a.^2/sqrt(a.^2.*(cos(lat1rad))+b.^2.*(sin(lat1rad)).^2);
N2(:,1)=a.^2/sqrt(a.^2.*(cos(lat2rad))+b.^2.*(sin(lat2rad)).^2);
N3(:,1)=a.^2/sqrt(a.^2.*(cos(lat3rad))+b.^2.*(sin(lat3rad)).^2);
N4(:,1)=a.^2/sqrt(a.^2.*(cos(lat4rad))+b.^2.*(sin(lat4rad)).^2);

%Height from earth surface
h=0; 

%Converting geodetic (lat,lon,h) to ECEF (x,y,z) coordinates
x1(:,1)=(N1+h).*cos(lat1rad).*cos(lon1rad);
y1(:,1)=(N1+h).*cos(lat1rad).*sin(lon1rad);
z1(:,1)=(N1+h).*sin(lat1rad);
x2(:,1)=(N2+h).*cos(lat2rad).*cos(lon2rad);
y2(:,1)=(N2+h).*cos(lat2rad).*sin(lon2rad);
z2(:,1)=(N2+h).*sin(lat2rad);
x3(:,1)=(N3+h).*cos(lat3rad).*cos(lon3rad);
y3(:,1)=(N3+h).*cos(lat3rad).*sin(lon3rad);
z3(:,1)=(N3+h).*sin(lat3rad);
x4(:,1)=(N4+h).*cos(lat4rad).*cos(lon4rad);
y4(:,1)=(N4+h).*cos(lat4rad).*sin(lon4rad);
z4(:,1)=(N4+h).*sin(lat4rad);

%Normal planes that contain great circles
NP1=cross([x1,y1,z1],[x2,y2,z2]);
NP2=cross([x3,y3,z3],[x4,y4,z4]);

%Line of intersection between two planes
L=cross(NP1,NP2);

%Intersection points
X1=zeros(length(lat1(:,1)),3); %Pre-assigning vector
for i=1:length(lat1(:,1))
    X1(i,:)=L(i,:)./sqrt(L(i,1)^2+L(i,2)^2+L(i,3)^2); %Normalizing
end
X2=-X1;

%Two path on sphere have two intesections on opposite side of sphere which are antipodes of each other
latint1=rad2deg(asin(X1(:,3)));
lonint1=rad2deg(atan2(X1(:,2),X1(:,1)));
latint2=rad2deg(asin(X2(:,3)));
lonint2=rad2deg(atan2(X2(:,2),X2(:,1)));

%--------------------------------------------------------------------------
%Intersection points based on their hemisphere locations

%Checking hemisphere location of intersections
SegInsctChkN1=(latint1>=0 & latint1<=90); %Checking if latint1 is in Northern Hemisphere
SegInsctChkS1=(latint1>=-90 & latint1<0); %Checking if latint1 is in Southern Hemisphere
SegInsctChkE1=(lonint1>=0 & lonint1<=180); %Checking if lonint1 is in Eastern Hemisphere
SegInsctChkW1=(lonint1>=-180 & lonint1<0); %Checking if lonint1 is in Western Hemisphere

SegInsctChkN2=(latint2>=0 & latint2<=90); %Checking if latint2 is in Northern Hemisphere
SegInsctChkS2=(latint2>=-90 & latint2<0); %Checking if latint2 is in Southern Hemisphere
SegInsctChkE2=(lonint2>=0 & lonint2<=180); %Checking if lonint2 is in Eastern Hemisphere
SegInsctChkW2=(lonint2>=-180 & lonint2<0); %Checking if lonint2 is in Western Hemisphere

%Pre-assigning vector
latintn=zeros(length(lat1(:,1)),1); %Pre-assigning vector
lonintn=zeros(length(lat1(:,1)),1); %Pre-assigning vector
latints=zeros(length(lat1(:,1)),1); %Pre-assigning vector
lonints=zeros(length(lat1(:,1)),1); %Pre-assigning vector
latinte=zeros(length(lat1(:,1)),1); %Pre-assigning vector
loninte=zeros(length(lat1(:,1)),1); %Pre-assigning vector
latintw=zeros(length(lat1(:,1)),1); %Pre-assigning vector
lonintw=zeros(length(lat1(:,1)),1); %Pre-assigning vector

%Intersection points based on their locations
%Intersection points in Northern Hemisphere
latintn(SegInsctChkN1==1)=latint1(SegInsctChkN1==1);
latintn(SegInsctChkN2==1)=latint2(SegInsctChkN2==1);
lonintn(SegInsctChkN1==1)=lonint1(SegInsctChkN1==1);
lonintn(SegInsctChkN2==1)=lonint2(SegInsctChkN2==1);

%Intersection points in Southern Hemisphere
latints(SegInsctChkS1==1)=latint1(SegInsctChkS1==1);
latints(SegInsctChkS2==1)=latint2(SegInsctChkS2==1);
lonints(SegInsctChkS1==1)=lonint1(SegInsctChkS1==1);
lonints(SegInsctChkS2==1)=lonint2(SegInsctChkS2==1);

%Intersection points in Eastern Hemisphere
latinte(SegInsctChkE1==1)=latint1(SegInsctChkE1==1);
latinte(SegInsctChkE2==1)=latint2(SegInsctChkE2==1);
loninte(SegInsctChkE1==1)=lonint1(SegInsctChkE1==1);
loninte(SegInsctChkE2==1)=lonint2(SegInsctChkE2==1);

%Intersection points in Western Hemisphere
latintw(SegInsctChkW1==1)=latint1(SegInsctChkW1==1);
latintw(SegInsctChkW2==1)=latint2(SegInsctChkW2==1);
lonintw(SegInsctChkW1==1)=lonint1(SegInsctChkW1==1);
lonintw(SegInsctChkW2==1)=lonint2(SegInsctChkW2==1);

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1

    for i=1:length(lat1(:,1))
        plot([lon1(i,1),lon2(i,1)],[lat1(i,1),lat2(i,1)])
        hold on
        plot([lon3(i,1),lon4(i,1)],[lat3(i,1),lat4(i,1)])
    end
    
    scatter(lonint1,latint1,'filled')
    scatter(lonint2,latint2,'filled')

    plot([-180,180],[0,0],'k')
    plot([0,0],[-90,90],'k')

    xlabel('Longitude (Degree)')
    ylabel('Latitude (Degree)')
    xlim([-180,180])
    ylim([-90,90])
    %legend('First Segment','Second Segment','Intersection Point')

end

%--------------------------------------------------------------------------
