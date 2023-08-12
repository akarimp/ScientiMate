function [wind_fetch, z_mean, z_point, bed_slope_mean, x_fetch, y_fetch, z_fetch, dist_along_fetch] = windfetch(xgrid, ygrid, zgrid, x_point, y_point, winddir, waterlevel, shorelinelevel, n_midpoints, distCalcMethod, interpMethod, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2020-12-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

windfetch
=========

.. code:: MATLAB

    [wind_fetch, z_mean, z_point, bed_slope_mean, x_fetch, y_fetch, z_fetch, dist_along_fetch] = windfetch(xgrid, ygrid, zgrid, x_point, y_point, winddir, waterlevel, shorelinelevel, n_midpoints, distCalcMethod, interpMethod, dispout)

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
winddir=0;
    | Meteorological wind direction in (Degree)
    | It represents a direction wind comes from and is measured counter-clockwise from the North
    | 0 (degree): from North, 90 (degree): from East, 180 (degree): from South, 270 (degree): from West
waterlevel=0;
    | Water surface level (Should have the same datum as zgrid)
    | Size of waterlevel should be either 1 or should be equal to size of winddir 
    | If water level is zero (waterlevel=0), it means that the water level is at a level of surface elevation (z=0) in bathymetry data  
    | If water surface level is positive (waterlevel>0), it means that the water level is above zero elevation (z=0) in bathymetry data  
    |     ,therefore water level will be subtracted from bathymetry elevation 
    | If water surface level is negative (waterlevel<0), it means that the water level is below zero elevation (z=0) in bathymetry data  
    |     ,therefore water level will be added to bathymetry elevation 
shorelinelevel=0;
    | Shoreline elevation
    | z>=shorelinelevel is considered land, z<shorelinelevel is considered water
n_midpoints=100;
    | Number of middle points generated between first point and last point
    | if n_midpoints=1 then 1 point between first point and last point is generated
distCalcMethod='gc';
    | Distance calculation method 
    | 'cart': Distances are calculated on cartesian coordinate
    | 'gc': Distances are calculated on Great Circle based on Vincenty formula, Vincenty (1975)
    | Earth radius considered as mean earth radius=6371000 m
interpMethod='nearest';
    | Interpolation method 
    | 'linear': Use default or 'linear' method to interpolate
    | 'nearest': Use nearest neighbor method to interpolate
dispout='no';
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

.. code:: MATLAB

    [xgrid,ygrid]=meshgrid(linspace(-5,5,100),linspace(-5,5,100));
    R=sqrt(xgrid.^2+ygrid.^2);
    zgrid=sin(R)./R;
    x_point=-4;
    y_point=-2;
    winddir=90;
    waterlevel=0;
    shorelinelevel=0;
    n_midpoints=100;
    [wind_fetch,z_mean,z_point,bed_slope_mean,x_fetch,y_fetch,z_fetch,dist_along_fetch]=windfetch(xgrid,ygrid,zgrid,x_point,y_point,winddir,waterlevel,shorelinelevel,n_midpoints,'cart','nearest','yes');

    [xgrid,ygrid]=meshgrid(linspace(-5,5,100),linspace(-5,5,100));
    zgrid=ygrid.*sin(xgrid)-xgrid.*cos(ygrid)-3;
    x_point=0;
    y_point=0;
    winddir=[0:30:360];
    waterlevel=ones(1,numel(winddir))*0.1;
    shorelinelevel=0;
    n_midpoints=100;
    [wind_fetch,z_mean,z_point,bed_slope_mean,x_fetch,y_fetch,z_fetch,dist_along_fetch]=windfetch(xgrid,ygrid,zgrid,x_point,y_point,winddir,waterlevel,shorelinelevel,n_midpoints,'cart','nearest','yes');

    %Download Persina Gulf and Gulf of Oman with coordinate exteneds xmin=56, xmax=63, ymin=19, ymax=28 from https://maps.ngdc.noaa.gov/viewers/grid-extract/index.html
    xyzfilename='PersianGulf.csv'; %e.g. xyzfilename='PersianGulf.csv'
    xyzfilelocation=pwd; %e.g. xyzfilelocation='C:\datafolder'
    [x,y,z]=readxyz(xyzfilename,xyzfilelocation,1,'all');
    [xgrid,ygrid,zgrid]=interpxyz2grid(x,y,z,100,'points',nanmin(x),nanmax(x),nanmin(y),nanmax(y),nanmin(z),nanmax(z),'all','nearest','no');
    x_point=58.0; %Or x_point=52.0
    y_point=24.5; %Or y_point=26.0
    winddir=[0;40;135;280];
    waterlevel=0;
    shorelinelevel=0;
    n_midpoints=100;
    [wind_fetch,z_mean,z_point,bed_slope_mean,x_fetch,y_fetch,z_fetch,dist_along_fetch]=windfetch(xgrid,ygrid,zgrid,x_point,y_point,winddir,waterlevel,shorelinelevel,n_midpoints,'gc','nearest','yes');

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
%}

%--------------------------------------------------------------------------
%CODE
%--------------------------------------------------------------------------
%Assign default values

switch nargin
    case 5
        winddir=0; waterlevel=0; shorelinelevel=0; n_midpoints=100; distCalcMethod='gc'; interpMethod='nearest'; dispout='no';
    case 6
        waterlevel=0; shorelinelevel=0; n_midpoints=100; distCalcMethod='gc'; interpMethod='nearest'; dispout='no';
    case 7
        shorelinelevel=0; n_midpoints=100; distCalcMethod='gc'; interpMethod='nearest'; dispout='no';
    case 8
        n_midpoints=100; distCalcMethod='gc'; interpMethod='nearest'; dispout='no';
    case 9
        distCalcMethod='gc'; interpMethod='nearest'; dispout='no';
    case 10
        interpMethod='nearest'; dispout='no';
    case 11
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(x_point)==1
    x_point=x_point';
end

if isrow(y_point)==1
    y_point=y_point';
end

if isrow(winddir)==1
    winddir=winddir';
end

if isrow(waterlevel)==1
    waterlevel=waterlevel';
end

if isrow(shorelinelevel)==1
    shorelinelevel=shorelinelevel';
end

%--------------------------------------------------------------------------
%Generate points on boundary of domain

%Domain extend
xmin=min(xgrid(:)); %Minimum of x value
xmax=max(xgrid(:)); %Maximum of x value
ymin=min(ygrid(:)); %Minimum of y value
ymax=max(ygrid(:)); %Maximum of y value

%Points on north, east, south, and west boundary rotates clockwise
y_boundary_north = ones(1,15)*ymax;
x_boundary_north = linspace(xmin,xmax,15);

y_boundary_east = linspace(ymax,ymin,15);
x_boundary_east = ones(1,15)*xmax;

y_boundary_south = ones(1,15)*ymin;
x_boundary_south = linspace(xmax,xmin,15);

y_boundary_west = linspace(ymin,ymax,15);
x_boundary_west = ones(1,15)*xmin;

%Points on domain boundary, starts from east and rotate counter-clockwise
%y_boundary = [y_boundary_east(y_boundary_east>=y_point),y_boundary_north(1,2:end),y_boundary_west(1,2:end),y_boundary_south(1,2:end-1),y_boundary_east(y_boundary_east<y_point)];
%x_boundary = [x_boundary_east(y_boundary_east>=y_point),x_boundary_north(1,2:end),x_boundary_west(1,2:end),x_boundary_south(1,2:end-1),x_boundary_east(y_boundary_east<y_point)];

%Points on domain boundary, starts from north and rotate clockwise
y_boundary = [y_boundary_north(x_boundary_north>=x_point),y_boundary_east(1,2:end),y_boundary_south(1,2:end),y_boundary_west(1,2:end-1),y_boundary_north(x_boundary_north<x_point)];
x_boundary = [x_boundary_north(x_boundary_north>=x_point),x_boundary_east(1,2:end),x_boundary_south(1,2:end),x_boundary_west(1,2:end-1),x_boundary_north(x_boundary_north<x_point)];

%Add the the last point before the first point and the first point after the last point to have a closed loop
y_boundary = [y_boundary(1,end),y_boundary,y_boundary(1,1)];
x_boundary = [x_boundary(1,end),x_boundary,x_boundary(1,1)];

%Change to column vector
y_boundary = y_boundary.';
x_boundary = x_boundary.';

%--------------------------------------------------------------------------
%Generating points on transect line starts from point of interest toward points on boundary 

%Calculate using cartesian formula
if strcmp(distCalcMethod,'cart')==1
    
    %Generating points from (x1,y1) to (x2,y2)
    n_midpoints(n_midpoints<1)=1;
    x_transect_360=linspace(x_point,x_boundary,n_midpoints+2);
    y_transect_360=linspace(y_point,y_boundary,n_midpoints+2);

%Calculate using great circle formula
elseif strcmp(distCalcMethod,'gc')==1

    %Converte to radian
    lat1_rad=deg2rad(y_point);
    lon1_rad=deg2rad(x_point);
    lat2_rad=deg2rad(y_boundary);
    lon2_rad=deg2rad(x_boundary);

    delta_lat_21_rad=lat2_rad-lat1_rad;
    delta_lon_21_rad=lon2_rad-lon1_rad;

    %Calculating azimuth (bearing) between start and end of the line
    azimuthdir_rad=atan2(sin(delta_lon_21_rad).*cos(lat2_rad),cos(lat1_rad).*sin(lat2_rad)-sin(lat1_rad).*cos(lat2_rad).*cos(delta_lon_21_rad)); %Azimuth (bearing) in radian

    R=6371000; %Earth radius in (m), mean earth radius=6371000 m
    deltasigma=atan2(sqrt((cos(lat2_rad).*sin(delta_lon_21_rad)).^2+(cos(lat1_rad).*sin(lat2_rad)-sin(lat1_rad).*cos(lat2_rad).*cos(delta_lon_21_rad)).^2),sin(lat1_rad).*sin(lat2_rad)+cos(lat1_rad).*cos(lat2_rad).*cos(delta_lon_21_rad)); %Central angle

    %Generating lat and lon between (lon1,lat1) and (lon2,lat2), first method
    n_midpoints(n_midpoints<1)=1;
    deltasigmai=linspace(0,deltasigma,n_midpoints+2); %Central angle for each point from first point, deltasigma=d/R
    lati_rad=asin(sin(lat1_rad).*cos(deltasigmai)+cos(lat1_rad).*sin(deltasigmai).*cos(azimuthdir_rad));
    loni_rad=lon1_rad+atan2(sin(azimuthdir_rad).*sin(deltasigmai).*cos(lat1_rad),cos(deltasigmai)-sin(lat1_rad).*sin(lati_rad));

    %Converte to degree
    y_transect_360=rad2deg(lati_rad);
    x_transect_360=rad2deg(loni_rad);

end

%--------------------------------------------------------------------------
%Calculate z (elevation) profile along transect path from boundary domain toward point of interest

%Interpolating the z (elevation) values for given (x,y) from input data
if strcmp(interpMethod,'linear')==1
    z_transect_360=interp2(xgrid,ygrid,zgrid,x_transect_360,y_transect_360);

    %Replacing NaN data point resulted from default method with ones from nearest method
    if sum(isnan(z_transect_360(:)))>0
        z_nearest=interp2(xgrid,ygrid,zgrid,x_transect_360,y_transect_360,'nearest');
        z_transect_360(isnan(z_transect_360)==1)=z_nearest(isnan(z_transect_360)==1); 
    end

%Interpolating data into grid using nearest neighbor method
elseif strcmp(interpMethod,'nearest')==1

    %First method
    z_transect_360=interp2(xgrid,ygrid,zgrid,x_transect_360,y_transect_360,'nearest');

end

%--------------------------------------------------------------------------
%Find intersection points between shoreline and transect lines

%Create land_water array
z_transect_360_land_water=z_transect_360;
z_transect_360_land_water(z_transect_360>=shorelinelevel)=1; %All points on land are changed to 1
z_transect_360_land_water(z_transect_360<shorelinelevel)=-1; %All points in water are changed to -1

%Pre-assign array
[M,N]=size(x_transect_360);
indx_beforeSHL=ones(M,1);
indx_afterSHL=ones(M,1);
x_shorelineintersect_360=zeros(M,1);
y_shorelineintersect_360=zeros(M,1);

%Find intersection points between shoreline and transect lines
for i=1:M

    %Check if point of interest (start of transect) is in water
    if z_transect_360_land_water(i,1)==-1

        %Check if transect pass through both water and land
        if any(z_transect_360_land_water(i,:)==1)==1
            
            %Find two points one before and one after shoreline
            indx_afterSHL(i,1)=find(z_transect_360_land_water(i,:)==1,1);
            indx_beforeSHL(i,1)=indx_afterSHL(i,1)-1;

            xpoint_beforeSHL=x_transect_360(i,indx_beforeSHL(i,1)); %Assign x of point before shoreline
            ypoint_beforeSHL=y_transect_360(i,indx_beforeSHL(i,1)); %Assign y of point before shoreline
            zpoint_beforeSHL=z_transect_360(i,indx_beforeSHL(i,1)); %Assign z of point before shoreline
            xpoint_afterSHL=x_transect_360(i,indx_afterSHL(i,1)); %Assign x of point after shoreline
            ypoint_afterSHL=y_transect_360(i,indx_afterSHL(i,1)); %Assign y of point after shoreline
            zpoint_afterSHL=x_transect_360(i,indx_afterSHL(i,1)); %Assign z of point after shoreline

            %Calculate intesction point between shoreline and trensect line
            m_Line_y=(zpoint_afterSHL-zpoint_beforeSHL)/(ypoint_afterSHL-ypoint_beforeSHL); %Slope of line between point before and point afer shoreline in latitudinal direction
            m_Line_x=(zpoint_afterSHL-zpoint_beforeSHL)/(xpoint_afterSHL-xpoint_beforeSHL); %Slope of line between point before and point after shoreline in longitudinal direction
            
            x_shorelineintersect_360(i,1)=(-zpoint_beforeSHL/m_Line_x)+xpoint_beforeSHL; %x of intesction point between shoreline and trensect line
            y_shorelineintersect_360(i,1)=(-zpoint_beforeSHL/m_Line_y)+ypoint_beforeSHL; %y of intesction point between shoreline and trensect line
            
            if isnan(y_shorelineintersect_360(i,1))==1 | isnan(x_shorelineintersect_360(i,1))==1 | isinf(y_shorelineintersect_360(i,1))==1 | isinf(x_shorelineintersect_360(i,1))==1
                %Assign mean velues to intesction point between shoreline and trensect line
                x_shorelineintersect_360(i,1)=(xpoint_beforeSHL+xpoint_afterSHL)/2; %x of intesction point between shoreline and trensect line
                y_shorelineintersect_360(i,1)=(ypoint_beforeSHL+ypoint_afterSHL)/2; %y of intesction point between shoreline and trensect line
            end

        %Transect line passes completely over water and wind fetch extend to the boundary
        else
            indx_afterSHL(i,1)=M;
            indx_beforeSHL(i,1)=M;
            x_shorelineintersect_360(i,1)=x_transect_360(i,end); %Assign the farest point (on boundary) x value as x value of intersection points
            y_shorelineintersect_360(i,1)=y_transect_360(i,end); %Assign the farest point (on boundary) y value as y value of intersection points

        end

    %Check if point of interest (start of transect) is on land
    elseif z_transect_360_land_water(i,1)==1
        indx_afterSHL(i,1)=NaN;
        indx_beforeSHL(i,1)=NaN;
        x_shorelineintersect_360(i,1)=NaN; %Assign NaN as x value of intersection points
        y_shorelineintersect_360(i,1)=NaN; %Assign NaN as y value of intersection points
    end

end

%--------------------------------------------------------------------------
%Generate points on transect line starts from shoreline toward the point of interest

%Calculate using cartesian formula
if strcmp(distCalcMethod,'cart')==1
    
    %Generating points from (x1,y1) to (x2,y2)
    n_midpoints(n_midpoints<1)=1;
    x_fetch_360=linspace(x_shorelineintersect_360,x_point,n_midpoints+2);
    y_fetch_360=linspace(y_shorelineintersect_360,y_point,n_midpoints+2);

%Calculate using great circle formula
elseif strcmp(distCalcMethod,'gc')==1

    %Converte to radian
    lat1_rad=deg2rad(y_shorelineintersect_360);
    lon1_rad=deg2rad(x_shorelineintersect_360);
    lat2_rad=deg2rad(y_point);
    lon2_rad=deg2rad(x_point);

    delta_lat_21_rad=lat2_rad-lat1_rad;
    delta_lon_21_rad=lon2_rad-lon1_rad;

    %Calculating azimuth (bearing) between start and end of the line
    azimuthdir_rad=atan2(sin(delta_lon_21_rad).*cos(lat2_rad),cos(lat1_rad).*sin(lat2_rad)-sin(lat1_rad).*cos(lat2_rad).*cos(delta_lon_21_rad)); %Azimuth (bearing) in radian

    R=6371000; %Earth radius in (m), mean earth radius=6371000 m
    deltasigma=atan2(sqrt((cos(lat2_rad).*sin(delta_lon_21_rad)).^2+(cos(lat1_rad).*sin(lat2_rad)-sin(lat1_rad).*cos(lat2_rad).*cos(delta_lon_21_rad)).^2),sin(lat1_rad).*sin(lat2_rad)+cos(lat1_rad).*cos(lat2_rad).*cos(delta_lon_21_rad)); %Central angle

    %Generating lat and lon between (lon1,lat1) and (lon2,lat2), first method
    n_midpoints(n_midpoints<1)=1;
    deltasigmai=linspace(0,deltasigma,n_midpoints+2); %Central angle for each point from first point, deltasigma=d/R
    lati_rad=asin(sin(lat1_rad).*cos(deltasigmai)+cos(lat1_rad).*sin(deltasigmai).*cos(azimuthdir_rad));
    loni_rad=lon1_rad+atan2(sin(azimuthdir_rad).*sin(deltasigmai).*cos(lat1_rad),cos(deltasigmai)-sin(lat1_rad).*sin(lati_rad));

    %Converte to degree
    y_fetch_360=rad2deg(lati_rad);
    x_fetch_360=rad2deg(loni_rad);

end

%--------------------------------------------------------------------------
%Calculate distace of wind fetch points starts from shoreline toward point of interest

%Calculate distance using cartesian formula
if strcmp(distCalcMethod,'cart')==1

    dist_xy_360=sqrt((x_fetch_360-x_point).^2+(y_fetch_360-y_point).^2); %Calculate distance from (x,y) to (x(1),y(1))

%Calculate distance using Vincenty formula
elseif strcmp(distCalcMethod,'gc')==1

    %Converte to radian
    lat1_rad=deg2rad(y_fetch_360);
    lon1_rad=deg2rad(x_fetch_360);
    lat2_rad=deg2rad(y_point);
    lon2_rad=deg2rad(x_point);

    delta_lat_21_rad=lat2_rad-lat1_rad;
    delta_lon_21_rad=lon2_rad-lon1_rad;

    R=6371000; %Earth radius in (m), mean earth radius=6371000 m
    deltasigma=atan2(sqrt((cos(lat2_rad).*sin(delta_lon_21_rad)).^2+(cos(lat1_rad).*sin(lat2_rad)-sin(lat1_rad).*cos(lat2_rad).*cos(delta_lon_21_rad)).^2),sin(lat1_rad).*sin(lat2_rad)+cos(lat1_rad).*cos(lat2_rad).*cos(delta_lon_21_rad)); %Central angle
    arclen=R.*deltasigma; %Total distance of the line 
    dist_xy_360=arclen;

end

%Total wind fetch distance from shoreline to the point of interest
wind_fetch_360=dist_xy_360(:,1);

%Wind fetch distance for each point from shoreline
dist_along_fetch_360=wind_fetch_360-dist_xy_360;

%--------------------------------------------------------------------------
%Calculate z (elevation) profile along transect path from shoreline toward point of interest

%Interpolating the z (elevation) values for given (x,y) from input data
if strcmp(interpMethod,'linear')==1
    z_fetch_360=interp2(xgrid,ygrid,zgrid,x_fetch_360,y_fetch_360);

    %Replacing NaN data point resulted from default method with ones from nearest method
    if sum(isnan(z_fetch_360(:)))>0
        z_nearest=interp2(xgrid,ygrid,zgrid,x_fetch_360,y_fetch_360,'nearest');
        z_fetch_360(isnan(z_fetch_360)==1)=z_nearest(isnan(z_fetch_360)==1); 
    end

%Interpolating data into grid using nearest neighbor method
elseif strcmp(interpMethod,'nearest')==1

    %First method
    z_fetch_360=interp2(xgrid,ygrid,zgrid,x_fetch_360,y_fetch_360,'nearest');

end

z_point=z_fetch_360(1,end); %Value of z at (x_point,y_point) location

%--------------------------------------------------------------------------
%Calculate mean z (elevation)

%First method
z_mean_360=nanmean(z_fetch_360,2);

%Second method
%Calculate zdx=z(i)*dx(i)
%zdx=(0.5.*(z_fetch_360(:,1:end-1)+z_fetch_360(:,2:end))).*diff(dist_along_fetch_360,1,2);

%Calculate mean z (elevation) as z_mean = sum(z(i)*dx(i)) / total_distance
%z_mean_360=sum(zdx,2)./dist_along_fetch_360(:,end);

%for i=1:length(lat_boundary(:,1))
%    if dist_along_fetch_360(i,end)==0
%        z_mean_360(i,1)=z_fetch_360(i,end);
%    end
%end

%--------------------------------------------------------------------------
%Calculate mean slope along wind fetch

%Calculate slope_along_fetch_360=dz(i)/dx(i)
slope_along_fetch_360=diff(z_fetch_360,1,2)./diff(dist_along_fetch_360,1,2);

%Calculate mean slope
bed_slope_mean_360=nanmean(slope_along_fetch_360,2);

%--------------------------------------------------------------------------
%Calculate fetch direction from shoreline toward point of interest

%Calculate using cartesian formula
if strcmp(distCalcMethod,'cart')==1
    
    %Generating direction from (x1,y1) to (x2,y2)
    delta_x_360 = x_point-x_shorelineintersect_360;
    delta_y_360 = y_point-y_shorelineintersect_360;
    winddir_360_trig_rad=atan2(delta_y_360,delta_x_360);

    winddir_360_trig=rad2deg(winddir_360_trig_rad); %Wind direction in degree

    %Converting trigonometric direction to meteorological direction
    winddir_360=270-winddir_360_trig;

    %Add 360 to all numbers to have them all positive
    %Use mod(360) to take care of the ones larger than 360
    winddir_360=mod((winddir_360+360),360); 

%Calculate using great circle formula
elseif strcmp(distCalcMethod,'gc')==1

    %Converte to radian
    lat1_rad=deg2rad(y_shorelineintersect_360);
    lon1_rad=deg2rad(x_shorelineintersect_360);
    lat2_rad=deg2rad(y_point);
    lon2_rad=deg2rad(x_point);

    delta_lat_21_rad=lat2_rad-lat1_rad;
    delta_lon_21_rad=lon2_rad-lon1_rad;

    %Calculating azimuth (bearing) between start and end of the line
    azimuthdir_rad=atan2(sin(delta_lon_21_rad).*cos(lat2_rad),cos(lat1_rad).*sin(lat2_rad)-sin(lat1_rad).*cos(lat2_rad).*cos(delta_lon_21_rad)); %Azimuth (bearing) in radian
    azimuthdir_360=rad2deg(azimuthdir_rad); %Azimuth (bearing) in degree

    %Converting azimuth (bearing) to meteorological direction
    winddir_360=azimuthdir_360+180;
    
    %Add 360 to all numbers to have them all positive
    %Use mod(360) to take care of the ones larger than 360
    winddir_360=mod((winddir_360+360),360); 
    
end

%Subtract 360 degree from the first point and add 360 degree to the last point to have a closed loop
winddir_360(1,1)=winddir_360(1,1)-360;
winddir_360(end,1)=winddir_360(end,1)+360;

%--------------------------------------------------------------------------
%Interpolate values for given winddir

wind_fetch = interp1(winddir_360,wind_fetch_360,winddir);
z_mean = interp1(winddir_360,z_mean_360,winddir);
bed_slope_mean = interp1(winddir_360,bed_slope_mean_360,winddir);
x_fetch = interp1(winddir_360,x_fetch_360,winddir);
y_fetch = interp1(winddir_360,y_fetch_360,winddir);
z_fetch = interp1(winddir_360,z_fetch_360,winddir);
dist_along_fetch = interp1(winddir_360,dist_along_fetch_360,winddir);

%--------------------------------------------------------------------------
%Adjust z (elevation) data for water surface elevation

%For numel(waterlevel)=1
if numel(waterlevel)==1
    z_mean = z_mean-waterlevel;
    z_point = z_point-waterlevel;
    z_fetch = z_fetch-waterlevel;

%For numel(waterlevel)>1
elseif numel(waterlevel)==numel(winddir)
    z_mean = z_mean-waterlevel;
    z_point = z_point-waterlevel;
    for i=1:length(winddir(:,1))
        z_fetch(i,:) = z_fetch(i,:)-waterlevel(i,1);
    end
end

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1
    
    subplot(2,1,1)
    %Plotting domain
    [Cont,hCont]=contourf(xgrid,ygrid,zgrid,[0,0]);
    %set(hCont,'LineColor','none');
    
    %Plotting line 
    hold on
    scatter(x_point,y_point,'filled')
    sct3 = scatter(x_shorelineintersect_360,y_shorelineintersect_360,'.');
    for i=1:length(winddir(:,1))
        %plot([x_point,x_shorelineintersect_360(i,1)],[y_point,y_shorelineintersect_360(i,1)],'k')
        %scatter(x_fetch_360(i,:),y_fetch_360(i,:),[],'k','*')
        sct1 = scatter(x_fetch(i,end),y_fetch(i,end),[],'k','*');
        sct2 = scatter(x_fetch(i,1),y_fetch(i,1),[],'k','x');
        plot([x_fetch(i,1),x_fetch(i,end)],[y_fetch(i,1),y_fetch(i,end)],'k')
    end
    
    xlabel('x')
    ylabel('y')
    legend([sct1, sct2, sct3],'(x-point,y-point)','Start of the Fetch','Shoreline')
    
    subplot(2,1,2)
    %Plotting depth along a line from point 2 to point1
    for i=1:length(winddir(:,1))
        sct1 = scatter(dist_along_fetch(i,end),z_fetch(i,end),[],'k','*');
        hold on
        sct2 = scatter(dist_along_fetch(i,1),z_fetch(i,1),[],'k','x');
        plot(dist_along_fetch(i,:),z_fetch(i,:))
    end
    
    xlabel('Distance from Shoreline')
    ylabel('z')
    legend([sct1, sct2],'(x-point,y-point)','Start of the Fetch')
    box on
        
end

%--------------------------------------------------------------------------
