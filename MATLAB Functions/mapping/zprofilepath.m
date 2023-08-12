function [z, zmean, distxy] = zprofilepath(xgrid, ygrid, zgrid, x, y, distCalcMethod, CalcMethod, dispout)
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

zprofilepath
============

.. code:: MATLAB

    [z, zmean, distxy] = zprofilepath(xgrid, ygrid, zgrid, x, y, distCalcMethod, CalcMethod, dispout)

Description
-----------

Calculate z (elevation, ...) profile along a path over a given 2d x-y domain (map, image, ...)

Inputs
------

xgrid
    x (longitude, pixel, ...) data as a [M*N] array
ygrid
    y (latitude, pixel, ...) data as a [M*N] array
zgrid
    z (elevation, ...) data as a [M*N] array
x
    x (longitude, pixel, ...) of the points along the line as (x,y)
y
    | y (latitude, pixel, ...) of the points along the line as (x,y)
    | If input data are latitude and longitude, they should be in Degree
distCalcMethod='cart';
    | Distance calculation method 
    | 'cart': Distances are calculated on cartesian coordinate
    | 'gc': Distances are calculated on Great Circle based on Vincenty formula, Vincenty (1975)
    | Earth radius coonsidered as mean earth radius=6371000 m
CalcMethod='nearest';
    | Interpolation method 
    | 'linear': Use default or 'linear' method to interpolate
    | 'nearest': Use nearest neighbor method to interpolate
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

z
    z (elevation, ...) data along a path at given points (x,y)
zmean
    Weighted mean of z (elevation) along a line calculated az zmean=1/(total_distance)*(sum(z(i)*dx(i)))
distxy
    | Distance at each points of (x,y) from the first point of (x(1),y(1))
    | If input data are latitude and longitude in Degree, distxy is in m

Examples
--------

.. code:: MATLAB

    [xgrid,ygrid]=meshgrid(linspace(-5,5,100),linspace(-5,5,100));
    zgrid=sin(xgrid.^2+ygrid.^2)./(xgrid.^2+ygrid.^2);
    x1=-3;
    y1=-3;
    x2=3;
    y2=3;
    x(:,1)=linspace(x1,x2,1000);
    y(:,1)=linspace(y1,y2,1000);
    [z,zmean,distxy]=zprofilepath(xgrid,ygrid,zgrid,x,y,'cart','nearest','yes');

References
----------

Vincenty, T. (1975). 
Direct and inverse solutions of geodesics on the ellipsoid with application of nested equations. 
Survey review, 23(176), 88-93.

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
        distCalcMethod='cart'; CalcMethod='nearest'; dispout='no';
    case 6
        CalcMethod='nearest'; dispout='no';
    case 7
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(x)==1
    x=x';
end

if isrow(y)==1
    y=y';
end

%--------------------------------------------------------------------------
%Calculating z (elevation) profile

%Interpolating the z (elevation) values for given (x,y) from input data
if strcmp(CalcMethod,'linear')==1
    z=interp2(xgrid,ygrid,zgrid,x,y);

    %Replacing NaN data point resulted from default method with ones from nearest method
    if sum(isnan(z(:)))>0
        znearest=interp2(xgrid,ygrid,zgrid,x,y,'nearest');
        z(isnan(z)==1)=znearest(isnan(z)==1); 
    end

%Interpolating data into grid using nearest neighbor method
elseif strcmp(CalcMethod,'nearest')==1

    %First method
    z=interp2(xgrid,ygrid,zgrid,x,y,'nearest');

    %Second method
    %Mapping values of x and y to associated index
    %[M,N]=size(xgrid);
    %xmap(:,1)=interp1([nanmin(xgrid(:));nanmax(xgrid(:))],[1;N],x,'extrap');
    %ymap(:,1)=interp1([nanmin(ygrid(:));nanmax(ygrid(:))],[1;M],y,'extrap');
    %xmap(xmap<1)=1;
    %xmap(xmap>N)=N;
    %ymap(ymap<1)=1;
    %ymap(ymap>M)=M;
    %z=diag(zgrid(fix(ymap),fix(xmap)));

end

%Calculating distance using cartesian formula
if strcmp(distCalcMethod,'cart')==1

    distxy(:,1)=sqrt((x-x(1,1)).^2+(y-y(1,1)).^2); %Calculating distance from (x,y) to (x(1),y(1))

%Calculating distance using Vincenty formula
elseif strcmp(distCalcMethod,'gc')==1

    %Converting to radian
    lat1rad=deg2rad(y(1,1));
    lon1rad=deg2rad(x(1,1));
    lat2rad=deg2rad(y);
    lon2rad=deg2rad(x);

    deltalatrad21=lat2rad-lat1rad;
    deltalonrad21=lon2rad-lon1rad;

    R=6371000; %Earth radius in (m), mean earth radius=6371000 m
    deltasigma=atan2(sqrt((cos(lat2rad).*sin(deltalonrad21)).^2+(cos(lat1rad).*sin(lat2rad)-sin(lat1rad).*cos(lat2rad).*cos(deltalonrad21)).^2),sin(lat1rad).*sin(lat2rad)+cos(lat1rad).*cos(lat2rad).*cos(deltalonrad21)); %Central angle
    arclen=R.*deltasigma; %Total distance of the line 
    distxy(:,1)=arclen;

end

%--------------------------------------------------------------------------
%Calculating mean z (elevation)

%Calculating zdx=z(i)*dx(i)
zdx=(0.5.*(z(1:end-1,1)+z(2:end,1))).*diff(distxy);

if distxy(end,1)~=0
    %Calculating mean z (elevation) as zmean=1/(total_distance)*(sum(z(i)*dx(i)))
    zmean=sum(zdx)/distxy(end,1);
else
    zmean=z(end,1);
end

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1
    
    subplot(2,1,1)
    %Plotting domain
    imagesc(xgrid,ygrid,zgrid);
    set(gca,'YDir','normal'); %Correct y axis direction
    colorbar

    %Plotting line 
    hold on
    plot(x,y);

    xlabel('x')
    ylabel('y')
    
    subplot(2,1,2)
    %Plotting depth along a line from point 2 to point1
    plot(distxy,z);
    xlabel('distance')
    ylabel('z')
    
end

%--------------------------------------------------------------------------
