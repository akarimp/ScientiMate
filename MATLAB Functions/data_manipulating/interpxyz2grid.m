function [xgrid, ygrid, zgrid] = interpxyz2grid(x, y, z, gridsize, gridsizetype, xmin, xmax, ymin, ymax, zmin, zmax, RetainRatio, interpMethod, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-12-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

interpxyz2grid
==============

.. code:: MATLAB

    [xgrid, ygrid, zgrid] = interpxyz2grid(x, y, z, gridsize, gridsizetype, xmin, xmax, ymin, ymax, zmin, zmax, RetainRatio, interpMethod, dispout)

Description
-----------

Interpolate x (longitude), y (latitude) and z (elevation) data into a defined mesh

Inputs
------

x
    x (longitude) data extracted from xyz file
y
    y (latitude) data extracted from xyz file
z
    z (elevation) data extracted from xyz file
gridsize=100;
    | Grid size in x (longitude) and y (latitude) directions to interpolate elevation data on them
    | if gridsizetype='length' then gridsize is a distance between grid points
    | if gridsizetype='points' then gridsize is number of grid points in each direction
gridsizetype='points';
    | Grid size type 
    | 'points': gridsize is considered as number of grid points in each direction
    | 'length': gridsize is considered as length between grid points
xmin=nanmin(x);
    Minimum x (longitude) of domain to be interpolated
xmax=nanmax(x);
    Maximum x (longitude) of domain to be interpolated
ymin=nanmin(y);
    Minimum y (latitude) of domain to be interpolated
ymax=nanmax(y);
    Maximum y (latitude) of domain to be interpolated
zmin=nanmin(z);
    | Minimum z (elevation) of domain to be interpolated
    | All z<zmin would be set to zmin
zmax=nanmax(z);
    | Maximum z (elevation) of domain to be interpolated
    | All z>zmax would be set to zmax
RetainRatio='all';
    | Define to down sample input data or not 
    | 'all': data are not down sampled
    | value between 0 and 1: percentage of retaining data
    | RetainRatio=0.8; : 80% of data are retained
interpMethod='nearest';
    | Interpolation method 
    | 'linear': Use default or 'linear' method to interpolate
    | 'nearest': Use nearest neighbor method to interpolate
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

xgrid
    Interpolated x (longitude) data on defined mesh
ygrid
    Interpolated y (latitude) data on defined mesh
zgrid
    Interpolated z (elevation) data on defined mesh

Examples
--------

.. code:: MATLAB

    x(:,1)=10.*rand(1000,1);
    y(:,1)=10.*rand(1000,1);
    z=x.^2+y.^2;
    [xgrid,ygrid,zgrid]=interpxyz2grid(x,y,z,100,'points',nanmin(x),nanmax(x),nanmin(y),nanmax(y),nanmin(z),nanmax(z),'all','nearest','yes');

    x(:,1)=(-90-(-91)).*rand(1000,1)+(-91);
    y(:,1)=(31-(30)).*rand(1000,1)+(30);
    z=x.^2+y.^2;
    [xgrid,ygrid,zgrid]=interpxyz2grid(x,y,z,0.005,'length',nanmin(x),nanmax(x),nanmin(y),nanmax(y),nanmin(z),nanmax(z),'all','linear','yes');

References
----------

Geospatial data

* https://www.mathworks.com/help/map/finding-geospatial-data.html
* https://maps.ngdc.noaa.gov/viewers/wcs-client/
* https://www.ngdc.noaa.gov/mgg/global/global.html
* https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/
* https://www.ngdc.noaa.gov/mgg/image/2minrelief.html
* https://www.ngdc.noaa.gov/mgg/coastal/crm.html
* https://viewer.nationalmap.gov/launch/
* https://earthexplorer.usgs.gov
* http://www.shadedrelief.com/cleantopo2/index.html

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
    case 3
        gridsize=100; gridsizetype='points'; xmin=nanmin(x(:)); xmax=nanmax(x(:)); ymin=nanmin(y(:)); ymax=nanmax(y(:)); zmin=nanmin(z(:)); zmax=nanmax(z(:)); RetainRatio='all'; interpMethod='nearest'; dispout='no';
    case 4
        gridsizetype='points'; xmin=nanmin(x(:)); xmax=nanmax(x(:)); ymin=nanmin(y(:)); ymax=nanmax(y(:)); zmin=nanmin(z(:)); zmax=nanmax(z(:)); RetainRatio='all'; interpMethod='nearest'; dispout='no';
    case 5
        xmin=nanmin(x(:)); xmax=nanmax(x(:)); ymin=nanmin(y(:)); ymax=nanmax(y(:)); zmin=nanmin(z(:)); zmax=nanmax(z(:)); RetainRatio='all'; interpMethod='nearest'; dispout='no';
    case 6
        xmax=nanmax(x(:)); ymin=nanmin(y(:)); ymax=nanmax(y(:)); zmin=nanmin(z(:)); zmax=nanmax(z(:)); RetainRatio='all'; interpMethod='nearest'; dispout='no';
    case 7
        ymin=nanmin(y(:)); ymax=nanmax(y(:)); zmin=nanmin(z(:)); zmax=nanmax(z(:)); RetainRatio='all'; interpMethod='nearest'; dispout='no';
    case 8
        ymax=nanmax(y(:)); zmin=nanmin(z(:)); zmax=nanmax(z(:)); RetainRatio='all'; interpMethod='nearest'; dispout='no';
    case 9
        zmin=nanmin(z(:)); zmax=nanmax(z(:)); RetainRatio='all'; interpMethod='nearest'; dispout='no';
    case 10
        zmax=nanmax(z(:)); RetainRatio='all'; interpMethod='nearest'; dispout='no';
    case 11
         RetainRatio='all'; interpMethod='nearest'; dispout='no';
    case 12
         interpMethod='nearest'; dispout='no';
    case 13
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

if isrow(z)==1
    z=z';
end

%--------------------------------------------------------------------------
%Defining domain area and associated data

%Assigning new x, y and z
xdata=x; %xdata
ydata=y; %ydata
zdata=z; %zdata

%Defining domain
%domain='all';
%    | Define a domain to be interpolated
%    | 'all': all data are interpolated 
%    | 'domain': only data within a defined domain are interpolated
if (xmin==nanmin(x(:)) & xmax==nanmax(x(:)) & ymin==nanmin(y(:)) & ymax==nanmax(y(:)) & zmin==nanmin(z(:)) & zmax==nanmax(z(:)))
    domain='all';
else
    domain='domain';
end

%Defining domain boundaries
if strcmp(domain,'all')==1

    xmin=nanmin(x(:)); 
    xmax=nanmax(x(:)); 
    ymin=nanmin(y(:)); 
    ymax=nanmax(y(:)); 
    zmin=nanmin(z(:)); 
    zmax=nanmax(z(:));

%Retaining data within a defined domain
elseif strcmp(domain,'domain')==1
    
    %Finding a location data that are within [xmin:xmax] range
    xIndx=find(x>=xmin & x<=xmax);
    
    %Retaining data that are within [xmin:xmax] range
    x1=xdata(xIndx,1); %x data
    y1=ydata(xIndx,1); %y data
    z1=zdata(xIndx,1); %z data
    
    %Finding a location data that are within [ymin:ymax] range
    yIndx=find(y1>=ymin & y1<=ymax);
    
    %Retaining data that are within [ymin:ymax] range
    x2=x1(yIndx,1); %xdata
    y2=y1(yIndx,1); %ydata
    z2=z1(yIndx,1); %zdata
    
    %Finding a location data that are within [zmin:zmax] range
    zIndx=find(z2>=zmin & z2<=zmax);
    
    %Retaining data that are within [zmin:zmax] range
    x3=x2(zIndx,1); %xdata
    y3=y2(zIndx,1); %ydata
    z3=z2(zIndx,1); %zdata

    %Assigning new x, y and z
    xdata=x3; %xdata
    ydata=y3; %ydata
    zdata=z3; %zdata

end

%--------------------------------------------------------------------------
%Down sampling the data

if strcmp(RetainRatio,'all')==1
    %Do nothing
    dummy=1;

elseif isnumeric(RetainRatio)==1

    %Defining index number to down sample data
    RetainRatio(RetainRatio<=0 | RetainRatio>1)=1;
    [M,N]=size(x);
    dsIndx(:,1)=fix(linspace(1,M,fix(M*RetainRatio)));

    x4=xdata; %xdata
    y4=ydata; %ydata
    z4=zdata; %zdata

    %Retaining down sampled data
    xdata=x4(dsIndx,1); %xdata
    ydata=y4(dsIndx,1); %ydata
    zdata=z4(dsIndx,1); %zdata

end

%--------------------------------------------------------------------------
%Interpolating elevation data

%Calculating grid size in each direction using grid length
if strcmp(gridsizetype,'length')==1

    %Checking x and y range to generate correct number of grid points
    if mod((xmax-xmin),gridsize)~=0
        xmax=xmax-mod((xmax-xmin),gridsize);
    end

    if mod((ymax-ymin),gridsize)~=0
        ymax=ymax-mod((ymax-ymin),gridsize);
    end

    gridsizex=gridsize; %Grid size in lat (x) direction
    gridsizey=gridsize; %Grid size in lon (y) direction

    ngridx=((xmax-xmin)/gridsizex)+1; %Number of grid points in each direction
    ngridy=((ymax-ymin)/gridsizey)+1; %Number of grid points in each direction

    %Generating grid
    [xgrid,ygrid]=meshgrid([xmin:gridsizex:xmax],[ymin:gridsizey:ymax]);

%Calculating grid size in each direction using number grid points 
elseif strcmp(gridsizetype,'points')==1

    ngridx=gridsize; %Number of grid points in each direction
    ngridy=gridsize; %Number of grid points in each direction

    gridsizex=(xmax-xmin)/(ngridx-1); %Grid size in x direction
    gridsizey=(ymax-ymin)/(ngridy-1); %Grid size in y direction

    %Generating grid
    [xgrid,ygrid]=meshgrid(linspace(xmin,xmax,ngridx),linspace(ymin,ymax,ngridy));

end

%Interpolating data into grid
[xgrid,ygrid,zgrid]=griddata(xdata,ydata,zdata,xgrid,ygrid,interpMethod);

%Replacing NaN data point resulted from interpolation by nearest point
if sum(isnan(zgrid(:)))>0
    [xgrid,ygrid,zgridnearest]=griddata(xdata,ydata,zdata,xgrid,ygrid,'nearest');
    zgrid(isnan(zgrid)==1)=zgridnearest(isnan(zgrid)==1);
end


%Makeing sure all z data are in range of [zmin,zmax]
zgrid(zgrid<zmin)=zmin;
zgrid(zgrid>zmax)=zmax;

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1
    
    %[Cont,hCont]=contourf(xgrid,ygrid,zgrid);
    %set(hCont,'LineColor','none');

    %pcol=pcolor(xgrid,ygrid,zgrid);
    %set(pcol,'edgecolor','none')

    imagesc([min(xgrid(:)),max(xgrid(:))],[min(ygrid(:)),max(ygrid(:))],zgrid)
    set(gca,'YDir','normal'); %Correct y axis direction

    xlabel('x')
    ylabel('y')
    colorbar
    
end

%--------------------------------------------------------------------------
