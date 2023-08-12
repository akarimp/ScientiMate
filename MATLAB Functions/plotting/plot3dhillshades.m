function plot3dhillshades(x, y, z, plottype)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2020-02-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

plot3dhillshades
================

.. code:: MATLAB

    plot3dhillshades(x, y, z, plottype)

Description
-----------

Plot hillshades (shaded relief) of x (longitude), y (latitude) and z (elevation) data

Inputs
------

x
    | x (longitude) data extracted from xyz file
    | Set x=[] if it is not available
    | It may be 1d or 2d array
y
    | y (latitude) data extracted from xyz file
    | Set y=[] if it is not available
    | It may be 1d or 2d array
z
    z (elevation) data extracted from xyz file
plottype='imagesc';
    | Plot type
    | 'imagesc': 2 dimensional plot using imagesc or imshow
    | 'pcolor': 2 dimensional plot using pcolor
    | 'contour': 2 dimensional contour plot, number of contour=ncolor
    | 'surface': 3 dimensional surface plot 

Outputs
-------


Examples
--------

.. code:: MATLAB

    [x,y]=meshgrid(linspace(-10,10,50),linspace(-10,10,50));
    r=sqrt(x.^2+y.^2)+1e-10; %Add 1e-10 to prevent divide by 0
    z=sin(r)./r;
    plot3dhillshades(x,y,z,'imagesc')

    [x,y]=meshgrid(linspace(-10,10,21),linspace(-10,10,21));
    z=(sin(x)+sin(y))./(x+y+1e-10); %Add 1e-10 to prevent divide by 0
    plot3dhillshades(x,y,z,'imagesc')

    x(:,1)=10.*rand(1000,1);
    y(:,1)=10.*rand(1000,1);
    z=x.^2+y.^2;
    plot3dhillshades(x,y,z,'imagesc')


References
----------

* https://pro.arcgis.com/en/pro-app/tool-reference/3d-analyst/how-hillshade-works.htm
* https://pro.arcgis.com/en/pro-app/tool-reference/3d-analyst/applying-a-z-factor.htm
* http://mike.teczno.com/img/hillshade.py
* http://geospatialpython.com/2013/12/python-and-elevation-data-creating.html
* https://github.com/ThomasLecocq/geophysique.be/blob/master/2014-02-25%20Shaded%20Relief%20Map%20in%20Python.ipynb
* https://www.geophysique.be/2014/02/25/shaded-relief-map-in-python/
* https://www.mathworks.com/matlabcentral/fileexchange/14863-hillshade
* http://www.reliefshading.com/

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
    plottype='imagesc';
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
%Check input data

%Create x if x=[]
if isempty(x)==1
    if ndims(z)==1
        x(:,1)=linspace(0,1,numel(z));
    else
        [x,~]=meshgrid(linspace(0,1,length(z(1,:))),linspace(0,1,length(z(:,1))));
    end
end

%Create y if y=[]
if isempty(y)==1
    if ndims(z)==1
        y(:,1)=linspace(0,1,numel(z));
    else
        [~,y]=meshgrid(linspace(0,1,length(z(1,:))),linspace(0,1,length(z(:,1))));
    end
end

%Input data shape
[Mx,Nx]=size(x);
[My,Ny]=size(y);
[Mz,Nz]=size(z);


%--------------------------------------------------------------------------

ncolor=32; %Number of colors to be used in colormap
plotfontsize=18; %Size of plot fonts
axisfontsize=18; %Size of axis fonts
xaxislabel='x'; %x axis label
yaxislabel='y'; %y axis label
cbarlabel='z'; %Colorbar label
dispcolorbar='no'; %Define to display colorbar or not ('yes': display, 'no': not display)
dispgrid='no'; %Define to display grid lines or not ('yes': display, 'no': not display)
gridsize=100; %Grid size in x (longitude) and y (latitude) directions to interpolate elevation data on them
gridsizetype='points'; %Grid size type 
xmin=nanmin(x(:)); %Minimum x (longitude) of domain to be plotted
xmax=nanmax(x(:)); %Maximum x (longitude) of domain to be plotted
ymin=nanmin(y(:)); %Minimum y (latitude) of domain to be plotted
ymax=nanmax(y(:)); %Maximum y (latitude) of domain to be plotted
zmin=nanmin(z(:)); %Minimum z (elevation) of domain to be plotted
zmax=nanmax(z(:)); %Maximum z (elevation) of domain to be plotted
interpMethod='nearest'; %Interpolation method 

%--------------------------------------------------------------------------
%Defining domain area and associated data

%Assigning new x, y and z
xdata=x; %xdata
ydata=y; %ydata
zdata=z; %zdata

%Defining domain boundaries
xmin=nanmin(x(:));
xmax=nanmax(x(:));
ymin=nanmin(y(:));
ymax=nanmax(y(:));
zmin=nanmin(z(:));
zmax=nanmax(z(:));


%--------------------------------------------------------------------------
%Interpolating elevation data

%Check if inputs are 1-d or 2-d arrays
if Mx==1 | Nx==1

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

    %Calculating grid size in each direction using number grid points 
    elseif strcmp(gridsizetype,'points')==1

        ngridx=gridsize; %Number of grid points in each direction
        ngridy=gridsize; %Number of grid points in each direction

        gridsizex=(xmax-xmin)/(ngridx-1); %Grid size in x direction
        gridsizey=(ymax-ymin)/(ngridy-1); %Grid size in y direction

    end

    %Generating grid
    %[xgrid,ygrid]=meshgrid(xmin:gridsizex:xmax,ymin:gridsizey:ymax);
    [xgrid,ygrid]=meshgrid(linspace(xmin,xmax,ngridx),linspace(ymin,ymax,ngridy));

    %Interpolating data into grid
    [xgrid,ygrid,zgrid]=griddata(xdata,ydata,zdata,xgrid,ygrid,interpMethod);

    %Replacing NaN data point resulted from interpolation by nearest point
    if sum(isnan(zgrid(:)))>0
        [xgrid,ygrid,zgridnearest]=griddata(xdata,ydata,zdata,xgrid,ygrid,'nearest');
        zgrid(isnan(zgrid)==1)=zgridnearest(isnan(zgrid)==1); 
    end

else
    xgrid=x;
    ygrid=y;
    zgrid=z;

end

%Makeing sure all z data are in range of [zmin,zmax]
zgrid(zgrid<zmin)=zmin;
zgrid(zgrid>zmax)=zmax;

%-----------------------------------------------------------------------------
%Calculate Hillshade using ArcGIS method

%Cell size
dx=xgrid(1,2)-xgrid(1,1);
dy=ygrid(2,1)-ygrid(1,1);

%Altitude and Azimuth of illumination source
altitude_deg=45; %Angle of illumination source above the horizon (default 45)
azimuth_deg=315; %Angular direction of the sun, measured from north in clockwise from 0 to 360 (default 315), Exp: azimuth of 90 degrees is east

%Calculate z-factor
%z-factor is a conversion factor that adjusts units of measure for vertical (elevation) units 
%when they are different from horizontal coordinate (x,y) units of the input surface
%If (x,y) units are decimal degrees and z-units are meters, then z-factor is:
%known_latitude=[0,10,20,30,40,50,60,70,80];
%known_z_factor=[0.00000898,0.00000912,0.00000956,0.00001036,0.00001171,0.00001395,0.00001792,0.00002619,0.00005156]; %z-factor if z is in meter
%z_factor=interp1(known_latitude,known_z_factor,ygrid,'linear',1);
z_factor=1;

%Compute illumination angle
zenith_deg=90-altitude_deg; %Angle from zenith point (directly overhead) to the direction of illumination source
zenith_rad=deg2rad(zenith_deg); %Convert to radian

%Compute illumination direction
azimuth_math=360-azimuth_deg+90; %Convert from geographic unit (compass direction) to a mathematic unit (right angle)
azimuth_math(azimuth_math>360)=azimuth_math-360; %Make sure azimuth_math is equal or less than 360
azimuth_rad=deg2rad(azimuth_math); %Convert to radian

%Compute slope and aspect
[dz_dx,dz_dy]=gradient(zgrid,dx,dy); %Rate of change in the x and y directions
slope_rad=atan(z_factor.*sqrt(dz_dx.^2+dz_dy.^2)); %Steepest downhill descent from each cell in the surface
aspect_rad=atan2(dz_dy,-dz_dx); %Direction that steepest downslope direction is facing
aspect_rad(aspect_rad<0)=2*pi+aspect_rad(aspect_rad<0); %Make sure aspect_rad is equal or larger than 0
aspect_rad(dz_dx==0 & dz_dy>0)=pi/2; %Check for dz_dx==0 and dz_dy>0
aspect_rad(dz_dx==0 & dz_dy<0)=2*pi-pi/2; %Check for dz_dx==0 and dz_dy<0

%Calculate hillshade
hillshade=255.*((cos(zenith_rad).*cos(slope_rad))+(sin(zenith_rad).*sin(slope_rad).*cos(azimuth_rad-aspect_rad)));
hillshade(hillshade<0)=0; %Check for hillshade<0

%-----------------------------------------------------------------------------
%Example

%zgrid=[[2450,2461,2483];[2452,2461,2483];[2447,2455,2477]];
%dx=5;
%dy=5;
%altitude_deg=45;
%azimuth_deg=315;
%z_factor=1;
%zenith_deg=90-altitude_deg; %Answer: 45
%zenith_rad=deg2rad(zenith_deg); %Answer: 0.7857142857
%azimuth_math=360-azimuth_deg+90; %Answer: 135
%azimuth_math(azimuth_math>360)=azimuth_math-360; %Answer: 135
%azimuth_rad=deg2rad(azimuth_math); %Answer: 2.3571428571
%[dz_dx,dz_dy]=gradient(zgrid,dx,dy); %Answer for center point: dz_dx=3.125, dz_dy=-0.525
%slope_rad=atan(z_factor.*sqrt(dz_dx.^2+dz_dy.^2)); %Answer for center point: 1.26511
%aspect_rad=atan2(dz_dy,-dz_dx); %Answer for center point: -2.9751469600412
%aspect_rad(aspect_rad<0)=2*pi+aspect_rad(aspect_rad<0); %Answer for center point: 3.310567
%aspect_rad(dz_dx==0 & dz_dy>0)=pi/2; %Answer for center point: 3.310567
%aspect_rad(dz_dx==0 & dz_dy<0)=2*pi-pi/2; %Answer for center point: 3.310567
%hillshade=255.*((cos(zenith_rad).*cos(slope_rad))+(sin(zenith_rad).*sin(slope_rad).*cos(azimuth_rad-aspect_rad))); %Answer for center point: 153.82
%hillshade(hillshade<0)=0; %Answer for center point: 153.82

%--------------------------------------------------------------------------
%Displaying results

if strcmp(plottype,'imagesc')==1 %2 dimensional plot using imagesc or imshow
    
    imagesc([min(xgrid(:)),max(xgrid(:))],[min(ygrid(:)),max(ygrid(:))],hillshade)
    set(gca,'YDir','normal'); %Correct y axis direction

elseif strcmp(plottype,'pcolor')==1 %2 dimensional plot using pcolor

    pcol=pcolor(xgrid,ygrid,hillshade);
    set(pcol,'edgecolor','none')

elseif strcmp(plottype,'contour')==1 %2 dimensional contour plot

    zmingrid=nanmin(hillshade(:));
    zmaxgrid=nanmax(hillshade(:));
    zlevels(:,1)=linspace(zmingrid,zmaxgrid,ncolor); %Deviding z levels between zmin and zmax to ncolor
    [Cont,hCont]=contourf(xgrid,ygrid,hillshade,zlevels);
    set(hCont,'LineColor','none');

elseif strcmp(plottype,'surface')==1 %Surface plot
    
    surf1=surf(xgrid,ygrid,hillshade);
    %set(surf1,'xscale','log','zscale','log')
    set(surf1,'FaceColor','interp','EdgeColor','none')

end

%Setting plot properties
if (strcmp(plottype,'imagesc')==1 | strcmp(plottype,'pcolor')==1 | strcmp(plottype,'contour')==1 | strcmp(plottype,'surface')==1)

    %Setting axis limits
    xlim([nanmin(xgrid(:)),nanmax(xgrid(:))])
    ylim([nanmin(ygrid(:)),nanmax(ygrid(:))])

    %Setting colormap
    colormap(gray)

    %Plotting colorbar
    if strcmp(dispcolorbar,'yes')==1
        cbar=colorbar;
        set(cbar,'fontsize',axisfontsize)
        xlabel(cbar,cbarlabel,'fontsize',axisfontsize)
        %c=max(abs([zmin,zmax]));
        %caxis([-c,c])
        %caxis([zmin,zmax])
    end

    %Plotting grid lines
    if strcmp(dispgrid,'yes')==1
        grid
        set(gca,'xtick',linspace(nanmin(xgrid(:)),nanmax(xgrid(:)),13));
        set(gca,'ytick',linspace(nanmin(ygrid(:)),nanmax(ygrid(:)),13));
        %set(gca,'xtick',linspace(-180,180,13));
        %set(gca,'ytick',linspace(-90,90,13));
    elseif isnumeric(dispgrid)==1
        grid
        set(gca,'xtick',linspace(nanmin(xgrid(:)),nanmax(xgrid(:)),dispgrid));
        set(gca,'ytick',linspace(nanmin(ygrid(:)),nanmax(ygrid(:)),dispgrid));    
    end

    %Setting figure font size
    set(gca,'fontsize',plotfontsize)

    %Setting axis label
    xlabel(xaxislabel,'fontsize',axisfontsize)
    ylabel(yaxislabel,'fontsize',axisfontsize)
    
end

%--------------------------------------------------------------------------
