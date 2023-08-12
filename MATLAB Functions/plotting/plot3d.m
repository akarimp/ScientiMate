function plot3d(x, y, z, plottype, cmapcolors)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2019-02-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

plot3d
======

.. code:: MATLAB

    plot3d(x, y, z, plottype, cmapcolors)

Description
-----------

Plot x , y, and z data in 2-d/3-d contour/surface plot

Inputs
------

x
    | x data
    | Set x=[] if it is not available
    | It may be 1d or 2d array
y
    | y data
    | Set y=[] if it is not available
    | It may be 1d or 2d array
z
    | z data
    | It may be 1d or 2d array
plottype='imagesc';
    | Plot type
    | 'imagesc': 2 dimensional plot using imagesc or imshow
    | 'pcolor': 2 dimensional plot using pcolor
    | 'contour': 2 dimensional contour plot, number of contour=32
    | 'contourf': 2 dimensional filled contour plot, number of contour=32
    | 'surface': 3 dimensional surface plot 
cmapcolors='blue';
    | Colormap style
    | 'blue': blue colormap
    | 'red': red colormap
    | 'green': green colormap
    | 'yellow': yellow colormap
    | 'purple': purple colormap
    | 'brown': brown colormap
    | 'gray': gray colormap
    | 'blue_red': blue-red colormap
    | 'red_blue': red-blue colormap
    | 'blue_green': blue-green colormap
    | 'green_blue': green-blue colormap
    | 'green_yellow': green-yellow colormap
    | 'yellow_green': yellow-green colormap
    | 'red_yellow': red-yellow colormap
    | 'yellow_red': yellow-red colormap
    | 'cyclic': cyclic/oscillation colormap 
    | 'seq': sequential colormap
    | User-defined colors may be used to generate colormap
    | User-defined colors should be defined as (M,3) array in RGB color format
    | At least two colors should be defined, i.e. M should be equal or larger than 2
    | User-defined colors values should be between 0 and 255
    | Any available colormap name such as 'cool', 'winter', etc also can be used

Outputs
-------


Examples
--------

.. code:: MATLAB

    [x,y]=meshgrid(linspace(-10,10,50),linspace(-10,10,50));
    r=sqrt(x.^2+y.^2)+1e-10; %Add 1e-10 to prevent divide by 0
    z=sin(r)./r;
    plot3d(x,y,z,'surface','purple')

    [x,y]=meshgrid(linspace(-10,10,21),linspace(-10,10,21));
    z=(sin(x)+sin(y))./(x+y+1e-10); %Add 1e-10 to prevent divide by 0
    plot3d(x,y,z,'pcolor','purple')

    x(:,1)=10.*rand(1000,1);
    y(:,1)=10.*rand(1000,1);
    z=x.^2+y.^2;
    plot3d(x,y,z,'pcolor','blue')

References
----------

Colormap

* https://matplotlib.org/tutorials/colors/colormaps.html
* https://www.mathworks.com/help/matlab/ref/colormap.html
* http://colorbrewer2.org
* http://matplotlib.org/cmocean/
* http://jdherman.github.io/colormap/

Color

* http://htmlcolorcodes.com

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
        plottype='imagesc'; cmapcolors='blue';
    case 4
        cmapcolors='blue';
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
%Set plot drawing size

plotfontsize=18; %Size of plot fonts
axisfontsize=18; %Size of axis fonts


%--------------------------------------------------------------------------
%Input variable names

%if isempty(x)==1 & isempty(y)==1
%    xaxislabel='x'; %x axis label
%    yaxislabel='y'; %y axis label
%    zaxislabel=inputname(1); %z axis label
%else
%    xaxislabel=inputname(1); %x axis label
%    yaxislabel=inputname(2); %y axis label
%    zaxislabel=inputname(3); %z axis label
%end 

xaxislabel='x'; %x axis label
yaxislabel='y'; %y axis label
zaxislabel='z'; %z axis label
cbarlabel='z'; %Colorbar label


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

ncolor=32; %Number of colors to be used in colormap
gridsize=100; %Grid size in x and y directions to interpolate elevation data on them
interpMethod='nearest'; %Interpolation method ('linear','nearest')

%--------------------------------------------------------------------------
%Interpolating z data

%Check if inputs are 1-d or 2-d arrays
if Mx==1 | Nx==1

    %Calculating grid size in each direction using number grid points 
    ngridx=gridsize; %Number of grid points in each direction
    ngridy=gridsize; %Number of grid points in each direction

    gridsizex=(xmax-xmin)/(ngridx-1); %Grid size in x direction
    gridsizey=(ymax-ymin)/(ngridy-1); %Grid size in y direction

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

%--------------------------------------------------------------------------
%Colors

light_blue=[212, 230, 241];
mid_blue=[41, 128, 185];
dark_blue=[21, 67, 96];

light_red=[242, 215, 213];
mid_red=[192, 57, 43];
dark_red=[100, 30, 22];

light_green=[200, 230, 201];
mid_green=[76, 175, 80];
dark_green=[27, 94, 32];

light_yellow=[255, 249, 196];
mid_yellow=[255, 235, 59];
dark_yellow=[245, 127, 23];

light_purple=[235, 222, 240];
mid_purple=[155, 89, 182];
dark_purple=[81, 46, 95];

light_brown=[246, 221, 204];
mid_brown=[211, 84, 0];
dark_brown=[110, 44, 0];

light_gray=[250, 250, 250];
dark_gray=[1, 1, 1];


%--------------------------------------------------------------------------
%Colormaps

%zclmap
%'blue': blue colormap
if strcmp(cmapcolors,'blue')==1
    zclmap=[light_blue;mid_blue;dark_blue]./255;

%'red': red colormap
elseif strcmp(cmapcolors,'red')==1
    zclmap=[light_red;mid_red;dark_red]./255;

%'green': green colormap
elseif strcmp(cmapcolors,'green')==1
    zclmap=[light_green;mid_green;dark_green]./255;

%'yellow': yellow colormap
elseif strcmp(cmapcolors,'yellow')==1
    zclmap=[light_yellow;mid_yellow;dark_yellow]./255;

%'purple': purple colormap
elseif strcmp(cmapcolors,'purple')==1
    zclmap=[light_purple;mid_purple;dark_purple]./255;

%'brown': brown colormap, 'cooper'
elseif strcmp(cmapcolors,'brown')==1
    zclmap=[light_brown;mid_brown;dark_brown]./255;

%'gray': gray colormap, 'gray'
elseif strcmp(cmapcolors,'gray')==1
    zclmap=[light_gray;dark_gray]./255;

%'blue_red': blue-red colormap, 'cool'
elseif strcmp(cmapcolors,'blue_red')==1
    zclmap=[mid_blue;mid_red]./255;

%'red_blue': red-blue colormap, 'cool-1'
elseif strcmp(cmapcolors,'red_blue')==1
    zclmap=[mid_red;mid_blue]./255;

%'blue_green': blue-green colormap, 'winter'
elseif strcmp(cmapcolors,'blue_green')==1
    zclmap=[mid_blue;mid_green]./255;

%'green_blue': green-blue colormap, 'winter-1'
elseif strcmp(cmapcolors,'green_blue')==1
    zclmap=[mid_green;mid_blue]./255;

%'green_yellow': green-yellow colormap, 'summer'
elseif strcmp(cmapcolors,'green_yellow')==1
    zclmap=[mid_green;mid_yellow]./255;

%'yellow_green': yellow-green colormap, 'summer-1'
elseif strcmp(cmapcolors,'yellow_green')==1
    zclmap=[mid_yellow;mid_green]./255;

%'red_yellow': red-yellow colormap, 'autumn'
elseif strcmp(cmapcolors,'red_yellow')==1
    zclmap=[mid_red;mid_yellow]./255;

%'yellow_red': yellow-red colormap, 'autumn-1'
elseif strcmp(cmapcolors,'yellow_red')==1
    zclmap=[mid_yellow;mid_red]./255;

%'cyclic': cyclic/oscillation colormap, 'hsv'
elseif strcmp(cmapcolors,'cyclic')==1
    %colormap('hsv')
    %zclmap=[mid_red;light_red;light_blue;mid_blue;mid_blue;light_blue;light_red;mid_red]./255;
    zclmap=[light_red;light_blue;mid_blue;mid_red]./255;

%'seq': sequential colormap, 'line'
elseif strcmp(cmapcolors,'seq')==1
    zclmap=lines(ncolor);

%User input colormap
elseif isnumeric(cmapcolors)==1
    zclmap=cmapcolors./255; %Converting RGB to values between 0 and 1
    zclmap(zclmap<0 )=0;
    zclmap(zclmap>1 )=1;

%Pre-defined colormap such as 'cool', 'winter', etc
else
    zclmap=colormap(cmapcolors);

end


%--------------------------------------------------------------------------
%Interpolating a colormap to include ncolor zclmap

if strcmp(cmapcolors,'seq')==1
    cmap_ncolor=lines(ncolor);
else
    cmap_ncolor=zeros(ncolor,3);
    cmap_ncolor(:,1:3)=interp1(linspace(1,ncolor,length(zclmap(:,1))),zclmap(:,1:3),linspace(1,ncolor,ncolor));
end


%--------------------------------------------------------------------------
%2 dimensional plot using imagesc or imshow

if strcmp(plottype,'imagesc')==1 | strcmp(plottype,'imagesc_grid')==1
    
    imagesc([min(xgrid(:)),max(xgrid(:))],[min(ygrid(:)),max(ygrid(:))],zgrid)
    set(gca,'YDir','normal'); %Correct y axis direction

end


%--------------------------------------------------------------------------
%2 dimensional plot using pcolor

if strcmp(plottype,'pcolor')==1 | strcmp(plottype,'pcolor_grid')==1

    pcol=pcolor(xgrid,ygrid,zgrid);
    set(pcol,'edgecolor','none')

end


%--------------------------------------------------------------------------
%2 dimensional contour plot

if strcmp(plottype,'contour')==1 | strcmp(plottype,'contour_grid')==1

    zmingrid=nanmin(zgrid(:));
    zmaxgrid=nanmax(zgrid(:));
    zlevels(:,1)=linspace(zmingrid,zmaxgrid,ncolor); %Deviding z levels between zmin and zmax to ncolor
    [Cont,hCont]=contour(xgrid,ygrid,zgrid,zlevels);

end


%--------------------------------------------------------------------------
%2 dimensional contour plot

if strcmp(plottype,'contourf')==1 | strcmp(plottype,'contourf_grid')==1

    zmingrid=nanmin(zgrid(:));
    zmaxgrid=nanmax(zgrid(:));
    zlevels(:,1)=linspace(zmingrid,zmaxgrid,ncolor); %Deviding z levels between zmin and zmax to ncolor
    [Cont,hCont]=contourf(xgrid,ygrid,zgrid,zlevels);
    set(hCont,'LineColor','none');

end


%--------------------------------------------------------------------------
%Surface plot

if strcmp(plottype,'surface')==1
    
    surf1=surf(xgrid,ygrid,zgrid);
    %set(surf1,'xscale','log','zscale','log')
    set(surf1,'FaceColor','interp','EdgeColor','none')

end


%--------------------------------------------------------------------------
%Set plot properties

%Draw outline around a plot
box on

%Set axis limits
xlim([nanmin(xgrid(:)),nanmax(xgrid(:))])
ylim([nanmin(ygrid(:)),nanmax(ygrid(:))])
if strcmp(plottype,'surface')==1
    zlim([nanmin(zgrid(:)),nanmax(zgrid(:))])
end

%Use colormap
colormap(cmap_ncolor)

%Plot colorbar
cbar=colorbar;
set(cbar,'fontsize',axisfontsize)
xlabel(cbar,cbarlabel,'fontsize',axisfontsize)
%c=max(abs([zmin,zmax]));
%caxis([-c,c])
%caxis([zmin,zmax])

%Plot grid lines
if strcmp(plottype,'imagesc_grid')==1 | strcmp(plottype,'pcolor_grid')==1 ...
    | strcmp(plottype,'contour_grid')==1 | strcmp(plottype,'contourf_grid')==1
    grid on
end

%Set figure font size
    set(gca,'fontsize',plotfontsize)

%Set axis label
xlabel(xaxislabel,'fontsize',axisfontsize)
ylabel(yaxislabel,'fontsize',axisfontsize)
if strcmp(plottype,'surface')==1
    zlabel(zaxislabel,'fontsize',axisfontsize)
end
    

%--------------------------------------------------------------------------
