function plot3ddem(x, y, z, plottype, cmapcolors)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2022-05-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

plot3ddem
=========

.. code:: MATLAB

    plot3ddem(x, y, z, plottype, cmapcolors)

Description
-----------

Plot x (longitude), y (latitude) and z (elevation) data into a defined mesh

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
cmapcolors='topocmap';
    | Colormap for z data
    | Topographic (Water/Land) colormaps:
    | 'topocmap': colormap developed by Arash Karimpour
    | 'topocmaprelief': colormap developed by Arash Karimpour
    | 'topocmapocean': colormap developed by Arash Karimpour
    | 'topocmapearth': colormap developed by Arash Karimpour
    | 'blueocean': colormap developed by Arash Karimpour
    | 'blueoceansea': colormap developed by Arash Karimpour
    | 'greenearth': colormap developed by Arash Karimpour
    | 'greenearthland': colormap developed by Arash Karimpour
    | 'blgrtopocmap': colormap developed by Arash Karimpour
    | 'blrdtopocmap': colormap developed by Arash Karimpour
    | 'grayearth': colormap developed by Arash Karimpour
    | 'etopo1': ETOPO1 colormap, https://www.ngdc.noaa.gov/mgg/global/global.html
    | 'gmtglobe': GMT_globe colormap, https://www.giss.nasa.gov/tools/panoply/colorbars/
    | 'gmtrelief': GMT_relief colormap, https://www.giss.nasa.gov/tools/panoply/colorbars/
    | 'aendekerk': Colormap from  Florian Aendekerk, http://www.mathworks.com/matlabcentral/fileexchange/63590-landseacolormap-m-
    | Any other available color map such as 'cool', 'winter', etc can be used
    | Colormap can be defined by user as [n*3] array in RGB color format between 0 and 255

Outputs
-------


Examples
--------

.. code:: MATLAB

    [x,y]=meshgrid(linspace(-10,10,50),linspace(-10,10,50));
    r=sqrt(x.^2+y.^2)+1e-10; %Add 1e-10 to prevent divide by 0
    z=sin(r)./r;
    plot3ddem(x,y,z,'pcolor','topocmap')

    [x,y]=meshgrid(linspace(-10,10,21),linspace(-10,10,21));
    z=(sin(x)+sin(y))./(x+y+1e-10); %Add 1e-10 to prevent divide by 0
    plot3ddem(x,y,z,'pcolor','topocmap')

    x(:,1)=10.*rand(1000,1);
    y(:,1)=10.*rand(1000,1);
    z=x.^2+y.^2;
    plot3ddem(x,y,z,'pcolor','topocmap')


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

Colormap

* http://colorbrewer2.org
* http://matplotlib.org/cmocean/
* https://matplotlib.org/users/colormaps.html
* http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml
* https://www.giss.nasa.gov/tools/panoply/colorbars/
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
    plottype='imagesc'; cmapcolors='topocmap';
case 4
    cmapcolors='topocmap';
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
waterlandcmapratio='auto'; %put z=0 at center of colormap
plotfontsize=18; %Size of plot fonts
axisfontsize=18; %Size of axis fonts
xaxislabel='x'; %x axis label
yaxislabel='y'; %y axis label
cbarlabel='z'; %Colorbar label
dispcolorbar='yes'; %Define to display colorbar or not ('yes': display, 'no': not display)
dispgrid='no'; %Define to display grid lines or not ('yes': display, 'no': not display)
gridsize=100; %Grid size in x (longitude) and y (latitude) directions to interpolate elevation data on them
gridsizetype='points'; %Grid size type 
xmin=nanmin(x(:)); %Minimum x (longitude) of domain to be plotted
xmax=nanmax(x(:)); %Maximum x (longitude) of domain to be plotted
ymin=nanmin(y(:)); %Minimum y (latitude) of domain to be plotted
ymax=nanmax(y(:)); %Maximum y (latitude) of domain to be plotted
zmin=nanmin(z(:)); %Minimum z (elevation) of domain to be plotted
zmax=nanmax(z(:)); %Maximum z (elevation) of domain to be plotted
RetainRatio='all'; %Define to down sample input data or not 
interpMethod='nearest'; %Interpolation method 

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

%--------------------------------------------------------------------------
%Colormaps

%'topocmap' colormap
%Colors from: http://htmlcolorcodes.com
if strcmp(cmapcolors,'topocmap')==1

    %Water colormap
    topocmap_water=[[21, 67, 96];[27, 79, 114];[40, 116, 166];[52, 152, 219];[133, 193, 233];[214, 234, 248]]./255;

    %Land colormap
    topocmapgreen=[[56, 142, 60]]./255; %Green: 1/30 land colormap
    topocmapyellow=[[252, 243, 207];[249, 231, 159];[247, 220, 111];[241, 196, 15]]./255; %Yellow: 3/30 land colormap
    topocmapbrown=[[243, 156, 18];[230, 126, 34];[211, 84, 0];[186, 74, 0];[160, 64, 0];[135, 54, 0];[110, 44, 0]]./255;
    topocmapwhite=[[236, 240, 241]]./255; %White: 1/30 land colormap

    %Interpolating brwon color from 7 colors to 24 colors
    topocmapbrown1(:,1:3)=interp1(linspace(1,24,length(topocmapbrown(:,1))),topocmapbrown(:,1:3),linspace(1,24,24)); %Brown: 24/30 land colormap

    %Assigning water and land colormap
    cmap_water1=topocmap_water;
    cmap_land1=[topocmapgreen;topocmapyellow;topocmapbrown1;topocmapwhite];

    %Convert to 0.55% sea and 0.45% land 
    %Mariana Trench depth= 10994 m
    %Mount Everest height= 8848 m
    %(max depth)/(max depth+max height)~=0.55
    cmap_water(:,1:3)=interp1(linspace(1,55,length(cmap_water1(:,1))),cmap_water1(:,1:3),linspace(1,55,55));
    cmap_land(:,1:3)=interp1(linspace(1,45,length(cmap_land1(:,1))),cmap_land1(:,1:3),linspace(1,45,45));

    %Assigning negative (water) and positive (land) colormap
    cmap_z=[cmap_water;cmap_land];

%'topocmaprelief' colormap
%Colors from globalrelief.m
elseif strcmp(cmapcolors,'topocmaprelief')==1

    %Water colormap
    topocmap_water=bone(128);
    topocmap_water=topocmap_water(65:128,:);

    %Land colormap
    topocmapsummer=summer(10);
    topocmapcopper=flipud(copper(54));

    %Assigning water and land colormap
    cmap_water1=topocmap_water;
    cmap_land1=[topocmapsummer;topocmapcopper];

    %Convert to 0.55% sea and 0.45% land
    %Mariana Trench depth= 10994 m
    %Mount Everest height= 8848 m
    %(max depth)/(max depth+max height)~=0.55
    cmap_water(:,1:3)=interp1(linspace(1,55,length(cmap_water1(:,1))),cmap_water1(:,1:3),linspace(1,55,55));
    cmap_land(:,1:3)=interp1(linspace(1,45,length(cmap_land1(:,1))),cmap_land1(:,1:3),linspace(1,45,45));

    %Assigning negative (water) and positive (land) colormap
    cmap_z=[cmap_water;cmap_land];

%'topocmapocean' colormap
%Colors from: http://htmlcolorcodes.com
elseif strcmp(cmapcolors,'topocmapocean')==1

    %Water colormap
    topocmap_water=[[21, 67, 96];[27, 79, 114];[40, 116, 166];[52, 152, 219];[133, 193, 233];[214, 234, 248]]./255;

    cmap_z=topocmap_water;

    %Assigning half of the colormap to negative (water) values and other half to positive (land) values 
    [M1,N1]=size(cmap_z);
    cmap_water=cmap_z(1:fix(M1/2),:);
    cmap_land=cmap_z(fix(M1/2)+1:end,:);

%'topocmapearth' colormap
%Colors from: http://htmlcolorcodes.com
elseif strcmp(cmapcolors,'topocmapearth')==1

    %Land colormap
    topocmapgreen=[[56, 142, 60]]./255; %Green: 1/30 land colormap
    topocmapyellow=[[252, 243, 207];[249, 231, 159];[247, 220, 111];[241, 196, 15]]./255; %Yellow: 3/30 land colormap
    topocmapbrown=[[243, 156, 18];[230, 126, 34];[211, 84, 0];[186, 74, 0];[160, 64, 0];[135, 54, 0];[110, 44, 0]]./255;
    topocmapwhite=[[236, 240, 241]]./255; %White: 1/30 land colormap

    %Interpolating brwon color from 7 colors to 24 colors
    topocmapbrown1(:,1:3)=interp1(linspace(1,24,length(topocmapbrown(:,1))),topocmapbrown(:,1:3),linspace(1,24,24)); %Brown: 24/30 land colormap

    %Assigning water and land colormap
    cmap_land1=[topocmapgreen;topocmapyellow;topocmapbrown1;topocmapwhite];

    cmap_z=cmap_land1;

    %Assigning half of the colormap to negative (water) values and other half to positive (land) values 
    [M1,N1]=size(cmap_z);
    cmap_water=cmap_z(1:fix(M1/2),:);
    cmap_land=cmap_z(fix(M1/2)+1:end,:);

%'blueocean' colormap
%Colors from: http://htmlcolorcodes.com
elseif strcmp(cmapcolors,'blueocean')==1

    %Water colormap
    blueoceanwater=[[52, 152, 219];[133, 193, 233]]./255;

    %Land colormap
    topocmapgreen=[[56, 142, 60]]./255; %Green: 1/30 land colormap
    topocmapyellow=[[252, 243, 207];[249, 231, 159];[247, 220, 111];[241, 196, 15]]./255; %Yellow: 3/30 land colormap
    topocmapbrown=[[243, 156, 18];[230, 126, 34];[211, 84, 0];[186, 74, 0];[160, 64, 0];[135, 54, 0];[110, 44, 0]]./255;
    topocmapwhite=[[236, 240, 241]]./255; %White: 1/30 land colormap

    %Interpolating brwon color from 7 colors to 24 colors
    topocmapbrown1(:,1:3)=interp1(linspace(1,24,length(topocmapbrown(:,1))),topocmapbrown(:,1:3),linspace(1,24,24)); %Brown: 24/30 land colormap

    %Assigning water and land colormap
    cmap_water1=blueoceanwater;
    cmap_land1=[topocmapgreen;topocmapyellow;topocmapbrown1;topocmapwhite];

    %Convert to 0.55% sea and 0.45% land 
    %Mariana Trench depth= 10994 m
    %Mount Everest height= 8848 m
    %(max depth)/(max depth+max height)~=0.55
    cmap_water(:,1:3)=interp1(linspace(1,55,length(cmap_water1(:,1))),cmap_water1(:,1:3),linspace(1,55,55));
    cmap_land(:,1:3)=interp1(linspace(1,45,length(cmap_land1(:,1))),cmap_land1(:,1:3),linspace(1,45,45));

    %Assigning negative (water) and positive (land) colormap
    cmap_z=[cmap_water;cmap_land];

%'blueoceansea' colormap
%Colors from: http://htmlcolorcodes.com
elseif strcmp(cmapcolors,'blueoceansea')==1

    %Water colormap
    blueoceanwater=[[52, 152, 219];[133, 193, 233]]./255;

    %At least 4 colors requred in colormap to guarantee at least two colors for each of cmap_water and cmap_land colormaps
    blueoceanwater4colors(:,1:3)=interp1(linspace(1,4,length(blueoceanwater(:,1))),blueoceanwater(:,1:3),linspace(1,4,4));

    cmap_z=blueoceanwater4colors;

    %Assigning half of the colormap to negative (water) values and other half to positive (land) values 
    [M1,N1]=size(cmap_z);
    cmap_water=cmap_z(1:fix(M1/2),:);
    cmap_land=cmap_z(fix(M1/2)+1:end,:);

%'greenearth' colormap
%Colors from: http://htmlcolorcodes.com
elseif strcmp(cmapcolors,'greenearth')==1

    %Mariana Trench depth= 10994 m
    %Mount Everest height= 8848 m
    %(max depth)/(max depth+max height)~=0.55

    %Water colormap
    topocmap_water=[[21, 67, 96];[27, 79, 114];[40, 116, 166];[52, 152, 219];[133, 193, 233];[214, 234, 248]]./255;

    %Land colormap
    greenearthgreen=[[67, 160, 71];[255, 241, 118]]./255; %Green

    %Assigning water and land colormap
    cmap_water1=topocmap_water;
    cmap_land1=greenearthgreen;

    %Convert to 0.55% sea and 0.45% land 
    %Mariana Trench depth= 10994 m
    %Mount Everest height= 8848 m
    %(max depth)/(max depth+max height)~=0.55
    cmap_water(:,1:3)=interp1(linspace(1,55,length(cmap_water1(:,1))),cmap_water1(:,1:3),linspace(1,55,55));
    cmap_land(:,1:3)=interp1(linspace(1,45,length(cmap_land1(:,1))),cmap_land1(:,1:3),linspace(1,45,45));

    %Assigning negative (water) and positive (land) colormap
    cmap_z=[cmap_water;cmap_land];

%'greenearthland' colormap
%Colors from: http://htmlcolorcodes.com
elseif strcmp(cmapcolors,'greenearthland')==1

    %Mariana Trench depth= 10994 m
    %Mount Everest height= 8848 m
    %(max depth)/(max depth+max height)~=0.55

    %Land colormap
    greenearthgreen=[[67, 160, 71];[255, 241, 118]]./255; %Green

    %At least 4 colors requred in colormap to guarantee at least two colors for each of cmap_water and cmap_land colormaps
    greenearthgreen4colors(:,1:3)=interp1(linspace(1,4,length(greenearthgreen(:,1))),greenearthgreen(:,1:3),linspace(1,4,4));

    cmap_z=greenearthgreen4colors;

    %Assigning half of the colormap to negative (water) values and other half to positive (land) values 
    [M1,N1]=size(cmap_z);
    cmap_water=cmap_z(1:fix(M1/2),:);
    cmap_land=cmap_z(fix(M1/2)+1:end,:);

%'blgrtopocmap' colormap
%Colors from: http://htmlcolorcodes.com
elseif strcmp(cmapcolors,'blgrtopocmap')==1

    %Mariana Trench depth= 10994 m
    %Mount Everest height= 8848 m
    %(max depth)/(max depth+max height)~=0.55

    %Water colormap
    blgrtopocmap_water=[[52, 152, 219];[133, 193, 233]]./255;

    %Land colormap
    blgrtopocmapgreen=[[67, 160, 71];[255, 241, 118]]./255; %Green

    %Assigning water and land colormap
    cmap_water1=blgrtopocmap_water;
    cmap_land1=blgrtopocmapgreen;

    %Convert to 0.55% sea and 0.45% land 
    %Mariana Trench depth= 10994 m
    %Mount Everest height= 8848 m
    %(max depth)/(max depth+max height)~=0.55
    cmap_water(:,1:3)=interp1(linspace(1,55,length(cmap_water1(:,1))),cmap_water1(:,1:3),linspace(1,55,55));
    cmap_land(:,1:3)=interp1(linspace(1,45,length(cmap_land1(:,1))),cmap_land1(:,1:3),linspace(1,45,45));

    %Assigning negative (water) and positive (land) colormap
    cmap_z=[cmap_water;cmap_land];

%'blrdtopocmap' colormap
%Colors from: http://htmlcolorcodes.com
elseif strcmp(cmapcolors,'blrdtopocmap')==1

    %Mariana Trench depth= 10994 m
    %Mount Everest height= 8848 m
    %(max depth)/(max depth+max height)~=0.55

    %Water colormap
    blrdtopocmap_water=[[41, 128, 185];[212, 230, 241]]./255; %Blue

    %Land colormap
    blrdtopocmap_land=[[242, 215, 213];[192, 57, 43]]./255; %Red

    %Assigning water and land colormap
    cmap_water1=blrdtopocmap_water;
    cmap_land1=blrdtopocmap_land;

    %Convert to 0.55% sea and 0.45% land 
    %Mariana Trench depth= 10994 m
    %Mount Everest height= 8848 m
    %(max depth)/(max depth+max height)~=0.55
    cmap_water(:,1:3)=interp1(linspace(1,55,length(cmap_water1(:,1))),cmap_water1(:,1:3),linspace(1,55,55));
    cmap_land(:,1:3)=interp1(linspace(1,45,length(cmap_land1(:,1))),cmap_land1(:,1:3),linspace(1,45,45));

    %Assigning negative (water) and positive (land) colormap
    cmap_z=[cmap_water;cmap_land];

%'grayearth' colormap
%Colors from: http://htmlcolorcodes.com
elseif strcmp(cmapcolors,'grayearth')==1

    %Mariana Trench depth= 10994 m
    %Mount Everest height= 8848 m
    %(max depth)/(max depth+max height)~=0.55

    %Water colormap
    topocmap_water=[[21, 67, 96];[27, 79, 114];[40, 116, 166];[52, 152, 219];[133, 193, 233];[214, 234, 248]]./255;

    %Land colormap
    grayearthgray=[[179, 182, 183];[66, 73, 73]]./255; %Gray

    %Assigning water and land colormap
    cmap_water1=topocmap_water;
    cmap_land1=grayearthgray;

    %Convert to 0.55% sea and 0.45% land 
    %Mariana Trench depth= 10994 m
    %Mount Everest height= 8848 m
    %(max depth)/(max depth+max height)~=0.55
    cmap_water(:,1:3)=interp1(linspace(1,55,length(cmap_water1(:,1))),cmap_water1(:,1:3),linspace(1,55,55));
    cmap_land(:,1:3)=interp1(linspace(1,45,length(cmap_land1(:,1))),cmap_land1(:,1:3),linspace(1,45,45));

    %Assigning negative (water) and positive (land) colormap
    cmap_z=[cmap_water;cmap_land];

%ETOPO1 colormap
%https://www.ngdc.noaa.gov/mgg/global/global.html
elseif strcmp(cmapcolors,'etopo1')==1

    %http://svn.idldev.com/vis/trunk/src/color/cpt-city/ngdc/ETOPO1.cpt
    % Colormap used in the ETOPO1 global relief map:
    % http://ngdc.noaa.gov/mgg/global/global.html
    %
    % Above sea level is a modified version of GMT_globe.cpt, 
    % Designed by Lester M. Anderson (CASP, UK) lester.anderson@casp.cam.ac.uk,
    % Modified by Jesse Varner and Elliot Lim (NOAA/NGDC) with a smaller band of white for the highest elevations.
    % The ocean colors are adapted from GMT_haxby.cpt, popularized by Bill Haxby, LDEO
    % COLOR_MODEL = RGB

    %ETOPO1 elavation range
    %min elevation= -10898.0
    %max elevation= 8271.0

    etopo1=[[10,0,121];[26,0,137];[38,0,152];[27,3,166];[16,6,180];[5,9,193];[0,14,203];[0,22,210];[0,30,216];[0,39,223];...
        [12,68,231];[26,102,240];[19,117,244];[14,133,249];[21,158,252];[30,178,255];[43,186,255];[55,193,255];[65,200,255];[79,210,255];...
        [94,223,255];[138,227,255];[188,230,255];[51,102,0];[51,204,102];[187,228,146];[255,220,185];[243,202,137];[230,184,88];[217,166,39];[168,154,31];...
        [164,144,25];[162,134,19];[159,123,13];[156,113,7];[153,102,0];[162,89,89];[178,118,118];[183,147,147];[194,176,176];[204,204,204];...
        [229,229,229];[255,255,255]]./255;

    etopo1water=[[10,0,121];[26,0,137];[38,0,152];[27,3,166];[16,6,180];[5,9,193];[0,14,203];[0,22,210];[0,30,216];[0,39,223];...
        [12,68,231];[26,102,240];[19,117,244];[14,133,249];[21,158,252];[30,178,255];[43,186,255];[55,193,255];[65,200,255];[79,210,255];...
        [94,223,255];[138,227,255];[188,230,255]]./255;

    etopo1land=[[51,102,0];[51,204,102];[187,228,146];[255,220,185];[243,202,137];[230,184,88];[217,166,39];[168,154,31];...
        [164,144,25];[162,134,19];[159,123,13];[156,113,7];[153,102,0];[162,89,89];[178,118,118];[183,147,147];[194,176,176];[204,204,204];...
        [229,229,229];[255,255,255]]./255;

    %Assigning negative (water) and positive (land) colormap
    cmap_water=etopo1water;
    cmap_land=etopo1land;
    cmap_z=etopo1;

%GMT_globe colormap
%https://www.giss.nasa.gov/tools/panoply/colorbars/
elseif strcmp(cmapcolors,'gmtglobe')==1

    %https://www.giss.nasa.gov/tools/panoply/colorbars/
    %https://www.giss.nasa.gov/tools/panoply/colorbars/gmt/GMT_globe.cpt
    % $Id: GMT_globe.cpt,v 1.1.1.1 2000/12/28 01:23:45 gmt Exp $
    %
    % Colormap using in global relief maps
    % Bathymetry colours manually redefined for blue-shade effect and
    % new topography colour scheme for use with Generic Mapping Tools.
    % Designed by Lester M. Anderson (CASP, UK) lester.anderson@casp.cam.ac.uk
    % COLOR_MODEL = RGB

    gmtglobe=[[153,0,255];[153,0,255];[153,0,255];[136,17,255];[119,34,255];[102,51,255];[85,68,255];[68,85,255];[51,102,255];[34,119,255];...
        [17,136,255];[0,153,255];[27,164,255];[54,175,255];[81,186,255];[108,197,255];[134,208,255];[161,219,255];[188,230,255];[215,241,255];...
        [241,252,255];[51,102,0];[51,204,102];[187,228,146];[255,220,185];[243,202,137];[230,184,88];[217,166,39];[168,154,31];[164,144,25];...
        [162,134,19];[159,123,13];[156,113,7];[153,102,0];[162,89,89];[178,118,118];[183,147,147];[194,176,176];[204,204,204];[229,229,229];...
        [242,242,242];[255,255,255]]./255;

    gmtglobewater=[[153,0,255];[153,0,255];[153,0,255];[136,17,255];[119,34,255];[102,51,255];[85,68,255];[68,85,255];[51,102,255];[34,119,255];...
        [17,136,255];[0,153,255];[27,164,255];[54,175,255];[81,186,255];[108,197,255];[134,208,255];[161,219,255];[188,230,255];[215,241,255];...
        [241,252,255]]./255;

    gmtglobeland=[[51,102,0];[51,204,102];[187,228,146];[255,220,185];[243,202,137];[230,184,88];[217,166,39];[168,154,31];[164,144,25];...
        [162,134,19];[159,123,13];[156,113,7];[153,102,0];[162,89,89];[178,118,118];[183,147,147];[194,176,176];[204,204,204];[229,229,229];...
        [242,242,242];[255,255,255]]./255;

    %Assigning negative (water) and positive (land) colormap
    cmap_water=gmtglobewater;
    cmap_land=gmtglobeland;
    cmap_z=gmtglobe;

%GMT_relief colormap
%https://www.giss.nasa.gov/tools/panoply/colorbars/
elseif strcmp(cmapcolors,'gmtrelief')==1

    %https://www.giss.nasa.gov/tools/panoply/colorbars/
    %https://www.giss.nasa.gov/tools/panoply/colorbars/gmt/GMT_relief.cpt
    %	$Id: GMT_relief.cpt,v 1.1.1.1 2000/12/28 01:23:45 gmt Exp $
    %
    % Colortable for whole earth relief used in Wessel topomaps
    % Designed by P. Wessel and F. Martinez, SOEST
    % COLOR_MODEL = RGB

    gmtrelief=[[0,0,0];[0,5,25];[0,10,50];[0,80,125];[0,150,200];[86,197,184];[172,245,168];[211,250,211];[250,255,255];[70,120,50];[120,100,50];[146,126,60];...
        [198,178,80];[250,230,100];[250,234,126];[252,238,152];[252,243,177];[253,249,216];[255,255,255]]./255;

    gmtreliefwater=[[0,0,0];[0,5,25];[0,10,50];[0,80,125];[0,150,200];[86,197,184];[172,245,168];[211,250,211];[250,255,255]]./255;

    gmtreliefland=[[70,120,50];[120,100,50];[146,126,60];[198,178,80];[250,230,100];[250,234,126];[252,238,152];[252,243,177];[253,249,216];[255,255,255]]./255;

    %Assigning negative (water) and positive (land) colormap
    cmap_water=gmtreliefwater;
    cmap_land=gmtreliefland;
    cmap_z=gmtrelief;

%Colormap from  Florian Aendekerk
%http://www.mathworks.com/matlabcentral/fileexchange/63590-landseacolormap-m-
elseif strcmp(cmapcolors,'aendekerk')==1

    aendekerkcmap=[[0.00000, 0.00000, 0.40000];[0.02524, 0.03810, 0.42857];[0.05048, 0.07619, 0.45714];[0.07571, 0.11429, 0.48571];[0.10095, 0.15238, 0.51429];...
        [0.12619, 0.19048, 0.54286];[0.15143, 0.22857, 0.57143];[0.17667, 0.26667, 0.60000];[0.20190, 0.30476, 0.62857];[0.22714, 0.34286, 0.65714];...
        [0.25238, 0.38095, 0.68571];[0.27762, 0.41905, 0.71429];[0.30286, 0.45714, 0.74286];[0.32810, 0.49524, 0.77143];[0.35333, 0.53333, 0.80000];...
        [0.37857, 0.57143, 0.82857];[0.40381, 0.60952, 0.85714];[0.42905, 0.64762, 0.88571];[0.45429, 0.68571, 0.91429];[0.47952, 0.72381, 0.94286];...
        [0.50476, 0.76190, 0.97143];[0.53000, 0.80000, 1.00000];[0.00000, 0.50000, 0.00000];[1.00000, 1.00000, 0.80000];[1.00000, 0.86667, 0.53333];...
        [1.00000, 0.73333, 0.26667];[1.00000, 0.60000, 0.00000];[1.00000, 0.60000, 0.00000];[0.95000, 0.56250, 0.00000];[0.90000, 0.52500, 0.00000];...
        [0.85000, 0.48750, 0.00000];[0.80000, 0.45000, 0.00000];[0.75000, 0.41250, 0.00000];[0.70000, 0.37500, 0.00000];[0.65000, 0.33750, 0.00000];...
        [0.60000, 0.30000, 0.00000];[0.55000, 0.26250, 0.00000];[0.50000, 0.22500, 0.00000];[0.45000, 0.18750, 0.00000];[0.40000, 0.15000, 0.00000]];

    aendekerkcmap_water=[[0.00000, 0.00000, 0.40000];[0.02524, 0.03810, 0.42857];[0.05048, 0.07619, 0.45714];[0.07571, 0.11429, 0.48571];[0.10095, 0.15238, 0.51429];...
        [0.12619, 0.19048, 0.54286];[0.15143, 0.22857, 0.57143];[0.17667, 0.26667, 0.60000];[0.20190, 0.30476, 0.62857];[0.22714, 0.34286, 0.65714];...
        [0.25238, 0.38095, 0.68571];[0.27762, 0.41905, 0.71429];[0.30286, 0.45714, 0.74286];[0.32810, 0.49524, 0.77143];[0.35333, 0.53333, 0.80000];...
        [0.37857, 0.57143, 0.82857];[0.40381, 0.60952, 0.85714];[0.42905, 0.64762, 0.88571];[0.45429, 0.68571, 0.91429];[0.47952, 0.72381, 0.94286];...
        [0.50476, 0.76190, 0.97143];[0.53000, 0.80000, 1.00000]];

    aendekerkcmap_land=[[0.00000, 0.50000, 0.00000];[1.00000, 1.00000, 0.80000];[1.00000, 0.86667, 0.53333];...
        [1.00000, 0.73333, 0.26667];[1.00000, 0.60000, 0.00000];[1.00000, 0.60000, 0.00000];[0.95000, 0.56250, 0.00000];[0.90000, 0.52500, 0.00000];...
        [0.85000, 0.48750, 0.00000];[0.80000, 0.45000, 0.00000];[0.75000, 0.41250, 0.00000];[0.70000, 0.37500, 0.00000];[0.65000, 0.33750, 0.00000];...
        [0.60000, 0.30000, 0.00000];[0.55000, 0.26250, 0.00000];[0.50000, 0.22500, 0.00000];[0.45000, 0.18750, 0.00000];[0.40000, 0.15000, 0.00000]];

    %Assigning negative (water) and positive (land) colormap
    cmap_water=aendekerkcmap_water;
    cmap_land=aendekerkcmap_land;
    cmap_z=aendekerkcmap;

%User input colormap
elseif isnumeric(cmapcolors)==1

    cmap_z=cmapcolors./255; %Converting RGB to values between 0 and 1
    cmap_z(cmap_z<0 )=0;
    cmap_z(cmap_z>1 )=1;

    %Assigning half of the colormap to negative (water) values and other half to positive (land) values 
    [M1,N1]=size(cmap_z);
    cmap_water=cmap_z(1:fix(M1/2),:);
    cmap_land=cmap_z(fix(M1/2)+1:end,:);

%Pre-defined colormap such as 'cool', 'winter', etc
else

    cmap_z=colormap(cmapcolors);

    %Assigning half of the colormap to negative (water) values and other half to positive (land) values 
    [M1,N1]=size(cmap_z);
    cmap_water=cmap_z(1:fix(M1/2),:);
    cmap_land=cmap_z(fix(M1/2)+1:end,:);

end

%--------------------------------------------------------------------------
%Interpolating a colormap to include ncolor colors

%Define contour levels
if strcmp(waterlandcmapratio,'auto')==1

    zmingrid=nanmin(zgrid(:));
    zmaxgrid=nanmax(zgrid(:));

    if zmingrid<0 & zmaxgrid>0
        num_z_levels=ncolor+1;
        z_range=max(zgrid(:))-min(zgrid(:));
        delta_z=(z_range/(num_z_levels));
        z_levels_neg=[-delta_z:-delta_z:min(zgrid(:))];
        z_levels_pos=[0:delta_z:max(zgrid(:))];
        z_levels=[fliplr(z_levels_neg),z_levels_pos];
    else
        z_levels=ncolor+1;
    end

end

%Interpolate colormap
if strcmp(waterlandcmapratio,'auto')==1

    %z data have both negative and positive values
    if zmingrid<0 & zmaxgrid>0

        %Defining negative (water) colormap
        ncolor_neg=length(z_levels_neg);
        cmap_ncolor_neg=zeros(ncolor_neg,3);
        cmap_ncolor_neg(:,1:3)=interp1(linspace(1,ncolor_neg,length(cmap_water(:,1))),cmap_water(:,1:3),linspace(1,ncolor_neg,ncolor_neg));

        %Defining positive (land) colormap
        ncolor_pos=length(z_levels_pos)-1;
        cmap_ncolor_pos=zeros(ncolor_pos,3);
        cmap_ncolor_pos(:,1:3)=interp1(linspace(1,ncolor_pos,length(cmap_land(:,1))),cmap_land(:,1:3),linspace(1,ncolor_pos,ncolor_pos));

        %Defining colormap that include ncolor colors
        cmap_ncolor=[cmap_ncolor_neg;cmap_ncolor_pos];

        %Scale z data
        zgrid(zgrid>max(z_levels(:)))=max(z_levels(:));
        zgrid(zgrid<min(z_levels(:)))=min(z_levels(:));

    %z data onle have positive values
    elseif zmingrid>=0 & zmaxgrid>=0
    
        ncolor_neg=0;
        ncolor_pos=ncolor;

        %Defining land colormap
        cmap_ncolor_pos=zeros(ncolor_pos,3);
        cmap_ncolor_pos(:,1:3)=interp1(linspace(1,ncolor_pos,length(cmap_land(:,1))),cmap_land(:,1:3),linspace(1,ncolor_pos,ncolor_pos));

        %Defining colormap that include ncolor colors
        cmap_ncolor=cmap_ncolor_pos;

    %z data only have negative values
    elseif zmingrid<=0 & zmaxgrid<=0

        ncolor_neg=ncolor;
        ncolor_pos=0;

        %Defining water colormap
        cmap_ncolor_neg=zeros(ncolor_neg,3);
        cmap_ncolor_neg(:,1:3)=interp1(linspace(1,ncolor_neg,length(cmap_water(:,1))),cmap_water(:,1:3),linspace(1,ncolor_neg,ncolor_neg));

        %Defining colormap that include ncolor colors
        cmap_ncolor=cmap_ncolor_neg;

    end

%Colormap interpolated to include ncolor colors without considering water_land ratio
elseif strcmp(waterlandcmapratio,'none')==1

    cmap_ncolor=zeros(ncolor,3);
    cmap_ncolor(:,1:3)=interp1(linspace(1,ncolor,length(cmap_z(:,1))),cmap_z(:,1:3),linspace(1,ncolor,ncolor));

%Assigning waterlandcmapratio of the colormap to negative (water) values and the rest to positive (land) values 
elseif isnumeric(waterlandcmapratio)==1

    %Checking waterlandcmapratio value
    waterlandcmapratio(waterlandcmapratio<=0 | waterlandcmapratio>1)=0.5;

    %Assigning water_land ratio to to negative (water) and positive (land) colormaps
    %[M1,N1]=size(cmap_z);
    %cmap_waterscaled=cmap_z(1:fix(M1/2),:);
    %cmap_landscaled=cmap_z(fix(M1/2)+1:end,:);

    %Defining water colormap
    ncolor_neg=(fix(ncolor*waterlandcmapratio));
    cmap_ncolor_neg=zeros(ncolor_neg,3);
    cmap_ncolor_neg(:,1:3)=interp1(linspace(1,ncolor_neg,length(cmap_water(:,1))),cmap_water(:,1:3),linspace(1,ncolor_neg,ncolor_neg));

    %Defining land colormap
    ncolor_pos=ncolor-ncolor_neg;
    %ncolor_pos=ncolor-(fix(M1/2)*waterlandcmapratio);
    cmap_ncolor_pos=zeros(ncolor_pos,3);
    cmap_ncolor_pos(:,1:3)=interp1(linspace(1,ncolor_pos,length(cmap_land(:,1))),cmap_land(:,1:3),linspace(1,ncolor_pos,ncolor_pos));

    %Defining colormap that include ncolor colors
    cmap_ncolor=[cmap_ncolor_neg;cmap_ncolor_pos];

end

%--------------------------------------------------------------------------
%Displaying results

if strcmp(plottype,'imagesc')==1 %2 dimensional plot using imagesc or imshow
    
    imagesc([min(xgrid(:)),max(xgrid(:))],[min(ygrid(:)),max(ygrid(:))],zgrid)
    set(gca,'YDir','normal'); %Correct y axis direction

elseif strcmp(plottype,'pcolor')==1 %2 dimensional plot using pcolor

    pcol=pcolor(xgrid,ygrid,zgrid);
    set(pcol,'edgecolor','none')

elseif strcmp(plottype,'contour')==1 %2 dimensional contour plot

    zmingrid=nanmin(zgrid(:));
    zmaxgrid=nanmax(zgrid(:));
    [Cont,hCont]=contourf(xgrid,ygrid,zgrid,z_levels);
    set(hCont,'LineColor','none');

elseif strcmp(plottype,'surface')==1 %Surface plot
    
    surf1=surf(xgrid,ygrid,zgrid);
    %set(surf1,'xscale','log','zscale','log')
    set(surf1,'FaceColor','interp','EdgeColor','none')

end

%Setting plot properties
if (strcmp(plottype,'imagesc')==1 | strcmp(plottype,'pcolor')==1 | strcmp(plottype,'contour')==1 | strcmp(plottype,'surface')==1)

    %Setting axis limits
    xlim([nanmin(xgrid(:)),nanmax(xgrid(:))])
    ylim([nanmin(ygrid(:)),nanmax(ygrid(:))])

    %Setting colormap
    colormap(cmap_ncolor)

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
