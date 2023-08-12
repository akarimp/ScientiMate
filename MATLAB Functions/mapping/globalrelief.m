function [x, y, z, xyz, xgrid, ygrid, zgrid] = globalrelief(filepath, xmin, xmax, ymin, ymax, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2020-09-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

globalrelief
============

.. code:: MATLAB

    [x, y, z, xyz, xgrid, ygrid, zgrid] = globalrelief(filepath, xmin, xmax, ymin, ymax, dispout)

Description
-----------

| Return x (longitude), y (latitude) and z (elevation) data from ETOPO1 Global Relief Model (Amante & Eakins, 2009) interpolated on 0.08 degree grid
| ETOPO1 is 1 arc-minute global relief, however, this data are interpolated on 4.8 arc-minute (0.08 degree)
| This data are obtained from ETOPO1 Global Relief Bedrock (grid-registered)
| ETOPO1 horizontal datum: WGS 84 geographic
| ETOPO1 vertical datum: sea level
| https://www.ngdc.noaa.gov/mgg/global/global.html

Inputs
------

filepath
    | Path of the folder that contains 'ETOPO1_XYZ_0_08_Deg_Matlab.mat' file
    | Download ETOPO1_XYZ_0_08_Deg_Matlab.mat file from
    | https://github.com/akarimp/ScientiMate/releases/download/2.0/ETOPO1_XYZ_0_08_Deg_Matlab.mat
    | Example: filepath = 'C:'
xmin=-180
    | Minimum x (longitude) of domain to be returned in degree 
    | It should be between -180 and 180 degree
xmax=180
    | Maximum x (longitude) of domain to be returned in degree
    | It should be between -180 and 180 degree
ymin=-90
    | Minimum y (latitude) of domain to be returned in degree
    | It should be between -90 and 90 degree
ymax=90
    | Maximum y (latitude) of domain to be returned in degree
    | It should be between -90 and 90 degree
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

x
    Interpolated x (longitude) data in degree
y
    Interpolated y (latitude) data in degree
z
    Interpolated z (elevation) data in degree
xyz
    | Interpolated xyz data
    | xyz is a 3-column array
    | 1st column contains longitude (x) data in (Degree)
    | 2nd column contains latitude (y) data in (Degree)
    | 3rd column contains elevation (z) data in (m)
xgrid
    Interpolated x (longitude) data on 2d mesh in degree
ygrid
    Interpolated y (latitude) data on 2d mesh in degree
zgrid
    Interpolated z (elevation) data on 2d mesh in degree

Examples
--------

.. code:: MATLAB

    %Path of ETOPO1_XYZ_0_08_Deg_Matlab.mat
    filepath = 'C:'

    %Globe
    [x,y,z,xyz,xgrid,ygrid,zgrid]=globalrelief(filepath,-180,180,-90,90,'yes');

    %Middle East
    [x,y,z,xyz,xgrid,ygrid,zgrid]=globalrelief(filepath,24,64,9,43,'yes');

    %North America
    [x,y,z,xyz,xgrid,ygrid,zgrid]=globalrelief(filepath,-169,-8,5,90,'yes');


References
----------

ETOPO1 Global Relief Model

    Amante, C. and B.W. Eakins, 2009.
    ETOPO1 1 Arc-Minute Global Relief Model: Procedures, Data Sources and Analysis.
    NOAA Technical Memorandum NESDIS NGDC-24. National Geophysical Data Center, NOAA.
    doi:10.7289/V5C8276M

* https://www.ngdc.noaa.gov/mgg/global/global.html
* https://maps.ngdc.noaa.gov/viewers/grid-extract/index.html
* https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/
* https://data.nodc.noaa.gov/cgi-bin/iso?id=gov.noaa.ngdc.mgg.dem:316

GEBCO Global ocean & land terrain models

* https://www.gebco.net/data_and_products/gridded_bathymetry_data/

Natural Earth 1:10m Raster Data

* https://www.naturalearthdata.com/downloads/10m-raster-data/

Geospatial data

* https://www.mathworks.com/help/map/finding-geospatial-data.html
* https://www.ngdc.noaa.gov/mgg/global/etopo2.html
* https://www.ngdc.noaa.gov/mgg/global/etopo5.HTML
* https://www.ngdc.noaa.gov/mgg/image/2minrelief.html
* https://www.ngdc.noaa.gov/mgg/coastal/crm.html
* https://www.ngdc.noaa.gov/mgg
* https://viewer.nationalmap.gov/launch
* https://earthexplorer.usgs.gov
* http://www.shadedrelief.com/cleantopo2/index.html
* https://www.usna.edu/Users/oceano/pguth/md_help/html/bathymetry.htm
* https://www.opentopodata.org
* https://en.wikipedia.org/wiki/Global_relief_model

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
case 1
    xmin=-180; xmax=180; ymin=-90; ymax=90; dispout='no';
case 2
    xmax=180; ymin=-90; ymax=90; dispout='no';
case 3
    ymin=-90; ymax=90; dispout='no';
case 4
    ymax=90; dispout='no';
case 5
    dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

%if isrow(x)==1
%x=x';
%end

%--------------------------------------------------------------------------
%Check input values

xmin(xmin<-180)=-180;
xmax(xmax>180)=180;
ymin(ymin<-90)=-90;
ymax(ymax>90)=90;

%--------------------------------------------------------------------------
%Read data

%Change folder to data folder
Current_Folder=pwd; %Current path
cd(filepath)

%The xyz data
xyz = importdata('ETOPO1_XYZ_0_08_Deg_Matlab.mat');
xyz = double(xyz);
x=xyz(:,1); %x data
y=xyz(:,2); %y data
z=xyz(:,3); %z data

%The z data
%z = importdata('ETOPO1_Z_0_125_Deg_Matlab.mat');

%Return to current folder
cd(Current_Folder)

%--------------------------------------------------------------------------
%xyz data

%The xy 2d and xy data
dx=0.08;
dy=0.08;
xm=[-180:dx:180];
yn=[-90:dy:90];
[Xnm_temp,Ynm_temp] = meshgrid(xm,yn);
%x=reshape(Xnm, [], 1); %Or x=Xnm(:);
%y=reshape(Ynm, [], 1); %Or y=Ynm(:);

%z 2d data
[M,N] = size(Xnm_temp);
Xnm = reshape(x, M,N);
Ynm = reshape(y, M,N);
Znm = reshape(z, M,N);
%Znm=(reshape(z, N, M)).'; %To get reshape similar to Numpy

%The xyz data
%xyz = [x,y,z];

%--------------------------------------------------------------------------
%Defining domain area and associated data for 1d data

%Retain data within a defined domain

%Finding a location data that are within [xmin:xmax] range
xIndx=find(x>=xmin & x<=xmax);

%Retaining data that are within [xmin:xmax] range
x1=x(xIndx,1); %x data
y1=y(xIndx,1); %y data
z1=z(xIndx,1); %z data

%Finding a location data that are within [ymin:ymax] range
yIndx=find(y1>=ymin & y1<=ymax);

%Retaining data that are within [ymin:ymax] range
x2=x1(yIndx,1); %x data
y2=y1(yIndx,1); %y data
z2=z1(yIndx,1); %z data

%--------------------------------------------------------------------------
%Defining domain area and associated data for 2d data

%Retain data within a 2d defined domain

%Finding a location data that are within [xmin:xmax] range (bounding box)
[xIndx_2d_row,xIndx_2d_col]=find(Xnm>=xmin & Xnm<=xmax);
i_start = min(xIndx_2d_row);
j_start = min(xIndx_2d_col);
i_end = max(xIndx_2d_row);
j_end = max(xIndx_2d_col);

%Retaining data that are within [xmin:xmax] range
Xnm_1 = Xnm(i_start:i_end, j_start:j_end); %x grid data
Ynm_1 = Ynm(i_start:i_end, j_start:j_end); %x grid data
Znm_1 = Znm(i_start:i_end, j_start:j_end); %x grid data

%Finding a location data that are within [ymin:ymax] range (bounding box)
[yIndx_2d_row,yIndx_2d_col]=find(Ynm_1>=ymin & Ynm_1<=ymax);
i_start = min(yIndx_2d_row);
j_start = min(yIndx_2d_col);
i_end = max(yIndx_2d_row);
j_end = max(yIndx_2d_col);

%Retaining data that are within [ymin:ymax] range
Xnm_2 = Xnm_1(i_start:i_end, j_start:j_end); %x grid data
Ynm_2 = Ynm_1(i_start:i_end, j_start:j_end); %x grid data
Znm_2 = Znm_1(i_start:i_end, j_start:j_end); %x grid data

%--------------------------------------------------------------------------
%Assigning new x, y and z

x=x2; %x data
y=y2; %y data
z=z2; %z data
xyz=[x,y,z]; %xyz data

xgrid=Xnm_2; %x data
ygrid=Ynm_2; %y data
zgrid=Znm_2; %z data

%--------------------------------------------------------------------------
%Colormap, developed by Arash Karimpour

cmap_water=bone(128);
cmap_water=cmap_water(65:128,:);
cmap_land_1=summer(10);
cmap_land_2=flipud(copper(54));
%cmap_z=[cmap_water;cmap_land_1;cmap_land_2];

%Assigning negative (water) and positive (land) colormap
if (min(z)<0 & max(z)>0)
    cmap_z=[cmap_water;cmap_land_1;cmap_land_2];
    z_lim=max([abs(max(z)),abs(min(z))]);
    z_lim_min=-z_lim;
    z_lim_max=z_lim;
elseif (min(z)>=0 & max(z)>0)
    %cmap_z=[cmap_land_1;cmap_land_2];
    cmap_z=cmap_land_2;
    z_lim_min=min(z);
    z_lim_max=max(z);
elseif (min(z)<0 & max(z)<=0)
    cmap_z=cmap_water;
    z_lim_min=min(z);
    z_lim_max=max(z);
end

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1
    
    %[Cont,hCont]=contourf(xgrid,ygrid,zgrid);
    %set(hCont,'LineColor','none');

    %pcol=pcolor(xgrid,ygrid,zgrid);
    %set(pcol,'edgecolor','none')

    imagesc([min(xgrid(:)),max(xgrid(:))],[min(ygrid(:)),max(ygrid(:))],zgrid)
    set(gca,'YDir','normal'); %Correct y axis direction

    xlabel('Longitude (degree)')
    ylabel('Latitude (degree)')
    colormap(cmap_z)
    colorbar
    caxis([z_lim_min, z_lim_max])
    
end
    
%--------------------------------------------------------------------------
