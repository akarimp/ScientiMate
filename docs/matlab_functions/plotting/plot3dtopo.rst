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

plot3dtopo
==========

.. code:: MATLAB

    [xgrid, ygrid, zgrid, cmap_ncolor] = plot3dtopo(x, y, z, ...
        zcolormap, ncolor, waterlandcmapratio, ...
        plotfontsize, axisfontsize, xaxislabel, yaxislabel, cbarlabel, dispcolorbar, dispgrid, ...
        gridsize, gridsizetype, xmin, xmax, ymin, ymax, zmin, zmax, RetainRatio, interpMethod, dispout)

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
zcolormap='topocmap';
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
ncolor=256;
    Number of colors to be used in colormap
waterlandcmapratio='none';
    | Scale of sea_colormap over land_colormap (between 0 to 1)
    | 'none': no scaling, original color map is used
    | 'auto': put z=0 at center of colormap
    | if equal to a number between 0 to 1:
    | assigning ratio equal to waterlandcmapratio of the colormap to negative (water) values
    | and the rest to positive (land) values
    | waterlandcmapratio=0.55; : assign 55% of colormap to negative (water) values
plotfontsize=12;
    Size of plot fonts
axisfontsize=12;
    Size of axis fonts
xaxislabel='x';
    x axis label
yaxislabel='y';
    y axis label
cbarlabel='z';
    Colorbar label
dispcolorbar='no';
    Define to display colorbar or not ('yes': display, 'no': not display)
dispgrid='no';
    | Define to display grid lines or not ('yes': display, 'no': not display)
    | 'yes': divide x axis and y axis to 13 intervals
    | if equal to a number:
    | divide x axis and y axis to dispgrid intervals
    | dispgrid=10; : divide x axis and y axis to 10 intervals
gridsize=100;
    | Grid size in x (longitude) and y (latitude) directions to interpolate elevation data on them
    | if gridsizetype='length' then gridsize is a distance between grid points
    | if gridsizetype='points' then gridsize is number of grid points in each direction
gridsizetype='points';
    | Grid size type
    | 'number': gridsize is considered as number of grid points in each direction
    | 'length': gridsize is considered as length between grid points
xmin=nanmin(x);
    Minimum x (longitude) of domain to be plotted
xmax=nanmax(x);
    Maximum x (longitude) of domain to be plotted
ymin=nanmin(y);
    Minimum y (latitude) of domain to be plotted
ymax=nanmax(y);
    Maximum y (latitude) of domain to be plotted
zmin=nanmin(z);
    | Minimum z (elevation) of domain to be plotted
    | All z<zmin would be set to zmin
zmax=nanmax(z);
    | Maximum z (elevation) of domain to be plotted
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
dispout='imagesc';
    | Define to display outputs or not
    | 'imagesc': 2 dimensional plot using imagesc or imshow
    | 'pcolor': 2 dimensional plot using pcolor
    | 'contour': 2 dimensional contour plot, number of contour=ncolor
    | 'surface': 3 dimensional surface plot
    | 'no': not display

Outputs
-------

xgrid
    Interpolated x (longitude) data on defined mesh
ygrid
    Interpolated y (latitude) data on defined mesh
zgrid
    Interpolated z (elevation) data on defined mesh
cmap_ncolor
    | Colormap for z levels with ncolor number of colors in RGB color format between 0 and 1
    | To convert 0-1 scale tp 0-255 scale, multiply cmap_ncolor values by 255

Examples
--------

.. code:: MATLAB

    x(:,1)=10.*rand(1000,1);
    y(:,1)=10.*rand(1000,1);
    z=x.^2+y.^2;
    [xgrid,ygrid,zgrid,cmap_ncolor]=plot3dtopo(x,y,z,...
        'topocmap',256,'none',...
        12,12,'x','y','z','no','no',...
        100,'points',nanmin(x(:)),nanmax(x(:)),nanmin(y(:)),nanmax(y(:)),nanmin(z(:)),nanmax(z(:)),'all','nearest','imagesc');

    x(:,1)=10.*rand(1000,1);
    y(:,1)=10.*rand(1000,1);
    z=x.^2+y.^2;
    [xgrid,ygrid]=meshgrid(linspace(nanmin(x),nanmax(x),100),linspace(nanmin(y),nanmax(y),100));
    zgrid=griddata(x,y,z,xgrid,ygrid);
    [xgrid,ygrid,zgrid,cmap_ncolor]=plot3dtopo(xgrid,ygrid,zgrid,...
        'topocmap',256,'none',...
        12,12,'x','y','z','no','no',...
        100,'points',nanmin(x(:)),nanmax(x(:)),nanmin(y(:)),nanmax(y(:)),nanmin(z(:)),nanmax(z(:)),'all','nearest','imagesc');

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
