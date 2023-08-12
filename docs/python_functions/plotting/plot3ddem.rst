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

scientimate.plot3ddem
=====================

.. code:: python

    scientimate.plot3ddem(x, y, z, plottype='imagesc', cmapcolors='topocmap')

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
plottype='imagesc'
    | Plot type
    | 'imagesc': 2 dimensional plot using imagesc or imshow
    | 'pcolor': 2 dimensional plot using pcolor
    | 'contour': 2 dimensional contour plot, number of contour=ncolor
    | 'surface': 3 dimensional surface plot
cmapcolors='topocmap'
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

.. code:: python

    import scientimate as sm
    import numpy as np

    x,y=np.meshgrid(np.linspace(-10,10,50),np.linspace(-10,10,50))
    r=np.sqrt(x**2+y**2)+1e-10 #Add 1e-10 to prevent divide by 0
    z=np.sin(r)/r
    sm.plot3ddem(x,y,z,'pcolor','topocmap')

    x,y=np.meshgrid(np.linspace(-10,10,21),np.linspace(-10,10,21))
    z=(np.sin(x)+np.sin(y))/(x+y+1e-10) #Add 1e-10 to prevent divide by 0
    sm.plot3ddem(x,y,z,'pcolor','topocmap')

    x=10*np.random.rand(1000)
    y=10*np.random.rand(1000)
    z=x**2+y**2
    sm.plot3ddem(x,y,z,'pcolor','topocmap')

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
