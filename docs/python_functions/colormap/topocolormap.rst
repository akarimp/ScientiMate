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

scientimate.topocolormap
========================

.. code:: python

    cmap_ncolor, cmap_ncolorplt, cmap_water, cmap_land = scientimate.topocolormap(zcolormap='topocmap', ncolor=256, dispout='no')

Description
-----------

Export a topographic colormap

Inputs
------

zcolormap='topocmap'
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
ncolor=256
    Number of colors to be used in colormap
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

cmap_ncolor
    | Colormap for z levels with ncolor number of colors in RGB color format between 0 and 1
    | To convert 0-1 scale to 0-255 scale, multiply cmap_ncolor values by 255
cmap_ncolorplt
    Colormap for z levels with ncolor number of colors in Matplotlib format
cmap_water
    Colormap for water in RGB color format between 0 and 1
cmap_land
    Colormap for land in RGB color format between 0 and 1

Examples
--------

.. code:: python

    import scientimate as sm
    cmap_ncolor,cmap_ncolorplt,cmap_water,cmap_land=sm.topocolormap('topocmap',256,'yes')
    cmap_ncolor,cmap_ncolorplt,cmap_water,cmap_land=sm.topocolormap('cool',256,'yes')

References
----------

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
