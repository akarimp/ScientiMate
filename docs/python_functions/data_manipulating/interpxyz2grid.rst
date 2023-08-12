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

scientimate.interpxyz2grid
==========================

.. code:: python

    xgrid, ygrid, zgrid = scientimate.interpxyz2grid(x, y, z, gridsize=100, gridsizetype='points', xmin=None, xmax=None, ymin=None, ymax=None, zmin=None, zmax=None, RetainRatio='all', interpMethod='nearest', dispout='no')

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
gridsize=100
    | Grid size in x (longitude) and y (latitude) directions to interpolate elevation data on them
    | if gridsizetype='length' then gridsize is a distance between grid points
    | if gridsizetype='points' then gridsize is number of grid points in each direction
gridsizetype='points'
    | Grid size type 
    | 'points': gridsize is considered as number of grid points in each direction
    | 'length': gridsize is considered as length between grid points
xmin=nanmin(x)
    Minimum x (longitude) of domain to be interpolated
xmax=nanmax(x)
    Maximum x (longitude) of domain to be interpolated
ymin=nanmin(y)
    Minimum y (latitude) of domain to be interpolated
ymax=nanmax(y)
    Maximum y (latitude) of domain to be interpolated
zmin=nanmin(z)
    | Minimum z (elevation) of domain to be interpolated
    | All z<zmin would be set to zmin
zmax=nanmax(z)
    | Maximum z (elevation) of domain to be interpolated
    | All z>zmax would be set to zmax
RetainRatio='all'
    | Define to down sample input data or not 
    | 'all': data are not down sampled
    | value between 0 and 1: percentage of retaining data
    | RetainRatio=0.8 : 80% of data are retained
interpMethod='nearest'
    | Interpolation method
    | 'linear': Use default or 'linear' method to interpolate
    | 'nearest': Use nearest neighbor method to interpolate
dispout='no'
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

.. code:: python

    import scientimate as sm
    import numpy as np

    x=10.*np.random.rand(1000)
    y=10.*np.random.rand(1000)
    z=x**2+y**2
    xgrid,ygrid,zgrid=sm.interpxyz2grid(x,y,z,100,'points',np.nanmin(x),np.nanmax(x),np.nanmin(y),np.nanmax(y),np.nanmin(z),np.nanmax(z),'all','nearest','yes')

    x=(-90-(-91))*np.random.rand(1000)+(-91)
    y=(31-(30))*np.random.rand(1000)+(30)
    z=x**2+y**2
    xgrid,ygrid,zgrid=sm.interpxyz2grid(x,y,z,0.005,'length',np.nanmin(x),np.nanmax(x),np.nanmin(y),np.nanmax(y),np.nanmin(z),np.nanmax(z),'all','linear','yes')

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
