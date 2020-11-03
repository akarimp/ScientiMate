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

scientimate.globalrelief
========================

.. code:: python

    x, y, z, xgrid, ygrid, zgrid = scientimate.globalrelief(xmin=-180, xmax=180, ymin=-90, ymax=90, dispout='no')

Description
-----------

| Return x (longitude), y (latitude) and z (elevation) data from ETOPO1 Global Relief Model (Amante & Eakins, 2009) interpolated on 0.125 degree grid
| ETOPO1 is 1 arc-minute global relief, however, this data are interpolated on 7.5 arc-minute (0.125 degree)
| This data are obtained from ETOPO1 Global Relief Bedrock (grid-registered)
| ETOPO1 horizontal datum: WGS 84 geographic
| ETOPO1 vertical datum: sea level
| https://www.ngdc.noaa.gov/mgg/global/global.html

Inputs
------

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
xgrid
    Interpolated x (longitude) data on 2d mesh in degree
ygrid
    Interpolated y (latitude) data on 2d mesh in degree
zgrid
    Interpolated z (elevation) data on 2d mesh in degree

Examples
--------

.. code:: python

    import scientimate as sm

    #Globe
    x,y,z,xgrid,ygrid,zgrid=sm.globalrelief(-180,-180,-90,90,'yes')

    #Middle East
    x,y,z,xgrid,ygrid,zgrid=sm.globalrelief(24,64,9,43,'yes')

    #North America
    x,y,z,xgrid,ygrid,zgrid=sm.globalrelief(-169,-8,5,90,'yes')

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
* https://viewer.nationalmap.gov/launch/
* https://earthexplorer.usgs.gov
* http://www.shadedrelief.com/cleantopo2/index.html
* https://www.usna.edu/Users/oceano/pguth/md_help/html/bathymetry.htm
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
