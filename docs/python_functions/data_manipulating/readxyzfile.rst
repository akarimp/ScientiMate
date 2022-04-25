.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-11-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

scientimate.readxyzfile
=======================

.. code:: python

    x, y, z = scientimate.readxyzfile(xyzfilename, xyzfilelocation=None, zscale=1, domain='all', xmin=-180, xmax=180, ymin=-90, ymax=90, savedata='no', outfilename='xyzdata.csv', outfilelocation=None)

Description
-----------

| Read and extract x (longitude), y (latitude) and z (elevation) data from ASCII gridded (tabular) xyz file
| Use readdatafile function for more options

Inputs
------

xyzfilename
    | Name of xyz file between ' ' mark, example: 'xyzfile.xyz'
    | xyz file should be in form of 3 coloumn format
xyzfilelocation=pwd
    Location of xyz file between ' ' mark, example: 'C:\'
zscale=1
    Scale z (elevation) data by factor of zscale
domain='all'
    | Define a domain to be extracted from data
    | 'all': all xyz data in input file are extracted 
    | 'domain': only data within a defined domain are extracted
xmin=-180
    Minimum x (longitude) of domain to be extracted
xmax=180
    Maximum x (longitude) of domain to be extracted
ymin=-90
    Minimum y (latitude) of domain to be extracted
ymax=90
    Maximum y (latitude) of domain to be extracted
savedata='no'
    | Define if save xyz data in a file or not in outfilelocation folder
    | 'no': does not save
    | 'yes': save xyz data as csv file
outfilename='xyzdata.csv'
    | Name of output file between ' ' mark, example: 'xyzdata.csv'
    | outfilename should have '.csv' extension
outfilelocation=pwd
    Location of output file between ' ' mark, example: 'C:\'

Outputs
-------

x
    x (longitude) data extracted from xyz file
y
    y (latitude) data extracted from xyz file
z
    z (elevation) data extracted from xyz file

Examples
--------

.. code:: python

    import scientimate as sm

    xyzfilename='xyzfile.xyz' #e.g. xyzfilename='PersianGulf_ETOPO1.xyz'
    xyzfilelocation='C:/' #e.g. xyzfilelocation='C:/datafolder'
    x,y,z=sm.readxyzfile(xyzfilename,xyzfilelocation)

    xyzfilename='xyzfile.xyz' #e.g. xyzfilename='PersianGulf_ETOPO1.xyz'
    xyzfilelocation='C:/' #e.g. xyzfilelocation='C:/datafolder'
    x,y,z=sm.readxyzfile(xyzfilename,xyzfilelocation,1,'all',-180,180,-90,90,'no')

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
