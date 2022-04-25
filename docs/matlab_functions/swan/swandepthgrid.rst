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

swandepthgrid
=============

.. code:: MATLAB

    [swandepth, swangrid, ncellx, ncelly, ngridx, ngridy] = swandepthgrid(xgrid, ygrid, zgrid, zmin, zmax, nanreplacement, savedata, xyoutfilename, zoutfilename, outfilelocation, dispout)

Description
-----------

Generate SWAN depth file and its associated x-y grid file

Inputs
------

xgrid
    x (longitude) of grid points as a [M*N] array
ygrid
    y (latitude) of grid points as a [M*N] array
zgrid
    | z (elevation) of grid points as a [M*N] array
    | z>0 is land, z<0 is water, z=0 is water surface
zmin=nanmin(zgrid);
    | Minimum z (elevation) to be considered
    | All z<zmin would be set to NaN
zmax=nanmax(zgrid);
    | Maximum z (elevation) to be considered
    | All z>zmax would be set to NaN
    | Example: to remove all land values from z data, set zmax=0
nanreplacement='no';
    | Replace NaN values with nanreplacement
    | By seting up an exception value equal to nanreplacement in SWAN input file, SWAN disregards these data points
    | If nanreplacement='no', then it does not replace NaN values
    | Example, nanreplacement=-999 will replace all NaN values with -999
    |     Then if -999 is set as an exception in SWAN input file, SWAN disregards all -999 values
savedata='no';
    | Define if save data in a file or not
    | 'no': does not save 
    | 'yes': save data as ascii file
xyoutfilename='swangrid.xy';
    | Name of output grid file between ' ' mark, example: 'swangrid.xy'
    | xyoutfilename should have '.xy' extension
zoutfilename='swandepth.dep';
    | Name of output depth file between ' ' mark, example: 'swandepth.dep'
    | zoutfilename should have '.dep' extension
outfilelocation=pwd;
    Location of output file between ' ' mark, example: 'C:\' in MATLAB, or 'C:/' in Python
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

swandepth
    | Depth data formated for SWAN
    | Note: SWAN (Delft3D) needs depth data not bathymetry data. 
    |     It means value below water surface should be positive and above surface negative.
    | Note: All NaN values are replaced with nanreplacement
    |     Set up an exception value equal to nanreplacement in SWAN to disregard these data points
swangrid
    Grid data formated for SWAN
ncellx
    | Number of cells in x direction (ncellx=ngridx-1)
    | In SWAN, ncellx is equal to a number of meshes in computational grid in x direction 
ncelly
    | Number of cells in y direction (ncelly=ngridy-1)
    | In SWAN, ncelly is equal to a number of meshes in computational grid in y direction 
ngridx
    Number of grid points in x direction
ngridy
    Number of grid points in y direction

Examples
--------

.. code:: MATLAB

    x(:,1)=10.*rand(1000,1);
    y(:,1)=10.*rand(1000,1);
    z=x.^2+y.^2-mean(x.^2+y.^2);
    [xgrid,ygrid]=meshgrid(linspace(min(x),max(x),100),linspace(min(y),max(y),100));
    zgrid=griddata(x,y,z,xgrid,ygrid);
    zmin=nanmin(zgrid(:));
    zmax=nanmax(zgrid(:));
    nanreplacement=-999;
    savedata='no';
    xyoutfilename='swangrid.xy';
    zoutfilename='swandepth.dep';
    outfilelocation=pwd;
    [swandepth,swangrid,ncellx,ncelly,ngridx,ngridy]=swandepthgrid(xgrid,ygrid,zgrid,zmin,zmax,nanreplacement,savedata,xyoutfilename,zoutfilename,outfilelocation,'yes');

References
----------

Booij, N. R. R. C., Ris, R. C., & Holthuijsen, L. H. (1999). 
A third‐generation wave model for coastal regions: 1. Model description and validation. 
Journal of geophysical research: Oceans, 104(C4), 7649-7666.

SWAN Team. (2007). S
WAN user manual. 
Delft University of Technology. The Netherlands.

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
