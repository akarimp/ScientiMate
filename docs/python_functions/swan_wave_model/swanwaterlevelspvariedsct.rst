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

scientimate.swanwaterlevelspvariedsct
=====================================

.. code:: python

    swanwaterlevel, waterlevelpoint = scientimate.swanwaterlevelspvariedsct(xgrid, ygrid, xpoint, ypoint, waterlevel, savedata='no', outfilename='swanwaterlevel.wl', outfilelocation=None, CalcMethod='linear')

Description
-----------

| Generate SWAN water level file for spatially varied water level from scattered input data
| This function can be used for any other scalar variable as well

Inputs
------

xgrid
    x (longitude) of output grid points as a [M*N] array
ygrid
    y (latitude) of output grid points as a [M*N] array
xpoint
    | x (longitude) of the location (such as meteorological station) that water level (or any scalar variable) is known in that location
    | as a [L] array
    | 1st element is 1st point, 2nd element is 2nd point, 3rd element is 3rd point, ...
ypoint
    | y (latitude) of the location (such as meteorological station) that water level (or any scalar variable) is known in that location
    | as a [L] array
    | 1st element is 1st point, 2nd element is 2nd point, 3rd element is 3rd point, ...
waterlevel
    | Water level (or any scalar variable) at (xpoint,ypoint) as a [K*L] array
    | K is number of time steps for a time series
    |     1st row (i.e. [1,:]) is 1st time step, 2nd row (i.e. [2,:]) is 2nd time step, ...
    |     Example: waterlevel(1,:) is the water level data for all points at the first time step
    | L is number of points (such as meteorological stations) that water level (or any scalar variable) is known in those locations
    |     L should be L>=3
    |     1st column (i.e. [:,1]) is 1st point, 2nd column (i.e. [:,2]) is 2nd point, 3rd column (i.e. [:,3]) is 3rd point, ...
    |     Example: waterlevel(:,1) is the water level data for all time steps at the first point
savedata='no'
    | Define if save data in a file or not
    | 'no': does not save 
    | 'yes': save data as ascii file
outfilename='swanwaterlevel.wl'
    | Name of output file between ' ' mark, example: 'swanwaterlevel.wl'
    | outfilename should have '.wl' extension
    | In case of using other scalar varable than water level, use proper name and extension
outfilelocation=pwd
    Location of output file between ' ' mark, example: 'C:\' in MATLAB, or 'C:/' in Python
CalcMethod='linear'
    | Interpolation method 
    | 'linear': Use default or 'linear' method to interpolate
    | 'nearest': Use nearest neighbor method to interpolate

Outputs
-------

swanwaterlevel
    | Spatially varied water level data (or any scalar variable) formated for SWAN
    | Water level data at each time step is assigned into the grid points
waterlevelpoint
    | Nearest interpolated water level (or any scalar variable) at (xpoint,ypoint) as a [K*L] array
    | K is number of time steps for a time series
    | L is number of points (such as meteorological stations) that water level (or any scalar variable) is known in those locations

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np

    xgrid,ygrid=np.meshgrid(np.linspace(-91,-90,100),np.linspace(28,30,100))
    xpoint=[-90.5,-90.3,-90.7] #Data are known at 3 locations
    ypoint=[29.2,29,28.8] #Data are known at 3 locations
    waterlevel=[[0.5,0.55,0.6],[0.64,0.69,0.74],[0.7,0.75,0.8],[0.4,0.45,0.5]] #Data for 4 time steps
    savedata='no'
    outfilename='swanwaterlevel.wl'
    outfilelocation=None
    CalcMethod='linear'
    swanwaterlevel,waterlevelpoint=sm.swanwaterlevelspvariedsct(xgrid,ygrid,xpoint,ypoint,waterlevel,savedata,outfilename,outfilelocation,CalcMethod)


    xgrid,ygrid=np.meshgrid(np.linspace(-91,-90,100),np.linspace(28,30,100))
    xpoint=[-90.5,-90.3,-90.7] #Data are known at 3 locations
    ypoint=[29.2,29,28.8] #Data are known at 3 locations
    waterlevel=[0.5,0.55,0.6] #Data for 1 time step
    savedata='no'
    outfilename='swanwaterlevel.wl'
    outfilelocation=None
    CalcMethod='linear'
    swanwaterlevel,waterlevelpoint=sm.swanwaterlevelspvariedsct(xgrid,ygrid,xpoint,ypoint,waterlevel,savedata,outfilename,outfilelocation,CalcMethod)

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
