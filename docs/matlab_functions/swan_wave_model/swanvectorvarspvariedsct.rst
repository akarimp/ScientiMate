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

swanvectorvarspvariedsct
========================

.. code:: MATLAB

    [swanvectorvariable, Vxpoint, Vypoint] = swanvectorvarspvariedsct(xgrid, ygrid, xpoint, ypoint, Vx, Vy, savedata, outfilename, outfilelocation, CalcMethod)

Description
-----------

Generate SWAN file for spatially varied vector variable from scattered input data

Inputs
------

xgrid
    x (longitude) of output grid points as a [M*N] array
ygrid
    y (latitude) of output grid points as a [M*N] array
xpoint
    | x (longitude) of the location (such as meteorological station) that vector variable is known in that location
    | as a [L] array
    | 1st element is 1st point, 2nd element is 2nd point, 3rd element is 3rd point, ...
ypoint
    | y (latitude) of the location (such as meteorological station) that vector variable is known in that location
    | as a [L] array
    | 1st element is 1st point, 2nd element is 2nd point, 3rd element is 3rd point, ...
Vx
    | Variable in x direction (x component of input variable) as a [K*L] array
    | K is number of time steps for a time series
    |     1st row (i.e. [1,:]) is 1st time step, 2nd row (i.e. [2,:]) is 2nd time step, ...
    |     Example: windvel(1,:) is the wind velocity data for all points at the first time step
    | L is number of points (such as meteorological stations) that wind velocity is known in those locations
    |     L should be L>=3
    |     1st column (i.e. [:,1]) is 1st point, 2nd column (i.e. [:,2]) is 2nd point, 3rd column (i.e. [:,3]) is 3rd point, ...
    |     Example: windvel(:,1) is the wind velocity data for all time steps at the first point
Vy
    | Variable in y direction (y component of input variable) as a [K*L] array
    | K is number of time steps for a time series
    |     1st row (i.e. [1,:]) is 1st time step, 2nd row (i.e. [2,:]) is 2nd time step, ...
    |     Example: windvel(1,:) is the wind velocity data for all points at the first time step
    | L is number of points (such as meteorological stations) that wind velocity is known in those locations
    |     L should be L>=3
    |     1st column (i.e. [:,1]) is 1st point, 2nd column (i.e. [:,2]) is 2nd point, 3rd column (i.e. [:,3]) is 3rd point, ...
    |     Example: windvel(:,1) is the wind velocity data for all time steps at the first point
savedata='no';
    | Define if save data in a file or not
    | 'no': does not save 
    | 'yes': save data as ascii file
outfilename='swanwind.wnd';
    | Name of output file between ' ' mark, example: 'swanwind.wnd'
    | outfilename should have proper name and extension
outfilelocation=pwd;
    Location of output file between ' ' mark, example: 'C:\' in MATLAB, or 'C:/' in Python
CalcMethod='linear';
    | Interpolation method 
    | 'linear': Use default or 'linear' method to interpolate
    | 'nearest': Use nearest neighbor method to interpolate

Outputs
-------

swanvectorvariable
    | Spatially varied vector variable data formated for SWAN
    | Vector variable data at each time step is assigned into the grid points
Vxpoint
    | Nearest interpolated vector variable in x direction at (xpoint,ypoint) as a [K*L] array
    | K is number of time steps for a time series
    | L is number of points (such as meteorological stations) that vector variable in x direction is known in those locations
Vypoint
    | Nearest interpolated vector variable in y direction at (xpoint,ypoint) as a [K*L] array
    | K is number of time steps for a time series
    | L is number of points (such as meteorological stations) that vector variable in y direction is known in those locations

Examples
--------

.. code:: MATLAB

    [xgrid,ygrid]=meshgrid(linspace(-91,-90,100),linspace(28,30,100));
    windvelx=[[10.5,10.55,10.6];[10.64,10.69,10.74];[10.7,10.75,10.8];[10.4,10.45,10.5]]; %Data for 4 time steps
    windvely=[[1.5,1.55,1.6];[1.64,1.69,1.74];[1.7,1.75,1.8];[1.4,1.45,1.5]]; %Data for 4 time steps
    xpoint=[-90.5;-90.3;-90.7]; %Data are known at 3 locations
    ypoint=[29.2;29;28.8]; %Data are known at 3 locations
    savedata='no';
    outfilename='swanwind.wnd';
    outfilelocation=pwd;
    CalcMethod='linear';
    [swanvectorvariable,windvelxpoint,windvelypoint]=swanvectorvarspvariedsct(xgrid,ygrid,xpoint,ypoint,windvelx,windvely,savedata,outfilename,outfilelocation,CalcMethod);


    [xgrid,ygrid]=meshgrid(linspace(-91,-90,100),linspace(28,30,100));
    windvelx=[10.5,10.55,10.6]; %Data for 1 time step
    windvely=[1.5,1.55,1.6]; %Data for 1 time step
    xpoint=[-90.5;-90.3;-90.7]; %Data are known at 3 locations
    ypoint=[29.2;29;28.8]; %Data are known at 3 locations
    savedata='no';
    outfilename='swanwind.wnd';
    outfilelocation=pwd;
    CalcMethod='linear';
    [swanvectorvariable,windvelxpoint,windvelypoint]=swanvectorvarspvariedsct(xgrid,ygrid,xpoint,ypoint,windvelx,windvely,savedata,outfilename,outfilelocation,CalcMethod);

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
