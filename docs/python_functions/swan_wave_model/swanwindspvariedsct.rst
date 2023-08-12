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

scientimate.swanwindspvariedsct
===============================

.. code:: python

    swanwind, windvelpoint, winddirpoint = scientimate.swanwindspvariedsct(xgrid, ygrid, xpoint, ypoint, windvel, winddir, winddirtype='mete', windvelmin=0, savedata='no', outfilename='swanwind.wnd', outfilelocation=None, CalcMethod='linear')

Description
-----------

Generate SWAN wind file for spatially varied wind from scattered input data

Inputs
------

xgrid
    x (longitude) of output grid points as a [M*N] array
ygrid
    y (latitude) of output grid points as a [M*N] array
xpoint
    | x (longitude) of the location (such as meteorological station) that wind is known in that location
    | as a [L] array
    | 1st element is 1st point, 2nd element is 2nd point, 3rd element is 3rd point, ...
ypoint
    | y (latitude) of the location (such as meteorological station) that wind is known in that location
    | as a [L] array
    | 1st element is 1st point, 2nd element is 2nd point, 3rd element is 3rd point, ...
windvel
    | Wind velocity as a [K*L] array
    | K is number of time steps for a time series
    |     1st row (i.e. [1,:]) is 1st time step, 2nd row (i.e. [2,:]) is 2nd time step, ...
    |     Example: windvel(1,:) is the wind velocity data for all points at the first time step
    | L is number of points (such as meteorological stations) that wind velocity is known in those locations
    |     L should be L>=3
    |     1st column (i.e. [:,1]) is 1st point, 2nd column (i.e. [:,2]) is 2nd point, 3rd column (i.e. [:,3]) is 3rd point, ...
    |     Example: windvel(:,1) is the wind velocity data for all time steps at the first point
winddir
    | Wind direction as a [K*L] array in (Degree)
    | K is number of time steps for a time series
    |     1st row (i.e. [1,:]) is 1st time step, 2nd row (i.e. [2,:]) is 2nd time step, ...
    |     Example: windvel(1,:) is the wind direction data for all points at the first time step
    | L is number of points (such as meteorological stations) that wind direction is known in those locations
    |     L should be L>=3
    |     1st column (i.e. [:,1]) is 1st point, 2nd column (i.e. [:,2]) is 2nd point, 3rd column (i.e. [:,3]) is 3rd point, ...
    |     Example: windvel(:,1) is the wind direction data for all time steps at the first point
winddirtype='mete'
    | Define wind direction type
    | 'mete': meteorological wind direction 
    |     Meteorological direction represents a direction wind comes from and is measured counter-clockwise from the North
    |     0 (degree): from North, 90 (degree): from East, 180 (degree): from South, 270 (degree): from West
    | 'trig': trigonometric wind direction
windvelmin=0
    Minimum allowed wind velocity
savedata='no'
    | Define if save data in a file or not
    | 'no': does not save 
    | 'yes': save data as ascii file
outfilename='swanwind.wnd'
    | Name of output file between ' ' mark, example: 'swanwind.wnd'
    | outfilename should have '.wnd' extension
outfilelocation=pwd
    Location of output file between ' ' mark, example: 'C:\' in MATLAB, or 'C:/' in Python
CalcMethod='linear'
    | Interpolation method 
    | 'linear': Use default or 'linear' method to interpolate
    | 'nearest': Use nearest neighbor method to interpolate

Outputs
-------

swanwind
    | Spatially varied wind velocity data formated for SWAN
    | Wind velocity data at each time step is assigned into the grid points
windvelpoint
    | Nearest interpolated wind velocity at (xpoint,ypoint) as a [K*L] array
    | K is number of time steps for a time series
    | L is number of points (such as meteorological stations) that wind velocity is known in those locations
winddirpoint
    | Nearest interpolated wind direction at (xpoint,ypoint) as a [K*L] array
    | K is number of time steps for a time series
    | L is number of points (such as meteorological stations) that wind direction is known in those locations

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np

    xgrid,ygrid=np.meshgrid(np.linspace(-91,-90,100),np.linspace(28,30,100))
    windvel=[[10.5,10.55,10.6],[10.64,10.69,10.74],[10.7,10.75,10.8],[10.4,10.45,10.5]] #Data for 4 time steps
    winddir=[[50,55,60],[64,69,74],[70,75,80],[40,45,50]] #Data for 4 time steps
    xpoint=[-90.5,-90.3,-90.7] #Data are known at 3 locations
    ypoint=[29.2,29,28.8] #Data are known at 3 locations
    winddirtype='mete'
    windvelmin=2.5
    savedata='no'
    outfilename='swanwind.wnd'
    outfilelocation=None
    CalcMethod='linear'
    swanwind,windvelpoint,winddirpoint=sm.swanwindspvariedsct(xgrid,ygrid,xpoint,ypoint,windvel,winddir,winddirtype,windvelmin,savedata,outfilename,outfilelocation,CalcMethod)


    xgrid,ygrid=np.meshgrid(np.linspace(-91,-90,100),np.linspace(28,30,100))
    windvel=[10.5,10.55,10.6] #Data for 1 time step
    winddir=[50,55,60] #Data for 1 time step
    xpoint=[-90.5,-90.3,-90.7] #Data are known at 3 locations
    ypoint=[29.2,29,28.8] #Data are known at 3 locations
    winddirtype='mete'
    windvelmin=2.5
    savedata='no'
    outfilename='swanwind.wnd'
    outfilelocation=None
    CalcMethod='linear'
    swanwind,windvelpoint,winddirpoint=sm.swanwindspvariedsct(xgrid,ygrid,xpoint,ypoint,windvel,winddir,winddirtype,windvelmin,savedata,outfilename,outfilelocation,CalcMethod)

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
