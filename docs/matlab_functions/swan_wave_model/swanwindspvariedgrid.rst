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

swanwindspvariedgrid
====================

.. code:: MATLAB

    [swanwind] = swanwindspvariedgrid(xgrid, ygrid, xpointgrid, ypointgrid, windvelgrid, winddirgrid, winddirtype, windvelmin, savedata, outfilename, outfilelocation, CalcMethod)

Description
-----------

Generate SWAN wind file for spatially varied wind from gridded input data

Inputs
------

xgrid
    x (longitude) of output grid points as a [M*N] array
ygrid
    y (latitude) of output grid points as a [M*N] array
xpointgrid
    x (longitude) of the locations that wind is known in those locations as a [K*L] array
ypointgrid
    y (latitude) of the locations that wind is known in those locations as a [K*L] array
windvelgrid
    | Wind velocity at (xpointgrid,ypointgrid) as a [K*L*P] array
    | P is number of time steps for a time series
winddirgrid
    | Wind direction at (xpointgrid,ypointgrid) as a [K*L*P] array in (Degree)
    | P is number of time steps for a time series
winddirtype='mete';
    | Define wind direction type
    | 'mete': meteorological wind direction 
    |     Meteorological direction represents a direction wind comes from and is measured counter-clockwise from the North
    |     0 (degree): from North, 90 (degree): from East, 180 (degree): from South, 270 (degree): from West
    | 'trig': trigonometric wind direction
windvelmin=0;
    Minimum allowed wind velocity
savedata='no';
    | Define if save data in a file or not
    | 'no': does not save 
    | 'yes': save data as ascii file
outfilename='swanwind.wnd';
    | Name of output file between ' ' mark, example: 'swanwind.wnd'
    | outfilename should have '.wnd' extension
outfilelocation=pwd;
    Location of output file between ' ' mark, example: 'C:\' in MATLAB, or 'C:/' in Python
CalcMethod='linear';
    | Interpolation method 
    | 'linear': Use default or 'linear' method to interpolate
    | 'nearest': Use nearest neighbor method to interpolate

Outputs
-------

swanwind
    | Spatially varied wind velocity data formated for SWAN
    | Wind velocity data at each time step is assigned into the grid points

Examples
--------

.. code:: MATLAB

    [xgrid,ygrid]=meshgrid(linspace(-91,-90,100),linspace(28,30,100));
    [xpointgrid,ypointgrid]=meshgrid(linspace(-92,-89,100),linspace(27,31,100));
    windvelgrid=10+(12-10).*rand(100,100,4); %Data for 4 time steps
    winddirgrid=60+(65-60).*rand(100,100,4); %Data for 4 time steps
    winddirtype='mete';
    windvelmin=2.5;
    savedata='no';
    outfilename='swanwind.wnd';
    outfilelocation=pwd;
    CalcMethod='linear';
    [swanwind]=swanwindspvariedgrid(xgrid,ygrid,xpointgrid,ypointgrid,windvelgrid,winddirgrid,winddirtype,windvelmin,savedata,outfilename,outfilelocation,CalcMethod);


    [xgrid,ygrid]=meshgrid(linspace(-91,-90,100),linspace(28,30,100));
    [xpointgrid,ypointgrid]=meshgrid(linspace(-92,-89,100),linspace(27,31,100));
    windvelgrid=10+(12-10).*rand(100,100); %Data for 1 time step
    winddirgrid=60+(65-60).*rand(100,100); %Data for 1 time step
    winddirtype='mete';
    windvelmin=2.5;
    savedata='no';
    outfilename='swanwind.wnd';
    outfilelocation=pwd;
    CalcMethod='linear';
    [swanwind]=swanwindspvariedgrid(xgrid,ygrid,xpointgrid,ypointgrid,windvelgrid,winddirgrid,winddirtype,windvelmin,savedata,outfilename,outfilelocation,CalcMethod);

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
