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

scientimate.swanwaterlevelspvariedgrid
======================================

.. code:: python

    swanwaterlevel = scientimate.swanwaterlevelspvariedgrid(xgrid, ygrid, xpointgrid, ypointgrid, waterlevelgrid, savedata='no', outfilename='swanwaterlevel.wl', outfilelocation=None, CalcMethod='linear')

Description
-----------

| Generate SWAN water level file for spatially varied water level from gridded input data
| This function can be used for any other scalar variable as well

Inputs
------

xgrid
    x (longitude) of output grid points as a [M*N] array
ygrid
    y (latitude) of output grid points as a [M*N] array
xpointgrid
    x (longitude) of the locations that water level (or any scalar variable) is known in those locations as a [K*L] array
ypointgrid
    y (latitude) of the locations that water level (or any scalar variable) is known in those locations as a [K*L] array
waterlevelgrid
    | Water level (or any scalar variable) at (xpointgrid,ypointgrid) as a [K*L*P] array
    | P is number of time steps for a time series
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

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np

    xgrid,ygrid=np.meshgrid(np.linspace(-91,-90,100),np.linspace(28,30,100))
    xpointgrid,ypointgrid=np.meshgrid(np.linspace(-92,-89,100),np.linspace(27,31,100))
    waterlevelgrid=0.5+(0.6-0.5)*np.random.rand(100,100,4) #Data for 4 time steps
    savedata='no'
    outfilename='swanwaterlevel.wl'
    outfilelocation=None
    CalcMethod='linear'
    swanwaterlevel=sm.swanwaterlevelspvariedgrid(xgrid,ygrid,xpointgrid,ypointgrid,waterlevelgrid,savedata,outfilename,outfilelocation,CalcMethod)


    xgrid,ygrid=np.meshgrid(np.linspace(-91,-90,100),np.linspace(28,30,100))
    xpointgrid,ypointgrid=np.meshgrid(np.linspace(-92,-89,100),np.linspace(27,31,100))
    waterlevelgrid=0.5+(0.6-0.5)*np.random.rand(100,100) #Data for 1 time step
    savedata='no'
    outfilename='swanwaterlevel.wl'
    outfilelocation=None
    CalcMethod='linear'
    swanwaterlevel=sm.swanwaterlevelspvariedgrid(xgrid,ygrid,xpointgrid,ypointgrid,waterlevelgrid,savedata,outfilename,outfilelocation,CalcMethod)

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
