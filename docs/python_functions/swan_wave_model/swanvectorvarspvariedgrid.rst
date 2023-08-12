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

scientimate.swanvectorvarspvariedgrid
=====================================

.. code:: python

    swanvectorvariable = scientimate.swanvectorvarspvariedgrid(xgrid, ygrid, xpointgrid, ypointgrid, Vxgrid, Vygrid, savedata='no', outfilename='swanwind.wnd', outfilelocation=None, CalcMethod='linear')

Description
-----------

Generate SWAN file for spatially varied vector variable from gridded input data

Inputs
------

xgrid
    x (longitude) of output grid points as a [M*N] array
ygrid
    y (latitude) of output grid points as a [M*N] array
xpointgrid
    x (longitude) of the locations that vector variable is known in those locations as a [K*L] array
ypointgrid
    y (latitude) of the locations that vector variable is known in those locations as a [K*L] array
Vxgrid
    | Variable in x direction (x component of input variable) at (xpointgrid,ypointgrid) as a [K*L*P] array
    | P is number of time steps for a time series
Vygrid
    | Variable in y direction (y component of input variable) at (xpointgrid,ypointgrid) as a [K*L*P] array in (Degree)
    | P is number of time steps for a time series
savedata='no'
    | Define if save data in a file or not
    | 'no': does not save 
    | 'yes': save data as ascii file
outfilename='swanwind.wnd'
    | Name of output file between ' ' mark, example: 'swanwind.wnd'
    | outfilename should have proper name and extension
outfilelocation=pwd
    Location of output file between ' ' mark, example: 'C:\' in MATLAB, or 'C:/' in Python
CalcMethod='linear'
    | Interpolation method 
    | 'linear': Use default or 'linear' method to interpolate
    | 'nearest': Use nearest neighbor method to interpolate

Outputs
-------

swanvectorvariable
    | Spatially varied vector variable data formated for SWAN
    | Vector variable data at each time step is assigned into the grid points

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np

    xgrid,ygrid=np.meshgrid(np.linspace(-91,-90,100),np.linspace(28,30,100))
    xpointgrid,ypointgrid=np.meshgrid(np.linspace(-92,-89,100),np.linspace(27,31,100))
    windvelxgrid=10+(12-10)*np.random.rand(100,100,4) #Data for 4 time steps
    windvelygrid=1+(2-1)*np.random.rand(100,100,4) #Data for 4 time steps
    savedata='no'
    outfilename='swanwind.wnd'
    outfilelocation=None
    CalcMethod='linear'
    swanvectorvariable=sm.swanvectorvarspvariedgrid(xgrid,ygrid,xpointgrid,ypointgrid,windvelxgrid,windvelygrid,savedata,outfilename,outfilelocation,CalcMethod)


    xgrid,ygrid=np.meshgrid(np.linspace(-91,-90,100),np.linspace(28,30,100))
    xpointgrid,ypointgrid=np.meshgrid(np.linspace(-92,-89,100),np.linspace(27,31,100))
    windvelxgrid=10+(12-10)*np.random.rand(100,100) #Data for 1 time step
    windvelygrid=1+(2-1)*np.random.rand(100,100) #Data for 1 time step
    savedata='no'
    outfilename='swanwind.wnd'
    outfilelocation=None
    CalcMethod='linear'
    swanvectorvariable=sm.swanvectorvarspvariedgrid(xgrid,ygrid,xpointgrid,ypointgrid,windvelxgrid,windvelygrid,savedata,outfilename,outfilelocation,CalcMethod)

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
