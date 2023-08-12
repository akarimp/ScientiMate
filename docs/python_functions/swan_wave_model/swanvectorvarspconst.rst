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

scientimate.swanvectorvarspconst
================================

.. code:: python

    swanvectorvariable = scientimate.swanvectorvarspconst(Vx, Vy, savedata='no', outfilename='swanwind.wnd', outfilelocation=None)

Description
-----------

Generate SWAN file for spatially constant vector variable

Inputs
------

Vx
    | Variable in x direction (x component of input variable)
    | If size of Vx>1, then it is considered as a time series
    | 1st element is 1st time step, 2nd element is 2nd time step, ...
Vy
    | Variable in y direction (y component of input variable)
    | Should have a same size as Vx
savedata='no'
    | Define if save data in a file or not
    | 'no': does not save 
    | 'yes': save data as ascii file
outfilename='swanwind.wnd'
    | Name of output file between ' ' mark, example: 'swanwind.wnd'
    | outfilename should have proper name and extension
outfilelocation=pwd
    Location of output file between ' ' mark, example: 'C:\' in MATLAB, or 'C:/' in Python

Outputs
-------

swanvectorvariable
    | Spatially constant vector variable formated for SWAN
    | Note: Vector variable at each time step is assigned into 4 points,
    |     assuming the vector variable domain is defined by 4 points, one at each corner

Examples
--------

.. code:: python

    import scientimate as sm

    windvelx=[10.5,10.6,10.55] #Data for 3 time steps
    windvely=[2.5,2.6,2.55] #Data for 3 time steps
    savedata='no'
    outfilename='swanwind.wnd'
    outfilelocation=None
    swanvectorvariable=sm.swanvectorvarspconst(windvelx,windvely,savedata,outfilename,outfilelocation)

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
