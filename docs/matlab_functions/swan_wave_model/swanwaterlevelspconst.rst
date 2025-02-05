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

swanwaterlevelspconst
=====================

.. code:: MATLAB

    [swanwaterlevel] = swanwaterlevelspconst(waterlevel, savedata, outfilename, outfilelocation)

Description
-----------

| Generate SWAN water level file for spatially constant water level
| This function can be used for any other scalar variable as well

Inputs
------

waterlevel
    | Water level (or any scalar variable)
    | If size of waterlevel>1, then it is considered as a time series
    | 1st element is 1st time step, 2nd element is 2nd time step, ...
savedata='no';
    | Define if save data in a file or not
    | 'no': does not save 
    | 'yes': save data as ascii file
outfilename='swanwaterlevel.wl';
    | Name of output file between ' ' mark, example: 'swanwaterlevel.wl'
    | outfilename should have '.wl' extension
    | In case of using other scalar varable than water level, use proper name and extension
outfilelocation=pwd;
    Location of output file between ' ' mark, example: 'C:\' in MATLAB, or 'C:/' in Python

Outputs
-------

swanwaterlevel
    | Spatially constant water level data (or any scalar variable) formated for SWAN
    | Note: Water level at each time step is assigned into 4 points, 
    |     assuming the water level domain is defined by 4 points, one at each corner

Examples
--------

.. code:: MATLAB

    waterlevel=[0.5;0.6;0.55]; %Data for 3 time steps
    savedata='no';
    outfilename='swanwaterlevel.wl';
    outfilelocation=pwd;
    [swanwaterlevel]=swanwaterlevelspconst(waterlevel,savedata,outfilename,outfilelocation);

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
