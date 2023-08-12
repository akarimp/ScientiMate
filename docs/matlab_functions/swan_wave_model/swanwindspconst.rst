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

swanwindspconst
===============

.. code:: MATLAB

    [swanwind] = swanwindspconst(windvel, winddir, winddirtype, windvelmin, savedata, outfilename, outfilelocation)

Description
-----------

Generate SWAN wind file for spatially constant wind

Inputs
------

windvel
    | Wind velocity
    | If size of windvel>1, then it is considered as a time series
    | 1st element is 1st time step, 2nd element is 2nd time step, ...
winddir
    Wind direction in (Degree)
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

Outputs
-------

swanwind
    | Spatially constant wind data formated for SWAN
    | Note: Wind at each time step is assigned into 4 points,
    |     assuming the wind domain is defined by 4 points, one at each corner

Examples
--------

.. code:: MATLAB

    windvel=[10.5;10.6;10.55]; %Data for 3 time steps
    winddir=[30;32;28]; %Data for 3 time steps
    winddirtype='mete';
    windvelmin=2.5;
    savedata='no';
    outfilename='swanwind.wnd';
    outfilelocation=pwd;
    [swanwind]=swanwindspconst(windvel,winddir,winddirtype,windvelmin,savedata,outfilename,outfilelocation);

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
