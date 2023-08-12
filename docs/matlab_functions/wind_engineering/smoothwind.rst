.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2021-02-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

smoothwind
==========

.. code:: MATLAB

    [windvel_smooth, winddir_smooth] = smoothwind(windvel, winddir, window_length, dispout)

Description
-----------

Smooth wind data using moving average window

Inputs
------

windvel=[];
    | Wind velocity time series data
    | Leave wind velocity empty if only winddir is available
winddir=[];
    | Wind direction time series data
    | Leave wind direction empty if only windvel is available
window_length=5;
    Number of data points in sliding window for calculating moving average values, must be odd
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

windvel_smooth
    Smoothed wind velocity
winddir_smooth
    | Smoothed wind direction
    | Wind directions with large variations are not smoothed

Examples
--------

.. code:: MATLAB

    %Data from https://tidesandcurrents.noaa.gov for Grand Isle, LA, USA (8761724), for June 1st, 2017, reported hourly
    windvel=[3;4.7;4.9;5.3;3.3;3.4;3.3;3.8;3.1;2;1.3;1.2;1.5;3.2;2.9;3;2.9;3.7;3.7;3.1;3.4;2.6;2.5;2.5]; %24 Hour wind velocity
    winddir=[78;86;88;107;131;151;163;163;153;150;148;105;105;75;95;103;97;103;108;111;124;183;171;113]; %24 Hour wind direction
    [windvel_smooth,winddir_smooth]=smoothwind(windvel,winddir,5,'yes');

References
----------

U.S. Army Corps of Engineers (2015). 
Coastal Engineering Manual. 
Engineer Manual 1110-2-1100, Washington, D.C.: U.S. Army Corps of Engineers.

Yamartino, R. J. (1984). 
A comparison of several "single-pass" estimators of the standard deviation of wind direction. 
Journal of Climate and Applied Meteorology, 23(9), 1362-1366.

.. License & Disclaimer
.. --------------------
..
.. Copyright (c) 2021 Arash Karimpour
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
