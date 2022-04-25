.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-07-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

windavg
=======

.. code:: MATLAB

    [windvelavg, winddiravg] = windavg(windvel, winddir, NPointsAvg, NPointsInterval, dispout)

Description
-----------

Average wind velocity and wind direction

Inputs
------

windvel
    Wind velocity time series data
winddir
    Wind direction time series data in (degree)
NPointsAvg=length(windvel(:,1));
    Number of data points from start of each section (interval) to be averaged
NPointsInterval=length(windvel(:,1));
    Number of points that each section (interval) has
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

windvelavg
    Averaged wind velocity data
winddiravg
    Averaged wind direction data in (degree)

Examples
--------

.. code:: MATLAB

    windvel(:,1)=10.*rand(10,1);
    winddir(:,1)=45.*rand(10,1);
    [windvelavg,winddiravg]=windavg(windvel,winddir);

    windvel(:,1)=10.*rand(5*60,1); %One data point every minute for 5 hours
    winddir(:,1)=225.*rand(5*60,1); %One data point every minute for 5 hours
    [windvelavg,winddiravg]=windavg(windvel,winddir,10,60,'yes');

References
----------

Yamartino, R. J. (1984). 
A comparison of several "single-pass" estimators of the standard deviation of wind direction. 
Journal of Climate and Applied Meteorology, 23(9), 1362-1366.

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
