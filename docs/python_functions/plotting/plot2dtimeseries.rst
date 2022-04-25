.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2019-02-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

scientimate.plot2dtimeseries
============================

.. code:: python

    scientimate.plot2dtimeseries(x, starttime='2000-10-20 00:00:00', endtime='2100-10-20 23:00:00', plottype='line', cmapcolors='blue', sizestyle='medium')

Description
-----------

Plot x data in 2-d timeseries

Inputs
------

x
    x data as a 2-d array of (M,N)
starttime='2000-10-20 00:00:00'
    | Start date and time for time series generation
    | Format: 'yyyy-mm-dd HH:MM:SS'
endtime='2100-10-20 23:00:00'
    | End date and time for time series generation
    | Format: 'yyyy-mm-dd HH:MM:SS'
plottype='bar'
    | Plot type
    | 'line': line plot
    | 'line_grid': line plot with grid lines
    | 'line_subplot': line plot with subplots
    | 'line_subplot_grid': line plot with subplots and grid lines
    | 'bar': bar plot
    | 'bar_grid': bar plot with grid lines
    | 'barh': horizontal bar plot
    | 'barh_grid': horizontal bar plot with grid lines
cmapcolors='seq'
    | Colormap style
    | 'blue': blue colormap
    | 'red': red colormap
    | 'green': green colormap
    | 'yellow': yellow colormap
    | 'purple': purple colormap
    | 'brown': brown colormap
    | 'gray': gray colormap
    | 'blue_red': blue-red colormap
    | 'red_blue': red-blue colormap
    | 'blue_green': blue-green colormap
    | 'green_blue': green-blue colormap
    | 'green_yellow': green-yellow colormap
    | 'yellow_green': yellow-green colormap
    | 'red_yellow': red-yellow colormap
    | 'yellow_red': yellow-red colormap
    | 'cyclic': cyclic/oscillation colormap 
    | 'seq': sequential colormap
    | User-defined colors may be used to generate colormap
    | User-defined colors should be defined as (M,3) array in RGB color format
    | At least two colors should be defined, i.e. M should be equal or larger than 2
    | User-defined colors values should be between 0 and 255
    | Any available colormap name such as 'cool', 'winter', etc also can be used
sizestyle='medium'
    | Plot drawing size style
    | 'small': small plot size
    | 'medium': medium plot size
    | 'large': large plot size

Outputs
-------


Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np

    x=np.random.rand(366)
    starttime='2020-01-01 00:00:00'
    endtime='2020-12-31 00:00:00' 
    sm.plot2dtimeseries(x,starttime,endtime,'line','purple','medium')

    x=np.random.rand(31,3)
    starttime='2020-01-01 00:00:00'
    endtime='2020-01-31 00:00:00'
    sm.plot2dtimeseries(x,starttime,endtime,'line_subplot','seq','medium')

References
----------

Colormap

* https://matplotlib.org/tutorials/colors/colormaps.html
* https://www.mathworks.com/help/PYTHON/ref/colormap.html
* http://colorbrewer2.org
* http://matplotlib.org/cmocean/
* http://jdherman.github.io/colormap/

Color

* http://htmlcolorcodes.com

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
