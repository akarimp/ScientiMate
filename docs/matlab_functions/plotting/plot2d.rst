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

plot2d
======

.. code:: MATLAB

    plot2d(x, y, plottype, cmapcolors, sizestyle)

Description
-----------

Plot x and y data in 2-d plot

Inputs
------

x
    x data as a 2-d array of (M,N)
y=[];
    y data as a 2-d array of (M,N)
plottype='line';
    | Plot type
    | 'line': line plot
    | 'line_grid': line plot with grid lines
    | 'line_ascend': line plot with ascending line width
    | 'line_ascend_grid': line plot with ascending line width and grid lines
    | 'line_confid': line plot with 95% confidence intervals band (approximated)
    | 'line_confid_grid': line plot with 95% confidence intervals band (approximated) and grid lines
    | 'scatter': scatter plot
    | 'scatter_grid': scatter plot with grid lines
    | 'scatter_ascend': scatter plot with ascending point size
    | 'scatter_ascend_grid': scatter plot with grid lines and ascending point size
    | 'bar': bar plot
    | 'bar_grid': bar plot with grid lines
    | 'bar_stacked': stacked bar plot
    | 'barh': horizontal bar plot
    | 'barh_grid': horizontal bar plot with grid lines
    | 'barh_stacked': horizontal stacked bar plot
    | 'histogram': histogram plot
    | 'histogram_grid': histogram plot with grid lines
cmapcolors='blue';
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
sizestyle='medium';
    | Plot drawing size style
    | 'small': small plot size
    | 'medium': medium plot size
    | 'large': large plot size

Outputs
-------


Examples
--------

.. code:: MATLAB

    x=[0;1];
    y(1,:)=linspace(2,51,50);
    y(2,:)=linspace(2,51,50);
    plot2d(x,y,'line','blue_red','large')

    x(:,1)=linspace(1,10,10);
    y(:,1)=1+rand(10,1);
    y(:,2)=2+rand(10,1);
    plot2d(x,y,'line_confid','blue_red','large')

    x(:,1)=linspace(0,2*pi,1000);
    x=repmat(x,1,10);
    s=[1:10];
    s=repmat(s,1000,1);
    y=s+sin(1.0*pi.*x);
    plot2d(x,y,'line','cool','large')

    x=rand(100,3);
    y=rand(100,3);
    y(:,1)=1+2.0.*x(:,1)+rand(100,1);
    y(:,2)=3+2.0.*x(:,2)+rand(100,1);
    y(:,3)=5+2.0.*x(:,3)+rand(100,1);
    plot2d(x,y,'scatter','seq','large')

    x=rand(100,1);
    y=rand(100,1);
    plot2d(x,y,'scatter_ascend','purple','large')

    x=[1,1,1;2,2,2;3,3,3;4,4,4];
    y=[2,3,8;2,5,6;5,7,9;1,2,3];
    plot2d(x,y,'bar','purple','medium')

    x=[1,3,5,7,9,11,13,15];
    y=[2,3,9,8,2,5,6,9];
    plot2d(x,y,'bar','purple','medium')

    x=randn(1000,1);
    plot2d(x,[],'histogram','purple','medium')

References
----------

Colormap

* https://matplotlib.org/tutorials/colors/colormaps.html
* https://www.mathworks.com/help/matlab/ref/colormap.html
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
