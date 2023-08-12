.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2020-02-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

plot3dhillshades
================

.. code:: MATLAB

    plot3dhillshades(x, y, z, plottype)

Description
-----------

Plot hillshades (shaded relief) of x (longitude), y (latitude) and z (elevation) data

Inputs
------

x
    | x (longitude) data extracted from xyz file
    | Set x=[] if it is not available
    | It may be 1d or 2d array
y
    | y (latitude) data extracted from xyz file
    | Set y=[] if it is not available
    | It may be 1d or 2d array
z
    z (elevation) data extracted from xyz file
plottype='imagesc';
    | Plot type
    | 'imagesc': 2 dimensional plot using imagesc or imshow
    | 'pcolor': 2 dimensional plot using pcolor
    | 'contour': 2 dimensional contour plot, number of contour=ncolor
    | 'surface': 3 dimensional surface plot 

Outputs
-------


Examples
--------

.. code:: MATLAB

    [x,y]=meshgrid(linspace(-10,10,50),linspace(-10,10,50));
    r=sqrt(x.^2+y.^2)+1e-10; %Add 1e-10 to prevent divide by 0
    z=sin(r)./r;
    plot3dhillshades(x,y,z,'imagesc')

    [x,y]=meshgrid(linspace(-10,10,21),linspace(-10,10,21));
    z=(sin(x)+sin(y))./(x+y+1e-10); %Add 1e-10 to prevent divide by 0
    plot3dhillshades(x,y,z,'imagesc')

    x(:,1)=10.*rand(1000,1);
    y(:,1)=10.*rand(1000,1);
    z=x.^2+y.^2;
    plot3dhillshades(x,y,z,'imagesc')


References
----------

* https://pro.arcgis.com/en/pro-app/tool-reference/3d-analyst/how-hillshade-works.htm
* https://pro.arcgis.com/en/pro-app/tool-reference/3d-analyst/applying-a-z-factor.htm
* http://mike.teczno.com/img/hillshade.py
* http://geospatialpython.com/2013/12/python-and-elevation-data-creating.html
* https://github.com/ThomasLecocq/geophysique.be/blob/master/2014-02-25%20Shaded%20Relief%20Map%20in%20Python.ipynb
* https://www.geophysique.be/2014/02/25/shaded-relief-map-in-python/
* https://www.mathworks.com/matlabcentral/fileexchange/14863-hillshade
* http://www.reliefshading.com/

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
