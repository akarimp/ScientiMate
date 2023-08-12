.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-10-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

interpxyz2xyz
=============

.. code:: MATLAB

    [zPoint] = interpxyz2xyz(x, y, z, xPoint, yPoint, interpMethod, dispout)

Description
-----------

| Interpolate 2d scattered data on given point(s) by down-sampling the input data 
| For down-sampling, the first point of the k-nearest neighbors calculated from Euclidean distance is used
| interpxyz2xyz is more suitable for large dataset
| griddata is more efficient than interpxyz2xyz for regular size dataset

Inputs
------

x
    Coordinate of data in x direction, as an 1D array
y
    Coordinate of data in y direction, as an 1D array
z
    Value of data at (x,y) as z(x,y), as an 1D array
xPoint
    Coordinate (in x direction) of point that nearest point to that is desired to be found 
yPoint
    Coordinate (in y direction) of point that nearest point to that is desired to be found 
interpMethod='nearest';
    | Interpolation method 
    | 'linear': Use default or 'linear' method to interpolate
    | 'nearest': Use nearest neighbor method to interpolate
dispout='no';
    | Define to display outputs or not ('yes': display, 'no': not display)
    | '2d': 2 dimensional scatter plot 
    | 'surface': 3 dimensional surface plot 
    | 'no': not display 

Outputs
-------

zPoint
    Value of interpolated data at (xPoint,yPoint) as z(xPoint,yPoint) 

Examples
--------

.. code:: MATLAB

    x=10.*rand(100,1);
    y=10.*rand(100,1);
    z=100.*rand(100,1);
    xPoint=[2.5;5;7.5];
    yPoint=[3;6;9];
    [zPoint]=interpxyz2xyz(x,y,z,xPoint,yPoint,'linear','2d');

    x=10.*rand(100,1);
    y=10.*rand(100,1);
    z=100.*rand(100,1);
    [xgrid,ygrid]=meshgrid(linspace(min(x(:)),max(x(:)),100),linspace(min(y(:)),max(y(:)),100));
    [zgrid]=interpxyz2xyz(x,y,z,xgrid,ygrid,'nearest','no');

    x=10.*rand(100,1);
    y=10.*rand(100,1);
    z=y.*sin(x)-x.*cos(y);
    xPoint=10.*rand(10,1);
    yPoint=10.*rand(10,1);
    [zPoint]=interpxyz2xyz(x,y,z,xPoint,yPoint,'nearest','no');

    x=10.*rand(100,1);
    y=10.*rand(100,1);
    z=y.*sin(x)-x.*cos(y);
    [xgrid,ygrid]=meshgrid(linspace(min(x(:)),max(x(:)),100),linspace(min(y(:)),max(y(:)),100));
    [zgrid]=interpxyz2xyz(x,y,z,xgrid,ygrid,'nearest','surface');

References
----------


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
