.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-08-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

zprofilepath
============

.. code:: MATLAB

    [z, zmean, distxy] = zprofilepath(xgrid, ygrid, zgrid, x, y, distCalcMethod, CalcMethod, dispout)

Description
-----------

Calculate z (elevation, ...) profile along a path over a given 2d x-y domain (map, image, ...)

Inputs
------

xgrid
    x (longitude, pixel, ...) data as a [M*N] array
ygrid
    y (latitude, pixel, ...) data as a [M*N] array
zgrid
    z (elevation, ...) data as a [M*N] array
x
    x (longitude, pixel, ...) of the points along the line as (x,y)
y
    | y (latitude, pixel, ...) of the points along the line as (x,y)
    | If input data are latitude and longitude, they should be in Degree
distCalcMethod='cart';
    | Distance calculation method 
    | 'cart': Distances are calculated on cartesian coordinate
    | 'gc': Distances are calculated on Great Circle based on Vincenty formula, Vincenty (1975)
    | Earth radius coonsidered as mean earth radius=6371000 m
CalcMethod='nearest';
    | Interpolation method 
    | 'linear': Use default or 'linear' method to interpolate
    | 'nearest': Use nearest neighbor method to interpolate
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

z
    z (elevation, ...) data along a path at given points (x,y)
zmean
    Weighted mean of z (elevation) along a line calculated az zmean=1/(total_distance)*(sum(z(i)*dx(i)))
distxy
    | Distance at each points of (x,y) from the first point of (x(1),y(1))
    | If input data are latitude and longitude in Degree, distxy is in m

Examples
--------

.. code:: MATLAB

    [xgrid,ygrid]=meshgrid(linspace(-5,5,100),linspace(-5,5,100));
    zgrid=sin(xgrid.^2+ygrid.^2)./(xgrid.^2+ygrid.^2);
    x1=-3;
    y1=-3;
    x2=3;
    y2=3;
    x(:,1)=linspace(x1,x2,1000);
    y(:,1)=linspace(y1,y2,1000);
    [z,zmean,distxy]=zprofilepath(xgrid,ygrid,zgrid,x,y,'cart','nearest','yes');

References
----------

Vincenty, T. (1975). 
Direct and inverse solutions of geodesics on the ellipsoid with application of nested equations. 
Survey review, 23(176), 88-93.

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
