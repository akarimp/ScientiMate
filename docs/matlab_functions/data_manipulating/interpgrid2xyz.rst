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

interpgrid2xyz
==============

.. code:: MATLAB

    [z] = interpgrid2xyz(xgrid, ygrid, zgrid, x, y, dispout)

Description
-----------

Interpolate 2d gridded data on given scatter point(s) using nearest neighbor method 

Inputs
------

xgrid
    x data as a [M*N] array
ygrid
    y data as a [M*N] array
zgrid
    z data as z(x,y) as a [M*N] array
x
    x of point that nearest point to that is desired to be found
y
    y of point that nearest point to that is desired to be found
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

z
    Value of interpolated data at (x,y) as z(x,y)

Examples
--------

.. code:: MATLAB

    [xgrid,ygrid]=meshgrid(linspace(0,10,100),linspace(0,10,100));
    zgrid=ygrid.*sin(xgrid)-xgrid.*cos(ygrid);
    x=10.*rand(100,1);
    y=10.*rand(100,1);
    [z]=interpgrid2xyz(xgrid,ygrid,zgrid,x,y,'yes');

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
