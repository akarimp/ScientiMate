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

scientimate.distancecart
========================

.. code:: python

    distxy, theta = scientimate.distancecart(x1, y1, x2, y2, CalcMethod='1d', dispout='no')

Description
-----------

Calculate distance from (x1,y1) to (x2,y2) on cartesian coordinate

Inputs
------

x1
    x of start point (first point)
y1
    y of start point (first point)
x2
    x of end point (last point) 
y2
    y of end point (last point) 
CalcMethod='1d'
    | Distance calculation method 
    | '1d': use 1d array
    | 'pdist2': Use 2d distance function
    | 'vector': Use vectorized distance 
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

distxy
    | Distance from (x1,y1) to (x2,y2)
    | returns M*N array where M=length(x1) and N=length(x2)
    | mth row associated with mth point in (x,y)
    | nth column is associated with nth point in (x2,y2)
theta
    | Angle from start point to end point in (Degree)
    | returns M*N array where M=length(x1) and N=length(x2)
    | mth row associated with mth point in (x,y)
    | nth column is associated with nth point in (x2,y2)

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np

    x1=10*np.random.rand(100)
    y1=10*np.random.rand(100)
    x2=[2.5,5,7.5]
    y2=[3,6,9]
    distxy,theta=sm.distancecart(x1,y1,x2,y2,'1d','yes')

    x1=10*np.random.rand(100)
    y1=10*np.random.rand(100)
    x2=100*np.random.rand(10)
    y2=100*np.random.rand(10)
    distxy,theta=sm.distancecart(x1,y1,x2,y2,'pdist2','yes')

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
