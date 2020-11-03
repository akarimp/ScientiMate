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

scientimate.findknn
===================

.. code:: python

    indxknn, distknn = scientimate.findknn(x, y, xpoint, ypoint, numofneighbors=1, distCalcMethod='1d', dispout='no')

Description
-----------

Find k-nearest neighbors using Euclidean distance

Inputs
------

x
    Coordinate of data in x direction
y
    Coordinate of data in y direction
xpoint
    Coordinate (in x direction) of point that nearest point to that is disired to be found 
ypoint
    Coordinate (in y direction) of point that nearest point to that is disired to be found 
numofneighbors=1
    Number of nearest neighbors to (xpoint,ypoint) that are desired to be found
distCalcMethod='1d'
    | Distance calculation method 
    | '1d': use 1d array
    | 'pdist2': Use 2d distance function
    | 'vector': Use vectorized distance 
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

indxknn
    | Index of nearest neighbors points
    | returns M*N array where M=length(xpoint) and N=numofneighbors
    | nth row associated with nth point in (xpoint,ypoint)
distknn
    | Distance of nearest neighbors points
    | returns M*N array where M=length(xpoint) and N=numofneighbors
    | mth row associated with mth point in (xpoint,ypoint)
    | nth column associated with nth nearest neighbors

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np

    x=10*np.random.rand(100)
    y=10*np.random.rand(100)
    xpoint=np.mean(x)
    ypoint=np.mean(y)
    indxknn,distknn=sm.findknn(x,y,xpoint,ypoint,1,'1d','yes')

    x=10*np.random.rand(100)
    y=10*np.random.rand(100)
    xpoint=[2.5,5,7.5]
    ypoint=[3,6,9]
    indxknn,distknn=sm.findknn(x,y,xpoint,ypoint,10,'1d','yes')

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
