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

replacemissing2d
================

.. code:: MATLAB

    [xReplaced, NaN_Indx] = replacemissing2d(x, what2replace, gridsize_x, gridsize_y, interpMethod, dispout)

Description
-----------

Replace missing data points in 2d array

Inputs
------

x
    Input data
what2replace='both';
    | What needs to be replaced
    | 'NaN': replacing NaN data points
    | 'Inf': replacing Inf data points
    | 'both': replacing NaN and Inf data points
    | Number: replacing data points equal to Number
gridsize_x=1;
    | Grid size (distance between grid points) in x direction
    | Leave gridsize_x=1 if you do not have it
gridsize_y=1;
    | Grid size (distance between grid points) in y direction
    | Leave gridsize_y=1 if you do not have it
interpMethod='nearest';
    | Interpolation method
    | 'linear': Use default or 'linear' method to interpolate
    | 'nearest': Use nearest neighbor method to interpolate
    | 'knn': Use nearest neighbor method to interpolate (Use 'knn' for large array)
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

xReplaced
    Replaced data
NaN_Indx
    Logical index of replaced points

Examples
--------

.. code:: MATLAB

    x=[[1,0,3];[2,5,NaN];[3,NaN,1];[5,7,2]];
    [xReplaced, NaN_Indx] = replacemissing2d(x, 'NaN', 1, 1, 'nearest', 'yes');

    xgrid=randn(100,50);
    xgrid(randi(100,20,1),randi(50,20,1))=NaN;
    [xReplaced, NaN_Indx] = replacemissing2d(xgrid, 'NaN', 1, 1, 'knn', 'yes');

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
