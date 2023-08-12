.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-12-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

gridgenerator
=============

.. code:: MATLAB

    [xgrid, ygrid] = gridgenerator(xmin, xmax, ymin, ymax, gridsize, gridsizetype, dispout)

Description
-----------

Generate 2d x-y grid

Inputs
------

xmin
    Minimum x of the domain to be generated
xmax
    Maximum x of the domain to be generated
ymin
    Minimum y of the domain to be generated
ymax
    Maximum y of the domain to be generated
gridsize=100;
    | Grid size in x (longitude) and y (latitude) directions to interpolate elevation data on them
    |     if gridsizetype='length' then gridsize is a distance between grid points
    |     if gridsizetype='points' then gridsize is number of grid points in each direction
gridsizetype='points';
    | Grid size type 
    |     'number': gridsize is considered as number of grid points in each direction
    |     'length': gridsize is considered as length between grid points
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

xgrid
    x of the defined mesh
ygrid
    y of the defined mesh

Examples
--------

.. code:: MATLAB

    [xgrid,ygrid]=gridgenerator(0,100,0,100,100,'points','yes');
    [xgrid,ygrid]=gridgenerator(0,100,0,100,20,'length','yes');

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
