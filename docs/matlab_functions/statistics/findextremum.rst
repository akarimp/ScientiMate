.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2018-05-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

findextremum
============

.. code:: MATLAB

    [xmin, ymin, xmax, ymax] = findextremum(x, y, winlen, dispout)

Description
-----------

Find local extremum (minimum and maximum) in data

Inputs
------

x
    x data
y
    y data
winlen=3;
    | Window length, defines a number of points in sliding window used for defining maximum and minimum
    | Example: winlen=5 means two points on each side of each data point is used in calculation
    | Using a larger value for winlen makes it less sensitive 
    | winlen should be an odd number equal or larger than 3
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

xmin
    x of minmum points
ymin
    y of minmum points
xmax
    x of maximum points
ymax
    y of maximum points

Examples
--------

.. code:: MATLAB

    x(:,1)=linspace(0,30,1000);
    y=2.*exp(-0.1*2*pi/5.*x).*sin(sqrt(1-0.1^2)*2*pi/5.*x);
    [xmin,ymin,xmax,ymax]=findextremum(x,y,3,'yes');

    x(:,1)=linspace(0,50,1000);
    y(:,1)=sin(x)+0.1.*rand(1000,1);
    [xmin,ymin,xmax,ymax]=findextremum(x,y,15,'yes');

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
