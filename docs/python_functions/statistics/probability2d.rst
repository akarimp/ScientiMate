.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-06-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

scientimate.probability2d
=========================

.. code:: python

    fdensityxy, fdensityx, fdensityy, fdensitycumulativex, fdensitycumulativey, bincenterx, bincentery, xmean, ymean, xstd, ystd = scientimate.probability2d(x, y, binedgex=None, binedgey=None, dispout='no')

Description
-----------

Calculate 2D (joint) probability density distribution for two given datasets

Inputs
------

x
    First input dataset 
y
    Second input dataset 
binedgex
    | Bin edges for x data  
    | length(binedgex)=number of bin +1   
    | If there are N bins in histogram/distribution, then values in binedgex are as:   
    | 1st value: left edge of 1st bin, 2nd value: left edge of 2nd bin, ...   
    | (N)th value: left edge of last bin, (N+1)th value: right edge of last bin   
binedgey
    | Bin edge for y data  
    | length(binedgey)=number of bin +1   
    | If there are N bins in histogram/distribution, then values in binedgey are as:   
    | 1st value: left edge of 1st bin, 2nd value: left edge of 2nd bin, ...   
    | (N)th value: left edge of last bin, (N+1)th value: right edge of last bin   
dispout='no'
    | Define to display outputs or not
    | 'no': not display 
    | 'bar': bar plot
    | 'imagesc': 2 dimensional plot using imagesc or imshow
    | 'pcolor': 2 dimensional plot using pcolor
    | 'contour': 2 dimensional contour plot, number of contour=32
    | 'contourf': 2 dimensional filled contour plot, number of contour=32
    | 'surface': 3 dimensional surface plot 
    | 'binscatter': 2 dimensional histogram plot using imagesc or imshow

Outputs
-------

fdensityxy
    2D (joint) probability density distribution for x-y data
fdensityx
    Probability density distribution for x data
fdensityy
    Probability density distribution for y data
fdensitycumulativex
    Cumulative probability density distribution for x data
fdensitycumulativey
    Cumulative probability density distribution for y data
bincenterx
    Bin center for x data
bincentery
    Bin center for y data
xmean
    Mean value of x data
ymean
    Mean value of y data
xstd
    Standard deviation of x data
ystd
    Standard deviation of y data

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np

    x=(-0.1+(0.1-(-0.1)))*np.random.randn(1024*2)
    y=(-0.2+(0.2-(-0.2)))*np.random.randn(1024*2)
    binedgex=np.linspace(min(x),max(x),11)
    binedgey=np.linspace(min(y),max(y),11)
    fdensityxy,fdensityx,fdensityy,fdensitycumulativex,fdensitycumulativey,bincenterx,bincentery,xmean,ymean,xstd,ystd=sm.probability2d(x,y,binedgex,binedgey,'surface')

    x=np.random.randn(100000)
    y=1.5*x+np.random.randn(100000)
    binedgex=np.linspace(min(x),max(x),101)
    binedgey=np.linspace(min(y),max(y),101)
    sm.probability2d(x,y,binedgex,binedgey,'binscatter')

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
