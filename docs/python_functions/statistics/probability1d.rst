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

scientimate.probability1d
=========================

.. code:: python

    fdensity, fdensitycumulative, bincenter, xmean, xstd = scientimate.probability1d(x, binedge=None, dispout='no')

Description
-----------

Calculate 1D probability density distribution for a given dataset

Inputs
------

x
    Input data 
binedge
    | Bin edges  
    | length(binedge)=number of bin +1   
    | If there are N bins in histogram/distribution, then values in binedge are as:   
    | 1st value: left edge of 1st bin, 2nd value: left edge of 2nd bin, ...   
    | (N)th value: left edge of last bin, (N+1)th value: right edge of last bin   
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

fdensity
    Probability density distribution
fdensitycumulative
    Cumulative probability density distribution
bincenter
    Bin center
xmean
    Mean value of input data
xstd
    Standard deviation of input data

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np

    x=(-0.1+(0.1-(-0.1)))*np.random.randn(1024*2)
    binedge=np.linspace(min(x),max(x),11)
    fdensity,fdensitycumulative,bincenter,xmean,xstd=sm.probability1d(x,binedge,'yes')

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
