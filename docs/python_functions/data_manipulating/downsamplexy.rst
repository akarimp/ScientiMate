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

scientimate.downsamplexy
========================

.. code:: python

    x_ds, y_ds = scientimate.downsamplexy(x, y, RetainRatio)

Description
-----------

Downsample x and y data and retain given ratio

Inputs
------

x
    x data
y
    y data
RetainRatio=0.5
    | Define percentage of data to retain, value between 0 and 1
    | Example: RetainRatio=0.8 means 80% of data are retained

Outputs
-------

x_ds
    Downsample x data
y_ds
    Downsample y data

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np
    from numpy import random

    rng = np.random.default_rng()
    x=10*rng.random((1000,1))
    y=x**2
    x_ds, y_ds=sm.downsamplexy(x, y, 0.7)

    rng = np.random.default_rng()
    xgrid=(-90-(-91))*rng.random((1000,500))+(-91)
    ygrid=xgrid**2
    x_ds, y_ds=sm.downsamplexy(xgrid, ygrid, 0.3)

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
