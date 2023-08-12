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

scientimate.dataoverview
========================

.. code:: python

    scientimate.dataoverview(x, y=[], z=[])

Description
-----------

Display an overview of the input data

Inputs
------

x
    x data
y=[]
    y data (Optional)
z=[]
    z data (Optional)

Outputs
-------


Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np

    x=(-0.1+(0.1-(-0.1)))*np.random.randn(1024*2)
    sm.dataoverview(x)

    x=(-0.1+(0.1-(-0.1)))*np.random.randn(1024*2)
    y=x+(-0.01+(0.01-(-0.01)))*np.random.randn(1024*2)
    sm.dataoverview(x,y)

    x=(-0.1+(0.1-(-0.1)))*np.random.randn(1024*2)
    y=(-0.2+(0.2-(-0.2)))*np.random.randn(1024*2)
    z=x**2+y**2+0.1+(-0.1+(0.1-(-0.1)))*np.random.rand(1024*2)
    sm.dataoverview(x,y,z)

References
----------

* https://docs.python.org/3/library/string.html

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
