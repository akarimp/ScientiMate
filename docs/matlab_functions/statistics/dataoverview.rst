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

dataoverview
============

.. code:: MATLAB

    dataoverview(x, y, z)

Description
-----------

Display an overview of the input data

Inputs
------

x
    x data
y=[];
    y data (Optional)
z=[];
    z data (Optional)

Outputs
-------


Examples
--------

.. code:: MATLAB

    x(:,1)=(-0.1+(0.1-(-0.1))).*randn(1024*2,1);
    dataoverview(x);

    x(:,1)=(-0.1+(0.1-(-0.1))).*randn(1024*2,1);
    y(:,1)=x+(-0.01+(0.01-(-0.01))).*randn(1024*2,1);
    dataoverview(x,y);

    x(:,1)=(-0.1+(0.1-(-0.1))).*randn(1024*2,1);
    y(:,1)=(-0.2+(0.2-(-0.2))).*randn(1024*2,1);
    z=x.^2+y.^2+0.1+(-0.1+(0.1-(-0.1))).*rand(1024*2,1);
    dataoverview(x,y,z);

References
----------

* https://www.mathworks.com/help/matlab/ref/fprintf.html

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
