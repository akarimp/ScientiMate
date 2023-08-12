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

downsamplexyz
=============

.. code:: MATLAB

    function [x_ds, y_ds, z_ds] = downsamplexyz(x, y, z, RetainRatio)

Description
-----------

Downsample x, y,  and z data and retain given ratio

Inputs
------

x
    x data
y
    y data
z
    z data
RetainRatio=0.5;
    | Define percentage of data to retain, value between 0 and 1
    | Example: RetainRatio=0.8; means 80% of data are retained

Outputs
-------

x_ds
    Downsample x data
y_ds
    Downsample y data
z_ds
    Downsample z data

Examples
--------

.. code:: MATLAB

    x(:,1)=10.*rand(1000,1);
    y(:,1)=10.*rand(1000,1);
    z=x.^2+y.^2;
    [x_ds, y_ds, z_ds]=downsamplexyz(x, y, z, 0.7);

    x(:,1)=(-90-(-91)).*rand(1000,1)+(-91);
    y(:,1)=(31-(30)).*rand(1000,1)+(30);
    [xgrid,ygrid]=meshgrid(linspace(min(x),max(x),1000),linspace(min(y),max(y),500));
    zgrid=xgrid.^2+ygrid.^2;
    [x_ds, y_ds, z_ds]=downsamplexyz(xgrid, ygrid, zgrid, 0.3);

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
