.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2020-08-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

scientimate.movingwindow
========================

.. code:: python

    moving_statistics = scientimate.movingwindow(x, WindowSize=3, StatisticalMethod='mean', dispout='no')

Description
-----------

Calculate statistics of moving window through 1-d x data

Inputs
------

x
    Input data
WindowSize=3
    | Window size (number of adjacent elements) that is used for moving window, should be equal or larger than 3
    | Window size should be an odd integer
StatisticalMethod='mean'
    | Statistical value of moving window to be reported:
    | 'mean': Moving mean
    | 'std': Moving standard deviation
    | 'min': Moving minimum
    | 'max': Moving maximum
    | 'sum': Moving sum
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

moving_statistics
    Statistical value of moving window

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np

    fs=128
    t=np.linspace(0,9.5,10*fs)
    x=np.sin(2*np.pi*0.3*t)+0.1*np.sin(2*np.pi*4*t)
    moving_statistics=sm.movingwindow(x,37,'mean','yes')

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
