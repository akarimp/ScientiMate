.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2020-01-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

smoothsignal
============

.. code:: MATLAB

    [xSmoothed] = smoothsignal(x, WindowSize, method, dispout)

Description
-----------

Smooth input data using a window function

Inputs
------

x
    Input data
WindowSize=33;
    Window size (number of adjacent elements)that is used for smoothing, should be equal or larger than 3
method='moveavg';
    | Smoothing method
    | 'moveavg': moving average
    | 'lowpass': Low-Pass filter
    | 'savgol': Savitzky-Golay filter
    | 'butter': Butterworth filter
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

xSmoothed
    Smoothed data

Examples
--------

.. code:: MATLAB

    x(:,1)=sin(linspace(0,6*pi,200))+rand(1,200);
    [xSmoothed]=smoothsignal(x,33,'moveavg','yes');

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
