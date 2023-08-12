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

fitgoodness
===========

.. code:: MATLAB

    [r, R2, RMSE, MAE, SI, NSE, d, Bias, NMBias, RE] = fitgoodness(x, y, dispout)

Description
-----------

Calculate goodness of fit parameters

Inputs
------

x
    Dataset with true (exact or expected) values, such as theoretical values
y
    | Dataset that needs to be evaluated, such as model results or estimated values
    | Accuracy of y dataset is evaluated against x dataset
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

r
    Pearson correlation coefficient
R2
    Coefficient of determination
RMSE
    Root mean square error
MAE
    Mean absolute error
SI
    Scatter index
NSE
    Nash Sutcliffe efficiency coefficient
d
    Index of agreement
Bias
    Bias
NMBias
    Normalized mean bias
RE
    Relative error

Examples
--------

.. code:: MATLAB

    x(:,1)=(-0.1+(0.1-(-0.1))).*randn(1024*2,1);
    y(:,1)=x+(-0.01+(0.01-(-0.01))).*randn(1024*2,1);
    [r,R2,RMSE,MAE,SI,NSE,d,Bias,NMBias,RE]=fitgoodness(x,y,'yes');

    x=[1;2;3;4;5;6;7;8;9;10];
    y=[1.1;1.98;3.3;4.2;4.8;5.95;7.5;7.7;8.99;10.5];
    [r,R2,RMSE,MAE,SI,NSE,d,Bias,NMBias,RE]=fitgoodness(x,y,'yes');

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
