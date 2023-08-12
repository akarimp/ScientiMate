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

scientimate.curvefit2d
======================

.. code:: python

    c, yPredicted = scientimate.curvefit2d(x, y, MathExpression, coefIniGuess, fitmethod='fmin', dispout='no')

Description
-----------

Fit curve to 2 dimensinal input dataset

Inputs
------

x
    x data
y
    y data
MathExpression
    | Right hand side of the relation y=f(x) as 'f(x)'
    | Example: y=c[0]*x**2+c[1]*x+c[2] then MathExpression='c[0]*x**2+c[1]*x+c[2]'
    | Example: y=c[0]*exp(x)+c[1] then MathExpression='c[0]*np.exp(x)+c[1]'
    | Example: y=c[0]*sin(x)+c[1] then MathExpression='c[0]*np.sin(x)+c[1]'
    | Desired coefficients to be found should be as : c[0], c[1],...,c[n]
coefIniGuess
    | Initial guess for desired coefficients to be found
    | coefIniGuess=[guess0,guess1,...,guessn]
    | guess0 is initial guess for c[0],...,guessn is initial guess for c[n] 
fitmethod
    | Fitting method: 
    | 'lsq': curve-fitting using nonlinear least-squares  
    | 'fmin': curve-fitting by minimizing a sum of squared errors
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

c
    Desired coefficients to be found
yPredicted
    Predicted value from fitted curve

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np

    x=np.linspace(0,10,100)
    y=0.5*x**2+2*x+10+(-10+(10-(-10)))*np.random.rand(100)
    c,yPredicted=sm.curvefit2d(x,y,'c[0]*x**2+c[1]*x+c[2]',[1,1,1],'fmin','yes')

    x=np.linspace(0,10,100)
    y=5*x**2+7+(-200+(200-(-200)))*np.random.rand(100)
    c,yPredicted=sm.curvefit2d(x,y,'c[0]*x**c[1]+c[2]',[1,2,10],'lsq','yes')

    x=np.linspace(0,10,100)
    y=0.5*np.exp(x)+100+(-200+(200-(-200)))*np.random.rand(100)
    c,yPredicted=sm.curvefit2d(x,y,'c[0]*np.exp(x)+c[1]',[1,1],'fmin','yes')

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
