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

curvefit2d
==========

.. code:: MATLAB

    [c, yPredicted] = curvefit2d(x, y, MathExpression, coefIniGuess, fitmethod, dispout)

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
    | Example: y=c(1).*x.^2+c(2).*x+c(3) then MathExpression='c(1).*x.^2+c(2).*x+c(3)'
    | Example: y=c(1).*exp(x)+c(2) then MathExpression='c(1).*exp(x)+c(2)'
    | Example: y=c(1).*sin(x)+c(2) then MathExpression='c(1).*sin(x)+c(2)'
    | Desired coefficients to be found should be as : c(1), c(2),...,c(n)
coefIniGuess
    | Initial guess for desired coefficients to be found
    | coefIniGuess=[guess1,guess2,...,guessn]
    | guess1 is initial guess for c(1),...,guessn is initial guess for c(n) 
fitmethod
    | Fitting method: 
    | 'lsq': curve-fitting using nonlinear least-squares  
    | 'fmin': curve-fitting by minimizing a sum of squared errors
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

c
    Desired coefficients to be found
yPredicted
    Predicted value from fitted curve

Examples
--------

.. code:: MATLAB

    x(:,1)=linspace(0,10,100);
    y(:,1)=0.5.*x.^2+2.*x+10+(-10+(10-(-10))).*rand(100,1);
    [c,yPredicted]=curvefit2d(x,y,'c(1).*x.^2+c(2).*x+c(3)',[1,1,1],'lsq','yes');

    x(:,1)=linspace(0,10,100);
    y(:,1)=5.*x.^2+7+(-200+(200-(-200))).*rand(100,1);
    [c,yPredicted]=curvefit2d(x,y,'c(1).*x.^c(2)+c(3)',[1,2,10],'lsq','yes');

    x(:,1)=linspace(0,10,100);
    y(:,1)=0.5.*exp(x)+100+(-200+(200-(-200))).*rand(100,1);
    [c,yPredicted]=curvefit2d(x,y,'c(1).*exp(x)+c(2)',[1,1],'fmin','yes');

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
