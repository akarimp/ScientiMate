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

scientimate.curvefit3d
======================

.. code:: python

    c, zPredicted = scientimate.curvefit3d(x, y, z, MathExpression, coefIniGuess, fitmethod='fmin', dispout='no')

Description
-----------

Fit curve to 3 dimensinal input dataset

Inputs
------

x
    x data
y
    y data
z
    z data
MathExpression
    | Right hand side of the relation z=f(x,y) as 'f(x,y)'
    | Example: z=c[0]*x**2+c[1]*y+c[2] then MathExpression='c[0]*x**2+c[1]*y+c[2]'
    | Example: z=c[0]*exp(x)+c[1]*sin(y) then MathExpression='c[0]*np.exp(x)+c[1]*np.sin(y)'
    | Desired coefficients to be found should be as : c[0], c[1],...,c[n]
coefIniGuess
    | Initial guess for desired coefficients to be found
    | coefIniGuess=[guess0,guess1,...,guessn]
    | guess0 is initial guess for c[0],...,guessn is initial guess for c(n) 
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
zPredicted
    Predicted value from fitted curve

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np

    xx,yy=np.meshgrid(np.linspace(0,10,100),np.linspace(0,10,100))
    zz=xx**2+yy**2+10+(-10+(10-(-10)))*np.random.rand(100,100)
    x=xx.reshape(-1) #Flatten array
    y=yy.reshape(-1) #Flatten array
    z=zz.reshape(-1) #Flatten array
    c,zPredicted=sm.curvefit3d(x,y,z,'c[0]*x**2+c[1]*y**2',[1,1],'fmin','yes')

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
