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

scientimate.similaritymeasure
=============================

.. code:: python

    d = scientimate.similaritymeasure(x, y, CalcMethod='euclidean', dispout='no')

Description
-----------

Measure similarity between two arrays

Inputs
------

x
    First array, its similarity is measured against y array
y
    Second array, its similarity is measured against x array
CalcMethod='euclidean'
    | Similarity  calculation method
    | 'euclidean': Euclidean distance
    | 'manhattan': Manhattan distance
    | 'minkowski': Minkowski distance (power=3)
    | 'cosine': Cosine distance
    | 'pearson': Pearson's correlation coefficient
    | 'spearman': spearman's correlation coefficient
    | 'norm': Absolute difference of norm
    | 'covariance': Covariance
    | 'inv_covariance': Euclidean distance of inverse covariance
    | 'histogram': Mean of absolute difference of histogram
    | t-test': Two-sample t-test statistic
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

d
    Arrays similarity measure

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np

    x=[0,2,4,6]
    y=[2,3,5,7]
    d=sm.similaritymeasure(x,y,'euclidean','yes')

    x=[1,2,3,4,5,6,7,8,9,10]
    y=[1.1,1.98,3.3,4.2,4.8,5.95,7.5,7.7,8.99,10.5]
    d=sm.similaritymeasure(x,y,'pearson','yes')

References
----------
Kianimajd, A., Ruano, M. G., Carvalho, P., Henriques, J., Rocha, T., Paredes, S., & Ruano, A. E. (2017).
Comparison of different methods of measuring similarity in physiologic time series.
IFAC-PapersOnLine, 50(1), 11005-11010.

* https://en.wikipedia.org/wiki/Similarity_measure
* https://en.wikipedia.org/wiki/Goodness_of_fit
* https://dataaspirant.com/2015/04/11/five-most-popular-similarity-measures-implementation-in-python/
* https://towardsdatascience.com/similarity-measures-e3dbd4e58660
* https://www.mathworks.com/matlabcentral/answers/377944-how-to-calculate-a-percentage-of-similarity-between-two-arrays
* https://en.wikipedia.org/wiki/Template_matching
* https://www.mathworks.com/help/images/ref/normxcorr2.html

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
