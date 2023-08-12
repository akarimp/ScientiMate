.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-01-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

scientimate.linearwavesuperposition
===================================

.. code:: python

    Eta, t, Etaij = scientimate.linearwavesuperposition(a, T, Phi, fs=32, duration=10, dispout='no')

Description
-----------

Superposition linear waves

Inputs
------

a
    Wave amplitude in (m)
T
    Wave mean period in (s)
Phi
    Phase (radian)
fs=32
    Sample generation frequency (Hz), number of data points in one second
duration=10
    Duration time that data will be generated in (s)
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

Eta
    Water Surface Level Time Series in (m)
t
    Time in (s)
Etaij
    Separated Water Surface Level Time Series in (m)

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np

    Eta,t,Etaij=sm.linearwavesuperposition([0.1,0.2,0.3,0.4],[1,1.5,2,2.5],[np.pi/2,np.pi/4,np.pi/16,np.pi/32],32,10,'yes')

    Eta,t,Etaij=sm.linearwavesuperposition(np.array([0.1,0.2,0.3,0.4]),np.array([1,1.5,2,2.5]),np.array([np.pi/2,np.pi/4,np.pi/16,np.pi/32]),32,10,'yes')

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
