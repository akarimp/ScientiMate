.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-02-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

scientimate.wavedispersionds
============================

.. code:: python

    k, L, C, Cg = scientimate.wavedispersionds(h, T, Uc=0)

Description
-----------

| Solve water wave dispersion relation with presence of current (Doppler shift)
| Calculate wave number (k), wave length (L), wave celereity (C), and wave group velocity (Cg) using linear wave theory

Inputs
------

h
    Water depth in (m)
T
    | Wave period in (s) 
    | If peak wave frequency (Tp) is used, calculated values represent peak wave 
Uc=0
    | Current velocity in (m/s), for Doppler shift
    | Uc should have a same size as h
    | Note: inputs can be as a single value or a 1-D vertical array

Outputs
-------

k
    Wave number in (radian/m)
L
    Wave length in (m)
C
    Wave celerity in (m/s)
Cg
    Wave group celerity in (m/s)

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np

    k,L,C,Cg=sm.wavedispersionds(1,3,1)

    k,L,C,Cg=sm.wavedispersionds([1,1.1],[3,3.1],[1,1])

    k,L,C,Cg=sm.wavedispersionds(np.array([1,1.1]),np.array([3,3.1]),np.array([1,1]))

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
