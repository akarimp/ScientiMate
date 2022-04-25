.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2018-03-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

scientimate.stormsurge1d
========================

.. code:: python

    Eta, x, maxsurgeheight, L = scientimate.stormsurge1d(h0, U10, m=0, dispout='no')

Description
-----------

Calculate one dimensional storm surge using Dean Dalrymple (1991) method

Inputs
------

h0
    Deep-water depth in (m), 
U10
    10-min averaged wind velocity at 10 meter above a surface in (m/s)
m=0
    | Bed slope
    | Note: m=h0/L, where L is a length of the continental shelf in (m)
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

Eta
    Surge height along a x axis in (m)
x
    | Points on x axis in (m) 
    | x=0 located at a far-end of the water boundary
    | x=L located at a coastline
maxsurgeheight
    Maximum surge height at a coastline (x=L) in (m)
L
    | Length of the continental shelf in (m)
    | L=h0/m if m is not zero
    | L=50000 if m is zero (L=50000 meter for the hurricane Katrina)

Examples
--------

.. code:: python

    import scientimate as sm

    h0=42 #Example: h0=42 m for the hurricane Katrina
    U10=40 #Example: U10=40 m/s for the hurricane Katrina 
    m=0.00084
    Eta,x,maxsurgeheight,L=sm.stormsurge1d(h0,U10,m,'yes')

References
----------

Dean, R. G., & Dalrymple, R. A. (1991). 
Water wave mechanics for engineers and scientists (Vol. 2). 
World Scientific Publishing Company.

Wu, J. (1982). 
Wind‐stress coefficients over sea surface from breeze to hurricane. 
Journal of Geophysical Research: Oceans, 87(C12), 9704-9706.

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
