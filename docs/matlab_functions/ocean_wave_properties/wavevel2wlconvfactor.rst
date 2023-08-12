.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-04-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

wavevel2wlconvfactor
====================

.. code:: MATLAB

    [Kuv, f] = wavevel2wlconvfactor(f, h, heightfrombed, kCalcMethod, dispout)

Description
-----------

Calculate a water particle horizontal orbital velocity to the water surface elevation conversion factor

Inputs
------

f
    Frequency in (Hz) 
h
    Water depth in (m)
heightfrombed
    Height from bed that Kp calculated at in (m)
kCalcMethod='beji';
    | Wave number calculation method 
    | 'hunt': Hunt (1979), 'beji': Beji (2013), 'vatankhah': Vatankhah and Aghashariatmadari (2013) 
    | 'goda': Goda (2010), 'exact': calculate exact value 
dispout='no';
    | Define to display outputs or not ('yes': display, 'no': not display)
    | Note: inputs can be as a single value or a 1-D vertical array

Outputs
-------

Kuv
    Horizontal orbital velocity to the water surface elevation conversion factor
f
    Frequency (Hz)

Examples
--------

.. code:: MATLAB

    [Kuv,f]=wavevel2wlconvfactor(0.2,1,0.2,'beji','no');

    [Kuv,f]=wavevel2wlconvfactor([0.2;0.25],[0.5;0.6],0.2,'exact','yes');

References
----------

Beji, S. (2013). 
Improved explicit approximation of linear dispersion relationship for gravity waves. 
Coastal Engineering, 73, 11-12.

Goda, Y. (2010). 
Random seas and design of maritime structures. 
World scientific.

Hunt, J. N. (1979). 
Direct solution of wave dispersion equation. 
Journal of the Waterway Port Coastal and Ocean Division, 105(4), 457-459.

Vatankhah, A. R., & Aghashariatmadari, Z. (2013). 
Improved explicit approximation of linear dispersion relationship for gravity waves: A discussion. 
Coastal engineering, 78, 21-22.

Wiberg, P. L., & Sherwood, C. R. (2008). 
Calculating wave-generated bottom orbital velocities from surface-wave parameters. 
Computers & Geosciences, 34(10), 1243-1262.

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
