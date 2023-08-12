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

waveorbitalvelocity
===================

.. code:: MATLAB

    [Uw, A, Uwb, Ab] = waveorbitalvelocity(h, H, T, z, kCalcMethod)

Description
-----------

Calculate maximum wave orbital velocity and maximum wave orbital excursion using linear wave theory

Inputs
------
 
h
    Water depth in (m)
H
    | Wave height in (m), H can be approximated and replaced by Hrms for random waves
    | Hrms: root mean square wave height, Hrms=Hm0/sqrt(2)=Hs/sqrt(2) 
    | Hm0: zero moment wave height, Hs: significant wave height  
T
    | Wave period in (s), T can be approximated and replaced by mean wave period for random waves 
    | If peak wave frequency (Tp) is used, calculated values represent peak wave 
z
    Vertical elevation from water surface (z=0 at water surface, z=-h at the bed)
kCalcMethod='beji';
    | Wave number calculation method 
    | 'hunt': Hunt (1979), 'beji': Beji (2013), 'vatankhah': Vatankhah and Aghashariatmadari (2013) 
    | 'goda': Goda (2010), 'exact': calculate exact value 
    | Note: inputs can be as a single value or a 1-D vertical array

Outputs
-------

Uw
    Maximum wave orbital velocity at elevation z in (m/s)
A
    Maximum wave orbital excursion at elevation z in (m)
Uwb
    Wave orbital velocity at bed in (m/s)
Ab
    Wave orbital excursion at bed in (m)

Examples
--------

.. code:: MATLAB

    [Uw,A,Uwb,Ab]=waveorbitalvelocity(1,0.5,3,0.2,'beji');

    [Uw,A,Uwb,Ab]=waveorbitalvelocity([1;1.1],[0.5;0.6],[3;3.1],[0.2;0.3],'exact');

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
