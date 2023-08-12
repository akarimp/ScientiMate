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

surfaceroughness
================

.. code:: MATLAB

    [ustar, z0, d] = surfaceroughness(z, u, delta, dispout)

Description
-----------

Calculate shear velocity and surface roughness from a given velocity profile using Karimpour et al. (2012) method

Inputs
------

z
    Distance from a surface (elevation, height) in (m)
u
    Velocity at z in (m/s)
delta=max(z);
    Boundary layer height in (m)
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

z0
    Surface roughness in (m)
ustar
    Shear velocity (u*) in (m/s)
d
    | Zero plane displacement distance in (m)
    | Note: Above values are for a logarithmic velocity profile as:
    |     u=(u*/K)*ln((z-d)/z0)

Examples
--------

.. code:: MATLAB

    z(:,1)=[0.1:0.05:1];
    u(:,1)=2/0.4*log((z-0.003)/0.002);
    [ustar,z0,d]=surfaceroughness(z,u,max(z),'yes');

References
----------

Karimpour, A., Kaye, N. B., & Baratian-Ghorghi, Z. (2012). 
Modeling the neutrally stable atmospheric boundary layer for laboratory scale studies of the built environment. 
Building and Environment, 49, 203-211.

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
