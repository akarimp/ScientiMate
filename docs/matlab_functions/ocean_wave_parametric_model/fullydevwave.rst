.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-09-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

fullydevwave
============

.. code:: MATLAB

    [EhatFullyDev, fphatFullyDev] = fullydevwave(CalcMethod, windvel)

Description
-----------

Calculate a fully developed condition for wind wave growth

Inputs
------

CalcMethod='pierson';
    | Parametric wave model to be used for calculation 
    | 'pierson': Use method by Pierson and Moskowitz (1964)
    | 'spm': Use method by Shore Protection Manual (SPM), U.S. Army Corps of Engineers (1984)
    | 'cem': Use method by Coastal Engineering Manual (CEM),
    |     U.S. Army Corps of Engineers (2015)
windvel=0;
    | Wind velocity in (m/s)
    | Wind velocity should be measured (or represents velocity) at 10 m above surface
    | For 'cem' and 'spm' methods, wind velocity should be converted to duration of sustained wind by using gust factor

Outputs
-------

EhatFullyDev
    Dimensionless wave energy for a fully developed condition, Ehat=g^2*m0/U10^4
fphatFullyDev
    | Dimensionless peak wave frequency for a fully developed condition, fphat=fp*U10/9.81 
    | Note, g=9.81: gravitational acceleration
    |     U10: wind velocity
    |     fp: peak wave frequency
    |     m0: Zero-moment of water surface elevation power spectral density

Examples
--------

.. code:: MATLAB

    [EhatFullyDev,fphatFullyDev]=fullydevwave('pierson')

References
----------

Department of the Army, Waterways Experiment Station, Corps of Engineers, 
and Coastal Engineering Research Center (1984), 
Shore Protection Manual, Washington, 
D.C., vol. 1, 4th ed., 532 pp.

Pierson, W. J., & Moskowitz, L. (1964). 
A proposed spectral form for fully developed wind seas based on the similarity theory of SA Kitaigorodskii. 
Journal of geophysical research, 69(24), 5181-5190.

U.S. Army Corps of Engineers (2015). 
Coastal Engineering Manual. 
Engineer Manual 1110-2-1100, Washington, D.C.: U.S. Army Corps of Engineers.

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
