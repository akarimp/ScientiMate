.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-10-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

scientimate.hurricanedpcpt
==========================

.. code:: python

    dPcdt, dPcdthPahr = scientimate.hurricanedpcpt(Pc, dt=6*3600, CalcMethod='backward', dispout='no')

Description
-----------

Calculate hurricane central pressure (Pc) intensity change over time (dPc/dt)

Inputs
------

Pc
    Hurricane central surface pressure in (Pa)
dt=6*3600
    | Time interval between pressure data points in (s)
    | National Hurricane Center reports data every 6 hours 
CalcMethod='backward'
    | Calculation method 
    | 'forward': Calculate hurricane central pressure intensity change over time using forward difference method
    |            If CalcMethod='forward'; then last element is zero
    | 'backward': Calculate hurricane central pressure intensity change over time using backward difference method
    |            If CalcMethod='backward'; then first element is zero
    | 'central': Calculate hurricane central pressure intensity change over time using central difference method
    |            If CalcMethod='central'; then first and last elements are zero
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

dPcdt
    Hurricane central pressure (Pc) intensity change (dPc/dt) in (Pa/s)
dPcdthPahr
    Hurricane central pressure (Pc) intensity change (dPc/dt) in (hPa/hr)

Examples
--------

.. code:: python

    import scientimate as sm

    #Hurricane Katrina centeral pressure (Pa)
    Pc=[100800,100700,100700,100600,100300,100000,99700,99400,98800,98400,98300,98700,\
    97900,96800,95900,95000,94200,94800,94100,93000,90900,90200,90500,91300,\
    92000,92300,92800,94800,96100,97800,98500,99000,99400,99600]

    dPcdt,dPcdthPahr=sm.hurricanedpcpt(Pc,6*3600,'backward','yes')

References
----------

Data

* www.nhc.noaa.gov/data/
* www.nhc.noaa.gov/data/hurdat/hurdat2-format-nencpac.pdf
* coast.noaa.gov/hurricanes
* www.aoml.noaa.gov/hrd/data_sub/re_anal.html

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
