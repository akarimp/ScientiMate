.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2020-08-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

windspectrum
============

.. code:: MATLAB

    [S, f] = windspectrum(U10, Sigma, L, CalcMethod, fmin, fmax, dispout)

Description
-----------

Calculate wind spectrum with 513 frequencies

Inputs
------

U10
    Wind velocity at 10 m above surface in (m/s)
Sigma
    Wind velocity standard deviation 
        IEC 61400-1:
            | Sigma_u=Sigma
            | Sigma_v=0.8*Sigma
            | Sigma_w=0.5*Sigma
L
    Wind velocity integral length scale
        Suggestion to calculate L (wind velocity integral length scale):
            For Kaimal (1972) suggested by DNV-RP-C205 (2010):
                L=300*(z/300)^(0.46+0.074*log(z0))
            For Kaimal (1972) suggested by DNV-RP-C205 (2010) based on IEC 61400-1:
                | L=3.33*z for 0<z<=60 
                | L=200 for z>=60
            For Kaimal (1972) suggested by IEC 61400-1:
                | Lambda=0.7*z for for z<=60
                | Lambda=42 m for for z>60
                | Lu=8.1*Lambda
                | Lv=2.7*Lambda
                | Lw=0.66*Lambda
            For Davenport (1961):
                Lu=1200 m
            For Harris (1971):
                The same as the one for Kaimal (1972) suggested by DNV-RP-C205 (2010)
            z
                Height from surface in (m)
            z0
                Surface roughness
CalcMethod='kaimal_dnv';
    | Wind spectrum calculation method 
    | 'dnv': DNV-RP-C205 (2010)
    | 'kaimal_dnv': Kaimal (1972) suggested by DNV-RP-C205 (2010)
    | 'kaimal_iec': Kaimal (1972) suggested by IEC 61400-1
    | 'davenport': Davenport (1961)
    | 'harris': Harris (1971)
fmin=0.002;
    Lower boundary of spectrum in (Hz)
fmax=0.3;
    | Upper boundary of spectrum in (Hz)
    | Kaimal is valid in range from 0.002 Hz to 0.3 Hz (from 7.5 to 1000 cycles/hour)
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

S
    Wind spectrum in ((m/s)^2/Hz)
f
    Frequency in (Hz)

Examples
--------

.. code:: MATLAB

    U10=15;
    Sigma=1;
    z=30;
    L=3.33*z;
    [S, f] = windspectrum(U10, Sigma, L, 'kaimal_dnv', 0, 1, 'yes');

    U10=15;
    Sigma=1;
    z=30;
    Lambda=0.7*z;
    L==8.1*Lambda;
    [S, f] = windspectrum(U10, Sigma, L, 'kaimal_iec', 0, 1, 'yes');

References
----------

Bhattacharya, S. (2019).
Design of foundations for offshore wind turbines.
Wiley.

Bec, J. (2010).
Influence of wind spectrum formula choice on footbridge response.
In 5th international symposium on computational wind engineering (pp. 23-27).

Branlard, E. (2010).
Generation of time series from a spectrum.
Technical University Denmark. National Laboratory for Sustainable Energy.

Davenport, A. G. (1961). 
The spectrum of horizontal gustiness near the ground in high winds. 
Quarterly Journal of the Royal Meteorological Society, 87(372), 194-211.

Harris, R. I. The Nature of Wind, 
Proc. of the ModernDesign of Wind Sensitive Structures, Construction,
Industry Research and Information Association, 1971, London, U. K

Kaimal, J. C., Wyngaard, J. C. J., Izumi, Y., & Coté, O. R. (1972). 
Spectral characteristics of surface‐layer turbulence. 
Quarterly Journal of the Royal Meteorological Society, 98(417), 563-589.

Rose, S., & Apt, J. (2012). 
Generating wind time series as a hybrid of measured and simulated data. 
Wind Energy, 15(5), 699-715.

Udoh, I. E., & Zou, J. (2018). 
Wind spectral characteristics on strength design of floating offshore wind turbines. 
Ocean Systems Engineering, 8(3), 281-312.

VERITAS, D. N. (2010). ENVIRONMENTAL CONDITIONS AND ENVIRONMENTAL LOADS.

https://www.mathworks.com/help/aeroblks/wind.html

https://www.mathworks.com/help/aeroblks/vonkarmanwindturbulencemodelcontinuous.html

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
