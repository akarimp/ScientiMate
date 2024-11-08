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

scientimate.wavedim2dimless
===========================

.. code:: python

    dimlessValue = scientimate.wavedim2dimless(windvel, dimValue, ValueType='waveheight')

Description
-----------

Calculate dimensionless numbers from dimensional numbers

Inputs
------

windvel
    Wind velocity in (m/s)
dimValue
    Dimensional value to be converted to dimensionless value
ValueType='waveheight'
    | Type of the dimensional value 
    | 'fetch': Wind fetch in (m)
    | 'depth': Water depth in (m)
    | 'frequency': Wave frequency in (Hz)
    | 'energy': Zero-moment of surface elevation power spectral density in (m^2)
    | 'waveheight': Wave height in (m)
    | 'period': Wave period in (s)
    | 'wavelength': Wave length in (m)
    | 'wavenumber': Wave number in (1/m) or (Radian/m)
    | 'time': Time in (s)
    | 'length': Length in (m)

Outputs
-------

dimlessValue
    | Dimensionless value calculated from dimensional value
    | For ValueType='fecth': dimlessValue is dimensionless wind fetch (Fetchhat): Fetchhat=g*Fetch/U10^2
    | For ValueType='depth': dimlessValue is dimensionless depth (hhat): hhat=g*h/U10^2
    | For ValueType='frequency': dimlessValue is dimensionless wave frequency (fhat): fhat=f*U10/9.81
    | For ValueType='energy': dimlessValue is dimensionless wave energy (Ehat): Ehat=g^2*m0/U10^4
    | For ValueType='waveheight': dimlessValue is dimensionless wave height (Hhat): Hhat=g*H/U10^2
    | For ValueType='period': dimlessValue is dimensionless wave period (That): That=g*T/U10
    | For ValueType='wavelength': dimlessValue is dimensionless wave length (Lhat): Lhat=g*L/U10^2
    | For ValueType='wavenumber': dimlessValue is dimensionless wave number (khat): khat=k*U10^2/g
    | For ValueType='time': dimlessValue is dimensionless time (that): that=g*t/U10
    | For ValueType='length': dimlessValue is dimensionless time (that): that=g*t/U10
    | Note, g=9.81: gravitational acceleration
    |     U10: wind velocity
    |     Fetch: Wind fetch in (m)
    |     h: Water depth in (m)
    |     f: Wave frequency in (Hz)
    |     m0: Zero-moment of water surface elevation power spectral density in (m^2)
    |     H: Wave height in (m)
    |     T: Wave period in (s)
    |     L: Wave length in (m)
    |     k: Wave number in (1/m) or (Radian/m)
    |     t: Time in (s)

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np

    windvel=10*np.random.rand(100)
    Fetch=10000*np.random.rand(100)
    dimValue=Fetch.copy()
    dimlessValue=sm.wavedim2dimless(windvel,dimValue,'fetch')

References
----------

Sverdrup, H. U., & Munk, W. H. (1947). 
Wind, sea, and swell: theory of relations for forecasting. 
U.S. Navy Department, Hydrographic Office, Publication No. 601, 44 pp. 

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
