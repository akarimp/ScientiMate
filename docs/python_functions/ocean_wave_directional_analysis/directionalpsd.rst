.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-05-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

scientimate.directionalpsd
==========================

.. code:: python

    Syy2d, f2d, theta = scientimate.directionalpsd(Syy, f, Wavedir=0, calcMethod='mitsuyasu', Windvel=0, dtheta=15, dispout='no')

Description
-----------

| Calculate wave directional spectrum using parametric directional spreading function
| such as Mitsuyasu et al. (1975), Hasselmann et al. (1980) and Donelan et al. (1985)

Inputs
------

Syy
    One dimensional power spectral density (m^2/Hz)
f
    Frequency (Hz)
Wavedir=0
    Mean wave direction between 0 and 360 (Degree)
CalcMethod='mitsuyasu'
    | Directional wave spectrum calculation method 
    | 'pierson': Pierson et al. (1955), 'cos2': D=1/pi*(cos((theta-theta_mean)/2))^2 
    | 'mitsuyasu': Mitsuyasu (1975), 'hasselmann': Hasselmann (1980), 'donelan': Donelan et al. (1985) 
    | 'flat': uniform distribution in all directions
Windvel=0
    | Wind velocity at 10 meter above surface level in (m/s)
    | Wind velocity is required for Mitsuyasu (1975) and Hasselmann (1980) methods
    | For other methods use Windvel=0
dtheta=15
    Direction interval at which directional spectrum calculated between 0 and 360 (Degree)
dispout='no'
    | Define to display outputs or not
    | '2d': 2 dimensional plot, 'surface': Surface plot, 'polar': Polar plot, 'no': not display 

Outputs
-------

Syy2d
    Directional wave power spectral density (m^2/Hz/Degree)
f2d
    Directional frequency (Hz)
theta
    Direction (Degree)

Examples
--------

.. code:: python

    import scientimate as sm
    import numpy as np

    N=2**11 #Total number of points
    fs=8 #Sampling frequency
    df=fs/N #Frequency difference 
    f=np.arange(0,fs/2+df,df) #Frequency vector 
    f[0]=f[1]/2 #Assign temporarily non-zero value to fisrt element of f to prevent division by zero
    Syy=0.016*9.81**2/((2*np.pi)**4*(f**5))*np.exp(-1.25*(0.33/f)**4) #Calculating Spectrum 
    f[0]=0
    Syy[0]=0
    Syy2d,f2d,theta=sm.directionalpsd(Syy,f,45,'mitsuyasu',10,15,'polar')

References
----------

Banner, M. L. (1990). 
Equilibrium spectra of wind waves. 
Journal of Physical Oceanography, 20(7), 966-984.

Donelan, M. A., Hamilton, J., & Hui, W. (1985). 
Directional spectra of wind-generated waves. 
Philosophical Transactions of the Royal Society of London A: Mathematical, Physical and Engineering Sciences, 315(1534), 509-562.

Ewans, K. C. (1998). 
Observations of the directional spectrum of fetch-limited waves. 
Journal of Physical Oceanography, 28(3), 495-512.

Goda, Y. (1999). 
A comparative review on the functional forms of directional wave spectrum. 
Coastal Engineering Journal, 41(01), 1-20.

Hasselmann, D. E., Dunckel, M., & Ewing, J. A. (1980). 
Directional wave spectra observed during JONSWAP 1973. 
Journal of physical oceanography, 10(8), 1264-1280.

Hwang, P. A., & Wang, D. W. (2001). 
Directional distributions and mean square slopes in the equilibrium and saturation ranges of the wave spectrum. 
Journal of physical oceanography, 31(5), 1346-1360.

Mitsuyasu, H., Tasai, F., Suhara, T., Mizuno, S., Ohkusu, M., Honda, T., & Rikiishi, K. (1975). 
Observations of the directional spectrum of ocean WavesUsing a cloverleaf buoy. 
Journal of Physical Oceanography, 5(4), 750-760.

Pierson, W. J., Neumann, G., & James, R. W. (1955). 
Practical methods for observing and forecasting ocean waves by means of wave spectra and statistics. 
Publication 603, U.S. Navy Hydrographic Office, 284 pp. 

Sorensen, R. M. (2005). 
Basic coastal engineering (Vol. 10). 
Springer Science & Business Media.

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
