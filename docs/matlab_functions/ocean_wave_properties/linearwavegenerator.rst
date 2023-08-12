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

linearwavegenerator
===================

.. code:: MATLAB

    [Eta, t, Etaij] = linearwavegenerator(amin, amax, Tmin, Tmax, Phimin, Phimax, fs, duration, NoOfWave, dispout)

Description
-----------

Generate linear waves

Inputs
------

amin
    Min wave amplitude in (m)
amax
    Max wave amplitude in (m)
Tmin
    Min wave mean period in (s)
Tmax
    Max wave mean period in (s)
Phimin=0;
    Min Phase (radian)
Phimax=2*pi;
    Max Phase (radian) 
fs=32;
    Sample generation frequency (Hz), number of data points in one second
duration=10;
    Duration time that data will be generated in (s)
NoOfWave=2;
    Number of waves to be combined with each other
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

Eta
    Water Surface Level Time Series in (m)
t
    Time in (s)
Etaij
    Separated Water Surface Level Time Series in (m)

Examples
--------

.. code:: MATLAB

    [Eta,t,Etaij]=linearwavegenerator(0.2,0.4,1,3,0,2*pi,32,10,2,'yes');

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
