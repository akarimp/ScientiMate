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

wavepropfrompsd
===============

.. code:: MATLAB

    [Hm0, fp, Tp, Tm01, Tm02, m0, m1, m2] = wavepropfrompsd(Syy, f, fcL, fcH, dispout)

Description
-----------

Calculate wave properties from a power spectral density

Inputs
------

Syy
    Power spectral density (m^2/Hz)
f
    Frequency (Hz)
fcL=0;
    Low cut-off frequency, between 0*fs to 0.5*fs (Hz)
fcH=max(f);
    High cut-off frequency, between 0*fs to 0.5*fs (Hz)
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

Hm0
    Zero-Moment Wave Height (m)
fp
    Peak wave frequency (Hz)
Tp
    Peak wave period (second)
Tm01
    Wave Period from m01 (second), Mean Wave Period
Tm02
    Wave Period from m02 (second), Mean Zero Crossing Period
m0
    Zero-Moment of the power spectral density (m^2)
m1
    First-Moment of the power spectral density (m^2)
m2
    Second-Moment of the power spectral density (m^2)

Examples
--------

.. code:: MATLAB

    N=2^11; %Total number of points
    fs=8; %Sampling frequency
    df=fs/N; %Frequency difference 
    f(:,1)=[0:df:fs/2]; %Frequency vector 
    f(1,1)=f(2,1)/2; %Assign temporarily non-zero value to fisrt element of f to prevent division by zero
    Syy=0.016.*9.81.^2./((2.*pi).^4.*(f.^5)).*exp(-1.25.*(0.33./f).^4); %Calculating Spectrum 
    f(1,1)=0;
    Syy(1,1)=0;
    [Hm0,fp,Tp,Tm01,Tm02,m0,m1,m2]=wavepropfrompsd(Syy,f,0,8/2,'yes');

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
