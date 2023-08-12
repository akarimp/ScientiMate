.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-02-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

replacespike3dps
================

.. code:: MATLAB

    [xDespiked, Indx] = replacespike3dps(x, nrpeat, DespikeScaleFactor, interpMethod, dispout)

Description
-----------

Remove spikes in the time series based on 3D phase space method by Goring and Nikora (2002)

Inputs
------

x
    Input data
nrpeat=1;
    Number of time despiking procedure is repeating
DespikeScaleFactor=1;
    Scaling a despiking threshold
interpMethod='linear';
    | Interpolation method for replacing spike points:
    | Matlab/Octave: 'linear', 'nearest', 'next', 'previous', 'pchip', 'cubic', 'spline'
    | Python/Scipy : 'linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic'
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

xDespiked
    Dispiked data
Indx
    Index of despiked points

Examples
--------

.. code:: MATLAB

    fs=128;
    t(:,1)=linspace(0,9.5,10*fs);
    x(:,1)=sin(2*pi*0.3*t)+0.1*sin(2*pi*4*t);
    spikeloc(:,1)=[10:100:length(t(:,1))];
    x(spikeloc+round(2*randn(length(spikeloc(:,1)),1)),1)=sign(randn(length(spikeloc(:,1)),1));
    x(220:225,1)=1.5;
    x=x+5;
    [xDespiked,Indx]=replacespike3dps(x,2,1,'linear','yes');

    fs=2;
    t(:,1)=linspace(0,1023.5,1024*fs);
    x(:,1)=detrend(0.5.*cos(2*pi*0.2*t)+(-0.1+(0.1-(-0.1))).*rand(1024*fs,1));
    spikeloc(:,1)=[10:100:length(t(:,1))];
    x(spikeloc+round(2*randn(length(spikeloc(:,1)),1)),1)=sign(randn(length(spikeloc(:,1)),1));
    [xDespiked,Indx]=replacespike3dps(x,1,1,'linear','yes');

References
----------

Goring, D. G., & Nikora, V. I. (2002). 
Despiking acoustic Doppler velocimeter data. 
Journal of Hydraulic Engineering, 128(1), 117-126.

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
