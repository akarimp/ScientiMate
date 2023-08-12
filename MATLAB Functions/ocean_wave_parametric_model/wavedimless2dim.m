function [dimValue] = wavedimless2dim(windvel, dimlessValue, ValueType)
%{
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

wavedimless2dim
===============

.. code:: MATLAB

    [dimValue] = wavedimless2dim(windvel, dimlessValue, ValueType)

Description
-----------

Calculate dimensional numbers from dimensionless numbers

Inputs
------

windvel
    Wind velocity in (m/s)
dimlessValue
    Dimensionless value to be converted to dimensional value
ValueType='waveheight';
    | Type of the dimensionless value 
    | 'fecth': Dimensionless wind fetch (Fetchhat): Fetchhat=g*Fetch/U10^2
    | 'depth':Dimensionless depth (hhat): hhat=g*h/U10^2
    | 'frequency':Dimensionless wave frequency (fhat): fhat=f*U10/9.81
    | 'energy':Dimensionless wave energy (Ehat): Ehat=g^2*m0/U10^4
    | 'waveheight':Dimensionless wave height (Hhat): Hhat=g*H/U10^2
    | 'period':Dimensionless wave period (That): That=g*T/U10
    | 'wavelength':Dimensionless wave length (Lhat): Lhat=g*L/U10^2
    | 'wavenumber':Dimensionless wave number (khat): khat=k*U10^2/g
    | 'time':Dimensionless time (that): that=g*t/U10
    | 'length':Dimensionless length (xhat): xhat=g*x/U10^2

Outputs
-------

dimValue
    | Dimensional value calculated from dimensionless value
    | For ValueType='fetch': dimValue is wind fetch in (m)
    | For ValueType='depth': dimValue is water depth in (m)
    | For ValueType='frequency': dimValue is wave frequency in (Hz)
    | For ValueType='energy': dimValue is zero-moment of surface elevation power spectral density in (m^2)
    | For ValueType='waveheight': dimValue is wave height in (m)
    | For ValueType='period': dimValue is wave period in (s)
    | For ValueType='wavelength': dimValue is wave length in (m)
    | For ValueType='wavenumber': dimValue is wave number in (1/m) or (Radian/m)
    | For ValueType='time': dimValue is time in (s)
    | For ValueType='length': dimValue is length in (m)
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
    |     x: Length in (m)

Examples
--------

.. code:: MATLAB

    windvel(:,1)=10.*rand(100,1);
    Fetchhat(:,1)=10000.*rand(100,1);
    dimlessValue=Fetchhat;
    [dimValue]=wavedimless2dim(windvel,dimlessValue,'fetch');

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
%}

%--------------------------------------------------------------------------
%CODE
%--------------------------------------------------------------------------
%Assign default values

switch nargin
    case 2
        ValueType='waveheight';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(windvel)==1
    windvel=windvel';
end

if isrow(dimlessValue)==1
    dimlessValue=dimlessValue';
end

%--------------------------------------------------------------------------
%Calculate dimensional numbers from dimensionless number

%Calculating dimensionless Fetch
if strcmp(ValueType,'fetch')==1
    Fetchhat=dimlessValue;
    Fetch=Fetchhat.*(windvel.^2)./9.81; %Wind fetch: Fetchhat=g*Fetch/U10^2
    dimValue=Fetch;

%Calculating dimensionless water depth
elseif strcmp(ValueType,'depth')==1
    hhat=dimlessValue;
    h=hhat.*(windvel.^2)./9.81; %Depth: hhat=g*h/U10^2
    dimValue=h;

%Calculating dimensionless wave frequency
elseif strcmp(ValueType,'frequency')==1
    fhat=dimlessValue;
    f=fhat.*9.81./(windvel); %Wave frequency, fhat=f*U10/g
    dimValue=f;

%Calculating dimensionless wave energy
elseif strcmp(ValueType,'energy')==1
    Ehat=dimlessValue;
    m0=Ehat.*(windvel.^4)./(9.81^2); %Zero-Moment, Ehat=g^2*m0/U10^4
    dimValue=m0;

%Calculating dimensionless wave height
elseif strcmp(ValueType,'waveheight')==1
    Hhat=dimlessValue;
    H=Hhat.*(windvel.^2)./9.81; %Wave height: Hhat=g*H/U10^2
    dimValue=H;

%Calculating dimensionless wave period
elseif strcmp(ValueType,'period')==1
    That=dimlessValue;
    T=That.*windvel./9.81; %Wave period: That=g*T/U10
    dimValue=T;

%Calculating dimensionless wave length
elseif strcmp(ValueType,'wavelength')==1
    Lhat=dimlessValue;
    L=Lhat.*(windvel.^2)./9.81; %Wave length: Lhat=g*L/U10^2
    dimValue=L;

%Calculating dimensionless wave number
elseif strcmp(ValueType,'wavenumber')==1
    khat=dimlessValue;
    k=khat.*9.81./(windvel.^2); %Wave number: khat=k*U10^2/g
    dimValue=k;

%Calculating dimensionless time
elseif strcmp(ValueType,'time')==1
    that=dimlessValue;
    t=that.*windvel./9.81; %Time: that=g*t/U10
    dimValue=t;

%Calculating dimensionless length
elseif strcmp(ValueType,'length')==1
    xhat=dimlessValue;
    x=xhat.*(windvel.^2)./9.81; %Length: xhat=g*x/U10^2
    dimValue=x;

end

%--------------------------------------------------------------------------
