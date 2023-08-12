function [dimlessValue] = wavedim2dimless(windvel, dimValue, ValueType)
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

wavedim2dimless
===============

.. code:: MATLAB

    [dimlessValue] = wavedim2dimless(windvel, dimValue, ValueType)

Description
-----------

Calculate dimensionless numbers from dimensional numbers

Inputs
------

windvel
    Wind velocity in (m/s)
dimValue
    Dimensional value to be converted to dimensionless value
ValueType='waveheight';
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

.. code:: MATLAB

    windvel(:,1)=10.*rand(100,1);
    Fetch(:,1)=10000.*rand(100,1);
    dimValue=Fetch;
    [dimlessValue]=wavedim2dimless(windvel,dimValue,'fetch');

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

if isrow(dimValue)==1
    dimValue=dimValue';
end

%--------------------------------------------------------------------------
%Calculate dimensionless numbers from dimensional number

%Calculating dimensionless Fetch
if strcmp(ValueType,'fetch')==1
    Fetch=dimValue;
    Fetchhat=9.81.*Fetch./(windvel.^2); %Dimensionless fetch: Fetchhat=g*Fetch/U10^2
    dimlessValue=Fetchhat;

%Calculating dimensionless water depth
elseif strcmp(ValueType,'depth')==1
    h=dimValue;
    hhat=9.81.*h./(windvel.^2); %Dimensionless water depth: hhat=g*h/U10^2
    dimlessValue=hhat;

%Calculating dimensionless wave frequency
elseif strcmp(ValueType,'frequency')==1
    f=dimValue;
    fhat=f.*windvel./9.81; %Dimensionless wave frequency: fhat=f*U10/g
    dimlessValue=fhat;

%Calculating dimensionless wave energy
elseif strcmp(ValueType,'energy')==1
    m0=dimValue;
    Ehat=9.81^2.*m0./(windvel.^4); %Dimensionless wave energy: Ehat=g^2*m0/U10^4
    dimlessValue=Ehat;

%Calculating dimensionless wave height
elseif strcmp(ValueType,'waveheight')==1
    H=dimValue;
    Hhat=9.81.*H./(windvel.^2); %Dimensionless wave height: Hhat=g*H/U10^2
    dimlessValue=Hhat;

%Calculating dimensionless wave period
elseif strcmp(ValueType,'period')==1
    T=dimValue;
    That=9.81.*T./windvel; %Dimensionless wave period: That=g*T/U10
    dimlessValue=That;

%Calculating dimensionless wave length
elseif strcmp(ValueType,'wavelength')==1
    L=dimValue;
    Lhat=9.81.*L./(windvel.^2); %Dimensionless wave length: Lhat=g*L/U10^2
    dimlessValue=Lhat;

%Calculating dimensionless wave number
elseif strcmp(ValueType,'wavenumber')==1
    k=dimValue;
    khat=k.*(windvel.^2)./9.81; %Dimensionless wave number: khat=k*U10^2/g
    dimlessValue=khat;

%Calculating dimensionless time
elseif strcmp(ValueType,'time')==1
    t=dimValue;
    that=9.81.*t./windvel; %Dimensionless time: that=g*t/U10
    dimlessValue=that;

%Calculating dimensionless length
elseif strcmp(ValueType,'length')==1
    x=dimValue;
    xhat=9.81.*x./(windvel.^2); %Dimensionless length: xhat=g*x/U10^2
    dimlessValue=xhat;

end

%--------------------------------------------------------------------------
