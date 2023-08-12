function [Uz2_PowerLaw, Uz2_LogLaw, Uz2_LogLaw_Charnock] = windvelz1toz2(Uz1, z1, z2, alpha, z0, alpha_Charnock)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-07-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

windvelz1toz2
=============

.. code:: MATLAB

    [Uz2_PowerLaw, Uz2_LogLaw, Uz2_LogLaw_Charnock] = windvelz1toz2(Uz1, z1, z2, alpha, z0, alpha_Charnock)

Description
-----------

Convert wind velocity from first height, z1 (m), to second height, z2 (m), (e.g. 10 (m)) above surface

Inputs
------

Uz1
    Wind velocity at z1 Height in (m/s)
z1
    First Height (e.g. sensor height) in (m)
z2=10.0;
    Second Height that wind velocity is converted to in (m)
alpha=1/7;
    | Exponent of power law velocity profile, Uz1/Uz2=(z1/z2)**(alpha)
    | water: alpha=1/11.5, open area: alpha=1/9.5, suburban area: alpha=1/7, urban area: alpha=1/5 
z0=0.005;
    | Surface roughness in logarithmic velocity profile in (m), U/Ustar=(1/0.41)*ln(z/z0)
    | water: z0=0.005 (m), open area: z0=0.02 (m), suburban area: z0=0.3 (m), urban area: z0=2 (m) 
alpha_Charnock=0.0145;
    | Charnock (1955) parameters to calculate ocean surface roughness, z0=alpha_Charnock*(Ustar^2/g) 
    | Charnock (1955): alpha_Charnock=0.012, Garratt (1977): alpha_Charnock=0.0145, Wu (1980): alpha_Charnock=0.0185

Outputs
-------

Uz2_PowerLaw
    | Wind velocity at height of z2 above surface in (m/s)
    | calculated from power law velocity profile
Uz2_LogLaw
    | Wind velocity at height of z2 above surface in (m/s)
    | calculated from logarithmic velocity profile
Uz2_LogLaw_Charnock
    | Wind velocity at height of z2 above surface in (m/s)
    | calculated from logarithmic velocity profile 
    | using Charnock (1955) relationship for wind over ocean

Examples
--------

.. code:: MATLAB

    [Uz2_PowerLaw,Uz2_LogLaw,Uz2_LogLaw_Charnock]=windvelz1toz2(8,15,10);

    Uz1(:,1)=10.*rand(100,1);
    [Uz2_PowerLaw,Uz2_LogLaw,Uz2_LogLaw_Charnock]=windvelz1toz2(Uz1,15,10,1/7,0.005,0.0145);

References
----------

Charnock, H. (1955). 
Wind stress on a water surface. 
Quarterly Journal of the Royal Meteorological Society, 81(350), 639-640.

Garratt, J. R. (1977). 
Review of drag coefficients over oceans and continents. 
Monthly weather review, 105(7), 915-929.

Wu, J. (1980). 
Wind-stress coefficients over sea surface near neutral conditions—A revisit. 
Journal of Physical Oceanography, 10(5), 727-740.

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
        z2=10; alpha=1/7; z0=0.005; alpha_Charnock=0.0145;
    case 3
        alpha=1/7; z0=0.005; alpha_Charnock=0.0145;
    case 4
        z0=0.005; alpha_Charnock=0.0145;
    case 5
        alpha_Charnock=0.0145;
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(Uz1)==1
    Uz1=Uz1';
end

%--------------------------------------------------------------------------
%Converting wind velocity from z1 (m) height to z2 (m) height using power law velocity profile

Uz2_PowerLaw=Uz1.*(z2/z1).^(alpha);

%--------------------------------------------------------------------------
%Converting wind velocity from z1 (m) height to z2 (m) height using logarithmic velocity profile

%Calculating shear velocity, Ustar
Ustar=zeros(length(Uz1(:,1)),1); %Pre-assigning array
for i=1:length(Uz1(:,1))
    fun = @(x)(Uz1(i,1)-(x/0.41)*log(z1/z0)); %Logarithmic velocity profile, U=(Ustar/0.41)*ln(z/z0)
    Ustar(i,1)=fzero(fun,Uz1(i,1));
end

Uz2_LogLaw=(Ustar./0.41).*log(z2./z0); %Calculating wind velocity at z2 m height 

%--------------------------------------------------------------------------
%Converting wind velocity from z1 (m) height to z2 (m) height using Logarithmic velocity profile by using Charnock (1955) relationship

%Calculating shear velocity, Ustar, based on Charnock (1955) relationship: z0=alpha*(Ustar^2/g)
UstarCh=zeros(length(Uz1(:,1)),1); %Pre-assigning array
for i=1:length(Uz1(:,1))
    fun = @(x)(Uz1(i,1)-(x/0.41)*log(z1/(alpha_Charnock*x^2/9.81))); %Logarithmic velocity profile, U=(Ustar/0.41)*ln(z/z0)
    UstarCh(i,1)=fzero(fun,Uz1(i,1));
end

z0Ch=alpha_Charnock.*UstarCh.^2/9.81; %Charnock (1955) relationship: z0=alpha_Charnock*(Ustar^2/g)

Uz2_LogLaw_Charnock=UstarCh./0.41.*log(z2./z0Ch); %Calculating wind velocity at z2 m height 

%--------------------------------------------------------------------------
