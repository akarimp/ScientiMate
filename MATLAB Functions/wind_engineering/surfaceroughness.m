function [ustar, z0, d] = surfaceroughness(z, u, delta, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2018-03-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

surfaceroughness
================

.. code:: MATLAB

    [ustar, z0, d] = surfaceroughness(z, u, delta, dispout)

Description
-----------

Calculate shear velocity and surface roughness from a given velocity profile using Karimpour et al. (2012) method

Inputs
------

z
    Distance from a surface (elevation, height) in (m)
u
    Velocity at z in (m/s)
delta=max(z);
    Boundary layer height in (m)
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

z0
    Surface roughness in (m)
ustar
    Shear velocity (u*) in (m/s)
d
    | Zero plane displacement distance in (m)
    | Note: Above values are for a logarithmic velocity profile as:
    |     u=(u*/K)*ln((z-d)/z0)

Examples
--------

.. code:: MATLAB

    z(:,1)=[0.1:0.05:1];
    u(:,1)=2/0.4*log((z-0.003)/0.002);
    [ustar,z0,d]=surfaceroughness(z,u,max(z),'yes');

References
----------

Karimpour, A., Kaye, N. B., & Baratian-Ghorghi, Z. (2012). 
Modeling the neutrally stable atmospheric boundary layer for laboratory scale studies of the built environment. 
Building and Environment, 49, 203-211.

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
        delta=max(z); dispout='no';
    case 3
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(u)==1
    u=u';
end

if isrow(z)==1
    z=z';
end

%--------------------------------------------------------------------------
%Calculating shear velocity (friction velocity)

umaxminusu=max(u)-u; %umax-u
z_delta=z/delta;

z_delta_min=0.1; %Lower limit that data used in curve fitting
z_delta_max=0.5; %Upper limit that data used in curve fitting

umaxminusu=umaxminusu(z_delta>=z_delta_min & z_delta<=z_delta_max);
z_delta=z_delta(z_delta>=z_delta_min & z_delta<=z_delta_max);

x=z_delta;
y=umaxminusu;
p=polyfit(log(x),y,1);

ustar=-0.4*p(1); %Shear velocity (friction velocity), u*

%--------------------------------------------------------------------------
%Calculating surface roughness and zero plane displacement distance

exp_u_k_ustar=exp(u.*(0.4/ustar));

x=z;
y=exp_u_k_ustar;
p=polyfit(x,y,1);

z0=1/p(1); %Surface roughness, z0
d=abs(p(2)*z0); %Zero plane displacement distance

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1

    ufit=ustar/0.4*log((z-d)/z0);

    scatter(u,z)
    hold on
    plot(ufit,z)
    xlabel('Velocity')
    ylabel('Elevation')
    legend('Data','Fitted Curve')

end

%--------------------------------------------------------------------------
