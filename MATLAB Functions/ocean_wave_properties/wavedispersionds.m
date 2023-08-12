function [k, L, C, Cg] = wavedispersionds(h, T, Uc)
%{
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

wavedispersionds
================

.. code:: MATLAB

    [k, L, C, Cg] = wavedispersionds(h, T, Uc)

Description
-----------

| Solve water wave dispersion relation with presence of current (Doppler shift)
| Calculate wave number (k), wave length (L), wave celereity (C), and wave group velocity (Cg) using linear wave theory

Inputs
------

h
    Water depth in (m)
T
    | Wave period in (s) 
    | If peak wave frequency (Tp) is used, calculated values represent peak wave 
Uc=0;
    | Current velocity in (m/s), for Doppler shift
    | Uc should have a same size as h
    | Note: inputs can be as a single value or a 1-D vertical array

Outputs
-------

k
    Wave number in (radian/m)
L
    Wave length in (m)
C
    Wave celerity in (m/s)
Cg
    Wave group celerity in (m/s)

Examples
--------

.. code:: MATLAB

    [k,L,C,Cg]=wavedispersionds(1,3,1);

    [k,L,C,Cg]=wavedispersionds([1;1.1],[3;3.1],[1;1]);

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
%}

%--------------------------------------------------------------------------
%CODE
%--------------------------------------------------------------------------
%Assign default values

switch nargin
    case 2
        Uc=0;
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(h)==1
    h=h';
end

if isrow(T)==1
    T=T';
end

if isrow(Uc)==1
    Uc=Uc';
end

%--------------------------------------------------------------------------
%Calculating wave number (k)

%Make Uc the same size as h
if length(Uc(:,1))==1
    Uc1=Uc;
    Uc=zeros(length(h(:,1)),1);
    Uc(:,1)=Uc1;
end

f=1./T; %Wave frequency
w=2.*pi.*f; %Wave angular frequency, w=2*pi/T=2*pi*f

%Deep water
k0=(w.^2)./9.81; %Deep water wave number
k0h=k0.*h;

%Calculating exact wave number (k)

%Estimation of wave number (k) from Beji (2013)
kh=k0h.*(1+k0h.^1.09.*exp(-(1.55+1.3.*k0h+0.216.*k0h.^2)))./sqrt(tanh(k0h)); %Calculating wave number from Beji (2013)
kini=kh./h; %Initial value for k (Wave number from Beji (2013))
kini(w==0)=0;

%Calculating exact wave number (k)
k=zeros(length(f(:,1)),1); %Pre-assigning array to make program run faster
for i=1:length(f(:,1))
    if length(h(:,1))==1
        fun = @(x)((w(i,1)-x*Uc)^2-(9.81*x*tanh(x*h))); %(w-kU)^2=g.k.tanh(kd)
    else
        fun = @(x)((w(i,1)-x*Uc(i,1))^2-(9.81*x*tanh(x*h(i,1)))); %(w-kU)^2=g.k.tanh(kd)
    end
    k(i,1)=fzero(fun,kini(i,1)); %Wave number
end

%--------------------------------------------------------------------------
%Calculating wave length (L), wave celereity (C), and wave group velocity (Cg)

L=2.*pi./k; %Wave length
C=9.81.*T./(2.*pi).*tanh(k.*h); %Wave celerity
Cg=1./2.*C.*(1+(2.*k.*h)./sinh(2.*k.*h)); %Wave group celirity

%--------------------------------------------------------------------------
