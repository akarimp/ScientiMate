function [Kuv, f] = wavevel2wlconvfactor(f, h, heightfrombed, kCalcMethod, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-04-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

wavevel2wlconvfactor
====================

.. code:: MATLAB

    [Kuv, f] = wavevel2wlconvfactor(f, h, heightfrombed, kCalcMethod, dispout)

Description
-----------

Calculate a water particle horizontal orbital velocity to the water surface elevation conversion factor

Inputs
------

f
    Frequency in (Hz) 
h
    Water depth in (m)
heightfrombed
    Height from bed that Kp calculated at in (m)
kCalcMethod='beji';
    | Wave number calculation method 
    | 'hunt': Hunt (1979), 'beji': Beji (2013), 'vatankhah': Vatankhah and Aghashariatmadari (2013) 
    | 'goda': Goda (2010), 'exact': calculate exact value 
dispout='no';
    | Define to display outputs or not ('yes': display, 'no': not display)
    | Note: inputs can be as a single value or a 1-D vertical array

Outputs
-------

Kuv
    Horizontal orbital velocity to the water surface elevation conversion factor
f
    Frequency (Hz)

Examples
--------

.. code:: MATLAB

    [Kuv,f]=wavevel2wlconvfactor(0.2,1,0.2,'beji','no');

    [Kuv,f]=wavevel2wlconvfactor([0.2;0.25],[0.5;0.6],0.2,'exact','yes');

References
----------

Beji, S. (2013). 
Improved explicit approximation of linear dispersion relationship for gravity waves. 
Coastal Engineering, 73, 11-12.

Goda, Y. (2010). 
Random seas and design of maritime structures. 
World scientific.

Hunt, J. N. (1979). 
Direct solution of wave dispersion equation. 
Journal of the Waterway Port Coastal and Ocean Division, 105(4), 457-459.

Vatankhah, A. R., & Aghashariatmadari, Z. (2013). 
Improved explicit approximation of linear dispersion relationship for gravity waves: A discussion. 
Coastal engineering, 78, 21-22.

Wiberg, P. L., & Sherwood, C. R. (2008). 
Calculating wave-generated bottom orbital velocities from surface-wave parameters. 
Computers & Geosciences, 34(10), 1243-1262.

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
    case 3
        kCalcMethod='beji'; dispout='no';
    case 4
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(f)==1
    f=f';
end

if isrow(h)==1
    h=h';
end

%--------------------------------------------------------------------------
%Calculating wave number (k)

w=2.*pi.*f; %Wave angular frequency, w=2*pi/T=2*pi*f

%Deep water
k0=(w.^2)./9.81; %Deep water wave number
k0h=k0.*h;

%Estimation of wave number (k) from Hunt (1979)
if strcmp(kCalcMethod,'hunt')==1
    % c2gh=(k0h+(1.0+0.6522.*k0h+0.4622.*k0h.^2+0.0864.*k0h.^4+0.0675.*k0h.^5).^(-1)).^(-1); %Calculating wave number from Hunt (1979)
    % k=w./(sqrt(c2gh.*9.81*h));
    d1=0.6666666667; d2=0.3555555556; d3=0.1608465608; d4=0.0632098765; d5=0.0217540484; d6=0.0065407983;
    kh=sqrt(k0h.^2+(k0h)./(1+d1.*k0h+d2.*k0h.^2+d3.*k0h.^3+d4.*k0h.^4+d5.*k0h.^5+d6.*k0h.^6)); %Calculating wave number from Hunt (1979)
    k=kh./h; %Calculating wave number from Hunt (1979)

%Estimation of wave number (k) from Beji (2013)
elseif strcmp(kCalcMethod,'beji')==1
    kh=k0h.*(1+k0h.^1.09.*exp(-(1.55+1.3.*k0h+0.216.*k0h.^2)))./sqrt(tanh(k0h)); %Calculating wave number from Beji (2013)
    k=kh./h; %Calculating wave number from Beji (2013)
    k(w==0)=0;

%Estimation of wave number (k) from Vatankhah and Aghashariatmadari (2013)
elseif strcmp(kCalcMethod,'vatankhah')==1
    kh=(k0h+k0h.^2.*exp(-(3.2+k0h.^1.65)))./sqrt(tanh(k0h))+k0h.*(1-exp(-(k0h.^0.132))).^(5.0532+2.158.*k0h.^1.505); %Calculating wave number from Vatankhah and Aghashariatmadari (2013)
    k=kh./h; %Calculating wave number from Vatankhah and Aghashariatmadari (2013)
    k(w==0)=0;

%Estimation of wave number (k) from Goad (2010)
elseif strcmp(kCalcMethod,'goda')==1
    kh=zeros(length(k0h(:,1)),1);
    kh(k0h>=1)=k0h(k0h>=1);
    kh(k0h<1)=(k0h(k0h<1)).^0.5;
    for i=1:3
        kh=kh-((kh-k0h.*coth(kh))./(1+k0h.*((coth(kh)).^2-1))); %Calculating wave number from Goda (2010)
    end
    k=kh./h; %Calculating wave number from Goda (2010)
    k(w==0)=0;

%Calculating exact wave number (k)
elseif strcmp(kCalcMethod,'exact')==1

    %Estimation of wave number (k) from Beji (2013)
    kh=k0h.*(1+k0h.^1.09.*exp(-(1.55+1.3.*k0h+0.216.*k0h.^2)))./sqrt(tanh(k0h)); %Calculating wave number from Beji (2013)
    kini=kh./h; %Initial value for k (Wave number from Beji (2013))
    kini(w==0)=0;

    %Calculating exact wave number (k)
    k=zeros(length(f(:,1)),1); %Pre-assigning array to make program run faster
    for i=1:length(f(:,1))
        if length(h(:,1))==1
            fun = @(x)(w(i,1)^2-(9.81*x*tanh(x*h))); %w^2=g*k*tanh(kh)
        else
            fun = @(x)(w(i,1)^2-(9.81*x*tanh(x*h(i,1)))); %w^2=g*k*tanh(kh)
        end
        k(i,1)=fzero(fun,kini(i,1)); %Wave number
    end

end

%--------------------------------------------------------------------------
%Calculating a horizontal orbital velocity to the water surface elevation conversion factor

Kuv=(2*pi.*f).*cosh(k.*heightfrombed)./sinh(k.*h);

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1
    plot(f(f~=0),Kuv(f~=0))
    title('Orbital Velocity Conversion Factor vs Frequency')
    xlabel('Frequency (Hz)')
    ylabel('Orbital Velocity Conversion Factor, Kuv')
end

%--------------------------------------------------------------------------

