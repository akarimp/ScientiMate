function [Eta, t, Etaij] = stokeswavegenerator(h, amin, amax, Tmin, Tmax, Phimin, Phimax, fs, duration, NoOfWave, kCalcMethod, dispout)
%{
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

stokeswavegenerator
===================

.. code:: MATLAB

    [Eta, t, Etaij] = stokeswavegenerator(h, amin, amax, Tmin, Tmax, Phimin, Phimax, fs, duration, NoOfWave, kCalcMethod, dispout)

Description
-----------

Generate second order stokes' waves

Inputs
------

h
    Mean water depth in (m)
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
kCalcMethod='beji';
    | Wave number calculation method 
    | 'hunt': Hunt (1979), 'beji': Beji (2013), 'vatankhah': Vatankhah and Aghashariatmadari (2013) 
    | 'goda': Goda (2010), 'exact': calculate exact value 
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

    [Eta,t,Etaij]=stokeswavegenerator(5,0.2,0.4,1,3,0,2*pi,32,10,2,'beji','yes');

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
    case 5
        Phimin=0; Phimax=2*pi; fs=32; duration=10; NoOfWave=2; kCalcMethod='beji'; dispout='no';
    case 6
        Phimax=2*pi; fs=32; duration=10; NoOfWave=2; kCalcMethod='beji'; dispout='no';
    case 7
        fs=32; duration=10; NoOfWave=2; kCalcMethod='beji'; dispout='no';
    case 8
        duration=10; NoOfWave=2; kCalcMethod='beji'; dispout='no';
    case 9
        NoOfWave=2; kCalcMethod='beji'; dispout='no';
    case 10
        kCalcMethod='beji'; dispout='no';
    case 10
        dispout='no';
end

%--------------------------------------------------------------------------

sample=fs*duration; %number of sample in input file
dt=1/fs;
t(:,1)=linspace(dt,duration,sample); %time
a(:,1)=linspace(amin,amax,NoOfWave);
T(:,1)=linspace(Tmin,Tmax,NoOfWave);
Phi(:,1)=linspace(Phimin,Phimax,NoOfWave);
H=2.*a; % Wave height

f=1./T; %Wave frequency
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

L=2*pi./k;

% Generating Waves
FirstOrder=zeros(length(t(:,1)),NoOfWave); %Pre-assigning array to make program run faster
SecondOrder=zeros(length(t(:,1)),NoOfWave); %Pre-assigning array to make program run faster
for i=1:NoOfWave
    FirstOrder(:,i)=(H(i,1)./2).*cos(-w(i,1).*t+Phi(i,1));
    SecondOrder(:,i)=(pi.*H(i,1)./8).*(H(i,1)./L(i,1)).*(cosh(k(i,1).*h)).*(2+cosh(2.*k(i,1).*h))./((sinh(k(i,1).*h)).^3).*(cos(-2.*w(i,1).*t));
    % SecondOrder(:,i)=((H(i,1)).^2.*k(i,1)./16).*(cosh(k(i,1).*h)).*(2+cosh(2.*k(i,1).*h))./((sinh(k(i,1).*h)).^3).*(cos(-2.*w(i,1).*t));
end

Etaij=FirstOrder+SecondOrder;

Eta=sum(Etaij,2);

%for i=1:length(t(:,1))
%    Eta(i,1)=sum(Etaij(i,1:end));
%end

% -------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1
    
    subplot(2,1,1)
    for i=1:NoOfWave
        plot(t,Etaij(:,i))
        hold on
    end
    title('Generated Waves')
    xlabel('Time (s)')
    ylabel('Water Level (m)')
    xlim([t(1,1) t(end,1)])
    
    subplot(2,1,2)
    plot(t,Eta)
    hold on
    title('Generated Combined Wave')
    xlabel('Time (s)')
    ylabel('Water Level (m)')
    % legend('Water Level (m)')
    xlim([t(1,1) t(end,1)])
    
end

% -------------------------------------------------------------------------
