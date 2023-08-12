function [Eta, x, maxsurgeheight, L] = stormsurge1d(h0, U10, m, dispout)
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

stormsurge1d
============

.. code:: MATLAB

    [Eta, x, maxsurgeheight, L] = stormsurge1d(h0, U10, m, dispout)

Description
-----------

Calculate one dimensional storm surge using Dean Dalrymple (1991) method

Inputs
------

h0
    Deep-water depth in (m), 
U10
    10-min averaged wind velocity at 10 meter above a surface in (m/s)
m=0;
    | Bed slope
    | Note: m=h0/L, where L is a length of the continental shelf in (m)
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

Eta
    Surge height along a x axis in (m)
x
    | Points on x axis in (m) 
    | x=0 located at a far-end of the water boundary
    | x=L located at a coastline
maxsurgeheight
    Maximum surge height at a coastline (x=L) in (m)
L
    | Length of the continental shelf in (m)
    | L=h0/m if m is not zero
    | L=50000 if m is zero (L=50000 meter for the hurricane Katrina)

Examples
--------

.. code:: MATLAB

    h0=42; %Example: h0=42 m for the hurricane Katrina
    U10=40; %Example: U10=40 m/s for the hurricane Katrina 
    m=0.00084;
    [Eta,x,maxsurgeheight,L]=stormsurge1d(h0,U10,m,'yes');

References
----------

Dean, R. G., & Dalrymple, R. A. (1991). 
Water wave mechanics for engineers and scientists (Vol. 2). 
World Scientific Publishing Company.

Wu, J. (1982). 
Wind‐stress coefficients over sea surface from breeze to hurricane. 
Journal of Geophysical Research: Oceans, 87(C12), 9704-9706.

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
        m=0; dispout='no';
    case 3
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

%if isrow(U10)==1
%    U10=U10';
%end

%--------------------------------------------------------------------------
%Initial values

n=1.3; %Shear stress coefficient
Row=1028; %Seawater density in (kg/m3)
Rhoa=1.204; %Air density in (kg/m3)

%Calculating a length of the continental shelf
if m==0
    L=50000; %L=50000 m for the hurricane Katrina
else
    L=h0/m;
end

%--------------------------------------------------------------------------
%Calculating wind drag coefficient based on Wu (1982), method used by ADCIRC

CD=(0.8+0.065*U10).*1e-3; %Wind drag coefficient
CD(U10<7.5)=1.2875.*1e-3; %Wind drag coefficient
%CD(CD>0.0034)=0.0034; %Wind drag coefficient
Tauwind=Rhoa.*CD.*U10.^2; %Wind shear stress in (N/m^2)
A=n*Tauwind*L/(Row*9.81*h0^2);

%--------------------------------------------------------------------------
%Calculating the surge

x(:,1)=linspace(0,L,1000); %Points on a x axis

Eta=zeros(length(x(:,1)),1); %Pre-assign array
%Bed without slope
if m==0
    
    %Claculating Eta for a bed with a zero slope
    Eta=h0.*(sqrt(1+2*A.*x/L)-1);
    
%Bed with slpoe
else
    
    h=-m.*x+h0; %h=h0.*(1-x./L);
    
    %Claculating Eta for a bed with a slope
    Eta(1,1)=0;
    
    for i=2:length(x(:,1))
        
        Etaini=Eta(i-1); %Initial value for Eta from previous step
        fun = @(Eta1)(x(i,1)/L-((1-(h(i,1)+Eta1)/h0)-A*log(((h(i,1)+Eta1)/h0-A)/(1-A)))); %x/L=(1-(h+Eta)/h0)-A*log(((h+Eta)/h0-A)/(1-A))
        Eta(i,1)=fzero(fun,Etaini); %Surge height
        
    end
    
end

maxsurgeheight=Eta(end,1);

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1

    %Colors: https://htmlcolorcodes.com/
    
    if m==0
        
        %Plotting a seabed
        x1=[0;L;L;L+L/10];
        y1=[-h0;-h0;0;0];
        hold on
        plot(x1,y1,'Color',[160,64,0]./255)
        
        %Plotting a water surface
        x2=[0;L];
        y2=[0;0];
        plot(x2,y2,'Color',[41,128,185]./255)
        
        %Plotting a storm surge
        plot(x,Eta,'Color',[231,76,60]./255)

        xlim([0 L+L/10])
        xlabel('x')
        ylabel('Elevation')
        leg1=legend('Seabed','Water Surface','Storm Surge');
        set(leg1,'Location','SouthEast')    
        
    else
        
        %Plotting a seabed
        x1=[0;L;L+L/10];
        y1=[-h0;m*L-h0;m*L-h0];
        hold on
        plot(x1,y1,'Color',[160,64,0]./255)
        
        %Plotting a water surface
        x2=[0;L];
        y2=[0;0];
        plot(x2,y2,'Color',[41,128,185]./255)

        %Plotting a storm surge
        plot(x,Eta,'Color',[231,76,60]./255)

        xlim([0 L+L/10])
        xlabel('x')
        ylabel('Elevation')
        leg1=legend('Seabed','Water Surface','Storm Surge');
        set(leg1,'Location','SouthEast')    

    end

end

%--------------------------------------------------------------------------
