function [CD, Tauwind, Ustar] = winddrag(U10, CalcMethod, Rhoa, dispout)
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

winddrag
========

.. code:: MATLAB

    [CD, Tauwind, Ustar] = winddrag(U10, CalcMethod, Rhoa, dispout)

Description
-----------

Calculate wind drag coefficient, wind shear stress, and wind shear velocity

Inputs
------

U10
    | Wind velocity in (m/s)
    | Wind velocity should be measured (or represents velocity) at 10 m above surface
kCalcMethod='beji';
    | Drag coefficient calculation method 
    | 'vandoren': Van Doren (1953)
    | 'garratt': Garratt (1977)
    | 'smith': Smith (1977)
    | 'large': Large & Pond (1981)
    | 'wu': Wu (1982)
    | 'hwang': Hwang (2011)
    | 'zijlema': Zijlema et al. (2012)
    | 'cem': Use method by Coastal Engineering Manual (CEM),
    |     U.S. Army Corps of Engineers (2015)
Rhoa=1.204;
    Air density in (kg/m3)
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

CD
    Wind drag coefficient
Tauwind
    Wind shear stress in (N/m^2), Tauwind=Rhoa*CD*U10^2;
Ustar
    Wind shear velocity in (m/s), Ustar=(CD*(U10^2))^0.5;

Examples
--------

.. code:: MATLAB

    %U10(:,1)=[4:0.1:40];
    %[CD,Tauwind,Ustar]=winddrag(U10,'wu',1.204,'yes');

References
----------

Bryant, K. M., & Akbar, M. (2016). 
An Exploration of Wind Stress Calculation Techniques in Hurricane Storm Surge Modeling. 
Journal of Marine Science and Engineering, 4(3), 58.

Chen, Y., & Yu, X. (2017). 
Sensitivity of storm wave modeling to wind stress evaluation methods. 
Journal of Advances in Modeling Earth Systems.

Dean, R. G., & Dalrymple, R. A. (1991). 
Water wave mechanics for engineers and scientists (Vol. 2). 
world scientific publishing Co Inc.

Garratt, J. R. (1977). 
Review of drag coefficients over oceans and continents. 
Monthly weather review, 105(7), 915-929.

Hwang, P. A. (2011). 
A note on the ocean surface roughness spectrum. 
Journal of Atmospheric and Oceanic Technology, 28(3), 436-443.

Large, W. G., & Pond, S. (1981). 
Open ocean momentum flux measurements in moderate to strong winds. 
Journal of physical oceanography, 11(3), 324-336.

Rogers, W. E., Babanin, A. V., & Wang, D. W. (2012). 
Observation-consistent input and whitecapping dissipation in a model for wind-generated surface waves: 
Description and simple calculations. 
Journal of Atmospheric and Oceanic Technology, 29(9), 1329-1346.

Smith, S. D. (1980). 
Wind stress and heat flux over the ocean in gale force winds. 
Journal of Physical Oceanography, 10(5), 709-726.

U.S. Army Corps of Engineers (2015). 
Coastal Engineering Manual. 
Engineer Manual 1110-2-1100, Washington, D.C.: U.S. Army Corps of Engineers.

VanDorn, W. G. (1953). 
Wind stress on an artificial pond. 
Journal of Marine Research, 12(3), 249-276.

Wu, J. (1982). 
Wind‐stress coefficients over sea surface from breeze to hurricane. 
Journal of Geophysical Research: Oceans, 87(C12), 9704-9706.

Zijlema, M., Van Vledder, G. P., & Holthuijsen, L. H. (2012). 
Bottom friction and wind drag for wave models. 
Coastal Engineering, 65, 19-26.

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
    case 1
        CalcMethod='wu'; Rhoa=1.204; dispout='no';
    case 2
        Rhoa=1.204; dispout='no';
    case 2
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(U10)==1
    U10=U10';
end

if isrow(Rhoa)==1
    Rhoa=Rhoa';
end

%--------------------------------------------------------------------------
%Calculating wind drag coefficient

%Calculating wind drag coefficient based on Van Doren (1953) from Dean & Dalrymple (1991) page 157
if strcmp(CalcMethod,'vandoren')==1

    CD=(1.2e-6+2.25e-6.*(1-5.6./U10).^2).*1e3; %Wind drag coefficient
    CD(U10<=5.6)=(1.2e-6)*1e3; %Wind drag coefficient
    Tauwind=Rhoa.*CD.*U10.^2; %Wind shear stress in (N/m^2)
    Ustar=(CD.*(U10.^2)).^0.5; %Wind shear velocity in (m/s)

%Calculating wind drag coefficient based on Garratt (1977), method used by SWAN
elseif strcmp(CalcMethod,'garratt')==1

    CD=(0.75+0.067.*U10).*1e-3; %Wind drag coefficient
    CD(U10>41)=3.5.*1e-3; %Wind drag coefficient
    %CD(CD>0.0034)=0.0034; %Wind drag coefficient
    Tauwind=Rhoa.*CD.*U10.^2; %Wind shear stress in (N/m^2)
    Ustar=(CD.*(U10.^2)).^0.5; %Wind shear velocity in (m/s)

%Calculating wind drag coefficient based on Smith (1980)
elseif strcmp(CalcMethod,'smith')==1

    CD=(0.61+0.063.*U10).*1e-3; %Wind drag coefficient
    Tauwind=Rhoa.*CD.*U10.^2; %Wind shear stress in (N/m^2)
    Ustar=(CD.*(U10.^2)).^0.5; %Wind shear velocity in (m/s)

%Calculating wind drag coefficient based on Large & Pond (1981)
elseif strcmp(CalcMethod,'large')==1

    CD=(0.49+0.065.*U10).*1e-3; %Wind drag coefficient
    CD(U10<=10)=1.14*1e-3; %Wind drag coefficient
    Tauwind=Rhoa.*CD.*U10.^2; %Wind shear stress in (N/m^2)
    Ustar=(CD.*(U10.^2)).^0.5; %Wind shear velocity in (m/s)

%Calculating wind drag coefficient based on Wu (1982), method used by ADCIRC
elseif strcmp(CalcMethod,'wu')==1

    CD=(0.8+0.065*U10).*1e-3; %Wind drag coefficient
    CD(U10<7.5)=1.2875.*1e-3; %Wind drag coefficient
    %CD(CD>0.0034)=0.0034; %Wind drag coefficient
    Tauwind=Rhoa.*CD.*U10.^2; %Wind shear stress in (N/m^2)
    Ustar=(CD.*(U10.^2)).^0.5; %Wind shear velocity in (m/s)

%Calculating wind drag coefficient based on Hwang (2011)
elseif strcmp(CalcMethod,'hwang')==1

    CD=(-0.16.*U10.^2+9.67.*U10+80.58).*1e-5; %Wind drag coefficient
    %CD(U10>=50.33)=(2.026./U10(U10>=50.33)).^2; %Upper limitation by Rogers et al. (2012) using U*^2=CD*U10^2 or CD=(U*/U10)^2
    Tauwind=Rhoa.*CD.*U10.^2; %Wind shear stress in (N/m^2)
    Ustar=(CD.*(U10.^2)).^0.5; %Wind shear velocity in (m/s)

%Calculating wind drag coefficient based on Zijlema et al. (2012)
elseif strcmp(CalcMethod,'zijlema')==1

    CD=(0.55+2.97.*(U10./31.5)-1.49.*(U10./31.5).^2).*1e-3; %Wind drag coefficient
    Tauwind=Rhoa.*CD.*U10.^2; %Wind shear stress in (N/m^2)
    Ustar=(CD.*(U10.^2)).^0.5; %Wind shear velocity in (m/s)

%Calculating fully developed condition based Coastal Engineering Manual (CEM), U.S. Army Corps of Engineers (2015)
elseif strcmp(CalcMethod,'cem')==1

    CD=0.001.*(1.1+0.035.*U10); %Wind drag coefficient
    Tauwind=Rhoa.*CD.*U10.^2; %Wind shear stress in (N/m^2)
    Ustar=(CD.*(U10.^2)).^0.5; %Wind shear velocity in (m/s)

end

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1

    %Plotting
    plot(U10,CD)

    xlabel('Wind Velocity')
    ylabel('Wind Drag Coefficient, CD')

end

%--------------------------------------------------------------------------
