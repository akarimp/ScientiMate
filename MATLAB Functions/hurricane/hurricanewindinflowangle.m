function [Vxgrid, Vygrid, thetaVgrid, thetaVgridtan, thetagrid, Rgrid] = hurricanewindinflowangle(xgrid, ygrid, Vgrid, xCenter, yCenter, RVmax, inflowdirCalcMethod, distCalcMethod, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-11-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

hurricanewindinflowangle
========================

.. code:: MATLAB

    [Vxgrid, Vygrid, thetaVgrid, thetaVgridtan, thetagrid, Rgrid] = hurricanewindinflowangle(xgrid, ygrid, Vgrid, xCenter, yCenter, RVmax, inflowdirCalcMethod, distCalcMethod, dispout)

Description
-----------

Calculate hurricane velocity tangential and inflow angle and inflow velocity in x (East) and y (North) directions

Inputs
------

xgrid
    | x (longitude) of points which outputs are calculated at as a [M*N] array 
    | xgrid can be a single point or 1d or 2d array 
ygrid
    | y (latitude) of points which outputs are calculated at as a [M*N] array 
    | ygrid can be a single point or 1d or 2d array
Vgrid
    | Resultant hurricane wind velocity (Vx^2+Vy^2)^0.5 on (xgrid,ygrid) as a [M*N*L] array in (m/s)
    | L is a number of time steps
    | If only angle values are required, then set Vgrid equal to an arbitary constant such as Vgrid=1
    | For demonstaration purpose, set Vgrid equal to an arbitary constant such as Vgrid=1
xCenter
    x (longitude) of hurricane center (track) as a [L] array
yCenter
    y (latitude) of hurricane center (track) as a [L] array
RVmax
    Distance (radius) from hurricane center to a location of maximum hurricane wind velocity as a [L] array in (m)
inflowdirCalcMethod='bretschneider';
    | Inflow angle calculation method 
    | 'no': Inflow angle are not calculated, thetaVgrid=thetaVgridtan
    | 'bretschneider': Inflow angle are calculated based on Bretschneider (1972)
    | 'sobey': Inflow angle are calculated based on Sobey et al. (1977)
distCalcMethod='gc';
    | Distance calculation method 
    | 'cart': Distances are calculated on cartesian coordinate
    | 'gc': Distances are calculated on Great Circle based on Vincenty formula, Vincenty (1975)
    | Earth radius coonsidered as mean earth radius=6371000 m
dispout='no';
    | Define to display outputs or not
    | 'imagesc': 2 dimensional plot using imagesc or imshow
    | 'pcolor': 2 dimensional plot using pcolor
    | 'contour': 2 dimensional contour plot, number of contour=ncolor
    | 'quiver': 2 dimensional vector plot 
    | 'no': not display 
    | Use dispout='no' if calculation mesh is not 2d array
    | if there is more than one time step, only the last one is plotted
    | if flattendata='yes'; then dispout is set as dispout='no';

Outputs
-------

Vxgrid
    Hurricane wind velocity after applying inflow angle in x (East) direction on defined mesh in (m/s)
Vygrid
    Hurricane wind velocity after applying inflow angle in y (North) direction on defined mesh in (m/s)
thetaVgrid
    | Inflow angle (trigonometric direction) of hurricane velocity at each grid point in (Degree)
    | Inflow angle: angle between the inwardly spiraling surface wind 
    |     and the circular isobars around the hurricane center (Boose et al., 2004)
thetaVgridtan
    | Angle (trigonometric direction) of hurricane velocity at each grid point in (Degree)
    | thetaVgridtan is tangential angle respect to radius. 
    | Note: Outputs has dimension of [M,N,L] where [M,N] is size of the x-y grid and [L] is number of time steps
thetagrid
    Angle from hurricane center to each point on the grid in (Degree)
Rgrid
    Distance (radius) from hurricane center to each point on the grid

Examples
--------

.. code:: MATLAB

    %Creating calculation mesh
    [xgrid,ygrid]=meshgrid(linspace(-98,-68,100),linspace(16,44,100));

    %Longitude of Hurricane Katrine center at max velocity
    longCenter=-88.6;

    %Latitude of Hurricane Katrine center at max velocity
    latCenter=26.3;

    %Hurricane Katrina translational velocity (m/s) at max velocity
    Vt=5.18467;

    %Hurricane Katrina velocity azimuth (bearing) in (Degree) at max velocity
    VtAzmdir=306.76219;

    %Hurricane Katrina 1-min sustained maximum velocity (m/s) at max velocity
    Vmax=76.5;
    Vmax=Vmax-Vt; %Removing hurricane translation velocity from Vmax
    Vgmax=Vmax/0.8; %Converting surface velocity to gradient velocity

    %Calculating distance using spherical law of cosines
    Rgrid=(acos(sin(deg2rad(latCenter)).*sin(deg2rad(ygrid))+cos(deg2rad(latCenter)).*cos(deg2rad(ygrid)).*cos(deg2rad(xgrid)-deg2rad(longCenter)))).*6371000; %Radius

    %Calculating hurricane velocity at each radius using SLOSH model
    RVmax=32197; %Radius from hurricane center to a location of maximum hurricane wind
    Vgrid=Vgmax.*(2.*RVmax.*Rgrid)./((RVmax)^2+(Rgrid).^2); %Hurricane wind velocity at radius R

    [Vxgrid,Vygrid,thetaVgrid,thetaVgridtan,thetagrid,Rgrid]=hurricanewindinflowangle(xgrid,ygrid,Vgrid,longCenter,latCenter,RVmax,'bretschneider','gc','quiver');

References
----------

Data

* www.nhc.noaa.gov/data/
* www.nhc.noaa.gov/data/hurdat/hurdat2-format-nencpac.pdf
* coast.noaa.gov/hurricanes
* www.aoml.noaa.gov/hrd/data_sub/re_anal.html

Boose, E. R., Serrano, M. I., & Foster, D. R. (2004). 
Landscape and regional impacts of hurricanes in Puerto Rico. 
Ecological Monographs, 74(2), 335-352.

Bretschneider, C. L. (1972, January). 
A non-dimensional stationary hurricane wave model. 
In Offshore Technology Conference. Offshore Technology Conference.

Phadke, A. C., Martino, C. D., Cheung, K. F., & Houston, S. H. (2003). 
Modeling of tropical cyclone winds and waves for emergency management. 
Ocean Engineering, 30(4), 553-578.

Sobey, R. J., Harper, B. A., & Stark, K. P. (1977). 
Numerical simulation of tropical cyclone storm surge. 
James Cook University of North Queensland, Department of Civil & Systems Engineering.

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
    case 6
        inflowdirCalcMethod='bretschneider'; distCalcMethod='gc'; dispout='no';
    case 7
        distCalcMethod='gc'; dispout='no';
    case 8
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(xCenter)==1
    xCenter=xCenter';
end

if isrow(yCenter)==1
    yCenter=yCenter';
end

if isrow(RVmax)==1
    RVmax=RVmax';
end

%--------------------------------------------------------------------------
%Pre-assigning array

[Mgrid,Ngrid]=size(xgrid);
if length(xCenter(:,1))>1
    Rgrid=zeros(Mgrid,Ngrid,length(xCenter(:,1))); %Pre-assigning array
    thetagrid=zeros(Mgrid,Ngrid,length(xCenter(:,1))); %Pre-assigning array
    thetaVgridtan=zeros(Mgrid,Ngrid,length(xCenter(:,1))); %Pre-assigning array
    thetaVgrid=zeros(Mgrid,Ngrid,length(xCenter(:,1))); %Pre-assigning array
    Beta=zeros(Mgrid,Ngrid,length(xCenter(:,1))); %Pre-assigning array
else
    Rgrid=zeros(Mgrid,Ngrid); %Pre-assigning array
    thetagrid=zeros(Mgrid,Ngrid); %Pre-assigning array
    thetaVgridtan=zeros(Mgrid,Ngrid); %Pre-assigning array
    thetaVgrid=zeros(Mgrid,Ngrid); %Pre-assigning array
    Beta=zeros(Mgrid,Ngrid); %Pre-assigning array
end

%--------------------------------------------------------------------------
%Calculating distance (radius) from hurricane center to each point

for i=1:length(xCenter(:,1))

    %Calculating distance using cartesian formula
    if strcmp(distCalcMethod,'cart')==1

        distxy=sqrt((xgrid-xCenter(i,1)).^2+(ygrid-yCenter(i,1)).^2); %Calculating distance from (x,y) to (x(1),y(1))

        %Calculating angle of the line between start and end points

        thetarad=atan2(ygrid-yCenter(i,1),xgrid-xCenter(i,1)); %Angle in radian
        theta=rad2deg(thetarad); %Angle in degree

        %Add 360 to all numbers to have them all positive
        %Use mod(360) to take care of the ones larger than 360
        theta=mod((theta+360),360); 

        Rgrid(:,:,i)=distxy; %Distance from hurricane center to each point in (m)
        thetagrid(:,:,i)=theta; %Angle from hurricane center to each point in (degree)

    %Calculating distance using Vincenty formula
    elseif strcmp(distCalcMethod,'gc')==1

        %Converting to radian
        lat1rad=deg2rad(yCenter(i,1));
        lon1rad=deg2rad(xCenter(i,1));
        lat2rad=deg2rad(ygrid);
        lon2rad=deg2rad(xgrid);

        deltalatrad21=lat2rad-lat1rad;
        deltalonrad21=lon2rad-lon1rad;

        REarth=6371000; %Earth radius in (m), mean earth radius=6371000 m
        deltasigma=atan2(sqrt((cos(lat2rad).*sin(deltalonrad21)).^2+(cos(lat1rad).*sin(lat2rad)-sin(lat1rad).*cos(lat2rad).*cos(deltalonrad21)).^2),sin(lat1rad).*sin(lat2rad)+cos(lat1rad).*cos(lat2rad).*cos(deltalonrad21)); %Central angle
        arclen=REarth.*deltasigma; %Total distance of the line 
        distxy=arclen;

        %Calculating azimuth (bearing) between start and end of the line

        azimuthrad=atan2(sin(deltalonrad21).*cos(lat2rad),cos(lat1rad).*sin(lat2rad)-sin(lat1rad).*cos(lat2rad).*cos(deltalonrad21)); %Azimuth (bearing) in radian
        azimuth=rad2deg(azimuthrad); %Azimuth (bearing) in degree

        %Add 360 to all numbers to have them all positive
        %Use mod(360) to take care of the ones larger than 360
        azimuth=mod((azimuth+360),360); 

        Rgrid(:,:,i)=distxy; %Distance from hurricane center to each point in (m)
        thetagrid(:,:,i)=azimuth; %Angle from hurricane center to each point in (degree)

    end

end

%--------------------------------------------------------------------------
%Calculating hurricane velocity tangential and inflow angle and inflow velocity in x (East) and y (North) directions
%Tangential velocity has a right angle respect to radius
%thetaVgridtan is tangential angle respect to radius 
%thetaVgrid is inflow angle calculated based on Bretschneider (1972)
%Inflow angle: angle between the inwardly spiraling surface wind 
%and the circular isobars around the hurricane center (Boose et al., 2004)

for i=1:length(xCenter(:,1))

    %Calculating inflow angle based on Bretschneider (1972), e.g Phadke et al.,(2003), 
    if strcmp(inflowdirCalcMethod,'bretschneider')==1
        %Jeong, C. K., Panchang, V., & Demirbilek, Z. (2012). 
        %Parametric adjustments to the rankine vortex wind model for Gulf of Mexico hurricanes. 
        %Journal of Offshore Mechanics and Arctic Engineering, 134(4), 041102.
        
        %Beta is an inflow angle inward from the tangential flow (Martino et al., 2001)
        Beta(Rgrid<RVmax(i,1))=10+(1+Rgrid(Rgrid<RVmax(i,1))./RVmax(i,1));
        Beta(Rgrid>=RVmax(i,1) & Rgrid<1.2*RVmax(i,1))=20+25.*(Rgrid(Rgrid>=RVmax(i,1) & Rgrid<1.2*RVmax(i,1))./RVmax(i,1)-1);
        Beta(Rgrid>=1.2*RVmax(i,1))=25;

    %Calculating inflow angle based on Sobey et al. (1977), e.g. Martino et al. (2001)
    elseif strcmp(inflowdirCalcMethod,'sobey')==1
        %Sobey, R. J., Harper, B. A., & Stark, K. P. (1977). 
        %Numerical simulation of tropical cyclone storm surge. 
        %James Cook University of North Queensland, Department of Civil & Systems Engineering.

        %Martino, C. D., Cheung, K. F., Phadke, A. C., & Houston, S. H. (2001, January). 
        %Modeling of hurricane waves in Hawaiian waters. 
        %In The Eleventh International Offshore and Polar Engineering Conference. International Society of Offshore and Polar Engineers.    
        
        %Beta is an inflow angle inward from the tangential flow (Martino et al., 2001)
        Beta(Rgrid<RVmax(i,1))=10.*(Rgrid(Rgrid<RVmax(i,1))./RVmax(i,1));
        Beta(Rgrid>=RVmax(i,1) & Rgrid<1.2*RVmax(i,1))=10+75.*(Rgrid(Rgrid>=RVmax(i,1) & Rgrid<1.2*RVmax(i,1))./RVmax(i,1)-1);
        Beta(Rgrid>=1.2*RVmax(i,1))=25;

    else
        Beta=0;
    
    end

end

%Calculating velocity vector direction using cartesian formula
if strcmp(distCalcMethod,'cart')==1

    %Hurricane is in northern hemisphere
     if mean(yCenter(:,1))>=0
        %In northern hemisphere, velocity vector rotate counter clockwise with right angle respect to radius
        thetaVgridtan=thetagrid+90; %Angle of velocity vector in degree
        thetaVgrid=thetaVgridtan+Beta; %Inflow angle of velocity vector in degree

    %Hurricane is in southern hemisphere
    elseif mean(yCenter(:,1))<0
        %In southern hemisphere, velocity vector rotate clockwise with right angle respect to radius
        thetaVgridtan=thetagrid-90; %Angle of velocity vector in degree
        thetaVgrid=thetaVgridtan-Beta; %Inflow angle of velocity vector in degree

    end

    %Add 360 to all numbers to have them all positive
    %Use mod(360) to take care of the ones larger than 360
    thetaVgridtan=mod((thetaVgridtan+360),360); 
    thetaVgrid=mod((thetaVgrid+360),360); 

%Calculating velocity vector direction using great circle formula
elseif strcmp(distCalcMethod,'gc')==1

    %Converting azimuth (bearing) to trigonometric direction
    thetagrid=-thetagrid+90;

    %Hurricane is in northern hemisphere
     if mean(yCenter(:,1))>=0
        %In northern hemisphere, velocity vector rotate counter clockwise with right angle respect to radius
        thetaVgridtan=thetagrid+90; %Angle of velocity vector in degree
        thetaVgrid=thetaVgridtan+Beta; %Inflow angle of velocity vector in degree

    %Hurricane is in southern hemisphere
    elseif mean(yCenter(:,1))<0
        %In southern hemisphere, velocity vector rotate clockwise with right angle respect to radius
        thetaVgridtan=thetagrid-90; %Angle of velocity vector in degree
        thetaVgrid=thetaVgridtan-Beta; %Inflow angle of velocity vector in degree

    end

    %Add 360 to all numbers to have them all positive
    %Use mod(360) to take care of the ones larger than 360
    thetaVgridtan=mod((thetaVgridtan+360),360); 
    thetaVgrid=mod((thetaVgrid+360),360); 

end

inflowangle=Beta; %Inflow angle of velocity vector in degree
Vxgrid=Vgrid.*cos(deg2rad(thetaVgrid)); %Hurricane surface velocity in x (East) direction
Vygrid=Vgrid.*sin(deg2rad(thetaVgrid)); %Hurricane surface velocity in y (North) direction

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'imagesc')==1 %2 dimensional plot using imagesc or imshow
    
    imagesc(xgrid,ygrid,Vgrid(:,:,end))
    set(gca,'YDir','normal'); %Correct y axis direction

elseif strcmp(dispout,'pcolor')==1 %2 dimensional plot using pcolor

    pcol=pcolor(xgrid,ygrid,Vgrid(:,:,end));
    set(pcol,'edgecolor','none')

elseif strcmp(dispout,'contour')==1 %2 dimensional contour plot

    [Cont,hCont]=contourf(xgrid,ygrid,Vgrid(:,:,end));
    set(hCont,'LineColor','none');

elseif strcmp(dispout,'quiver')==1 %Surface plot
    
    [Cont,hCont]=contour(xgrid,ygrid,Vgrid(:,:,end));
    hold on
    qvr1=quiver(xgrid,ygrid,Vxgrid(:,:,end),Vygrid(:,:,end));
    %set(qvr1,'FaceColor','interp','EdgeColor','none')

end

%Setting plot properties
if (strcmp(dispout,'imagesc')==1 | strcmp(dispout,'pcolor')==1 | strcmp(dispout,'contour')==1 | strcmp(dispout,'quiver')==1)

    %Setting axis limits
    %xlim([min(xgrid(:)),max(xgrid(:))])
    %ylim([min(ygrid(:)),max(ygrid(:))])

    %Plotting colorbar
    colorbar;

    %Setting label
    xlabel('x')
    ylabel('y')
    title('Velocity (m/s)')
    
end

%--------------------------------------------------------------------------
