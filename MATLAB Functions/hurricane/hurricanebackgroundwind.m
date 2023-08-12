function [Vxgridbackwind, Vygridbackwind, Vgridbackwind, Rgrid] = hurricanebackgroundwind(xgrid, ygrid, Vxgrid, Vygrid, xCenter, yCenter, Vt, VtAzmdir, RVmax, backwindCalcMethod, Rmax, distCalcMethod, dispout)
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

hurricanebackgroundwind
=======================

.. code:: MATLAB

    [Vxgridbackwind, Vygridbackwind, Vgridbackwind, Rgrid] = hurricanebackgroundwind(xgrid, ygrid, Vxgrid, Vygrid, xCenter, yCenter, Vt, VtAzmdir, RVmax, backwindCalcMethod, Rmax, distCalcMethod, dispout)

Description
-----------

Calculate and add background wind velocity due to hurricane front motion to hurricane rotational wind velocity

Inputs
------

xgrid
    | x (longitude) of points which outputs are calculated at as a [M*N] array 
    | xgrid can be a single point or 1d or 2d array 
ygrid
    | y (latitude) of points which outputs are calculated at as a [M*N] array 
    | ygrid can be a single point or 1d or 2d array
Vxgrid
    | Hurricane wind velocity in x (East) direction as a [M*N*L] array in (m/s)
    | L is a number of time steps
    | If only background wind is required, set Vxgrid=0 and Vygrid=0
Vygrid
    | Hurricane wind velocity in y (North) direction as a [M*N*L] array in (m/s)
    | L is a number of time steps
    | If only background wind is required, set Vxgrid=0 and Vygrid=0
xCenter
    x (longitude) of hurricane center (track) as a [L] array
yCenter
    y (latitude) of hurricane center (track) as a [L] array
Vt
    Hurricane central translational velocity as a [L] array in (m/s)
VtAzmdir
    | Hurricane center velocity azimuth (bearing) direction as a [L] array in (Degree)
    | azimuth (bearing) direction which is measured clockwise from the north:
    | 0 (degree): toward North, 90 (degree): toward East, 180 (degree): toward South, 270 (degree): toward West 
RVmax
    Distance (radius) from hurricane center to a location of maximum hurricane wind velocity as a [L] array in (m)
backwindCalcMethod='slosh';
    | Calculation method for adding background wind velocity due to hurricane motion
    | background wind velocity is added to points with R<=Rmax, e.g. Chavas & Emanuel (2010), Lin & Chavas (2012)
    | 'no': background wind velocity is not added to hurricane velocities
    | 'constant': background wind velocity is added as constant value to hurricane velocities
    | 'slosh': Background wind velocity is calculated and added using SLOSH model by Jelesnianski et al. (1992)
    | 'lin': Background wind velocity is calculated and added based on Lin & Chavas (2012)
Rmax=423e3;
    | Maximum radius of hurricane from hurricane center in (m)
    | Outputs for points with R>Rmax is set to zero
    | Median values of Rmax is 423 km (e.g. Chavas & Emanuel 2010; Lin & Chavas, 2012)
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

Vxgridbackwind
    Hurricane wind velocity with background wind in x (East) direction on defined mesh in (m/s)
Vygridbackwind
    Hurricane wind velocity with background wind in y (North) direction on defined mesh in (m/s)
Vgridbackwind
    | Resultant hurricane wind velocity (Vx^2+Vy^2)^0.5 with background wind on defined mesh in (m/s)
    | Gradient wind velocity converted to wind velocity at 10 m above surface by V=VgToVCoeff*Vg
Rgrid
    | Distance (radius) from hurricane center to each point on the grid
    | Note: Outputs has dimension of [M,N,L] where [M,N] is size of the x-y grid and [L] is number of time steps

Examples
--------

.. code:: MATLAB

    %EXAMPLE 1

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

    %Calculating distance and angle using spherical law of cosines
    Rgrid=(acos(sin(deg2rad(latCenter)).*sin(deg2rad(ygrid))+cos(deg2rad(latCenter)).*cos(deg2rad(ygrid)).*cos(deg2rad(xgrid)-deg2rad(longCenter)))).*6371000; %Radius
    thetagrid=atan2(sin(deg2rad(xgrid)-deg2rad(longCenter)).*cos(deg2rad(ygrid)),cos(deg2rad(latCenter)).*sin(deg2rad(ygrid))-sin(deg2rad(latCenter)).*cos(deg2rad(ygrid)).*cos(deg2rad(xgrid)-deg2rad(longCenter))); %Azimuth in radian
    thetagrid=-thetagrid+pi/2; %Converting azimuth to trigonometric direction
    thetagrid=thetagrid+pi/2; %Angle of velocity vector in degree (right angle respect to radius)

    %Calculating hurricane velocity at each radius using SLOSH model
    RVmax=32197; %Radius from hurricane center to a location of maximum hurricane wind
    Vgrid=Vgmax.*(2.*RVmax.*Rgrid)./((RVmax)^2+(Rgrid).^2); %Hurricane wind velocity at radius R
    Vxgrid=Vgrid.*cos(thetagrid); %Hurricane velocity in x (East) direction
    Vygrid=Vgrid.*sin(thetagrid); %Hurricane velocity in y (North) direction

    [Vxgridbackwind,Vygridbackwind,Vgridbackwind,Rgrid]=hurricanebackgroundwind(xgrid,ygrid,Vxgrid,Vygrid,longCenter,latCenter,Vt,VtAzmdir,RVmax,'slosh',423e3,'gc','quiver');


    %EXAMPLE 2

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

    RVmax=32197; %Radius from hurricane center to a location of maximum hurricane wind

    Vxgrid=0; %Hurricane velocity in x (East) direction
    Vygrid=0; %Hurricane velocity in y (North) direction

    [Vxgridbackwind,Vygridbackwind,Vgridbackwind,Rgrid]=hurricanebackgroundwind(xgrid,ygrid,Vxgrid,Vygrid,longCenter,latCenter,Vt,VtAzmdir,RVmax,'slosh',423e3,'gc','quiver');


References
----------

Data

* www.nhc.noaa.gov/data/
* www.nhc.noaa.gov/data/hurdat/hurdat2-format-nencpac.pdf
* coast.noaa.gov/hurricanes
* www.aoml.noaa.gov/hrd/data_sub/re_anal.html

Chavas, D. R., & Emanuel, K. A. (2010). 
A QuikSCAT climatology of tropical cyclone size. 
Geophysical Research Letters, 37(18).

Jelesnianski, C. P., Chen, J., & Shaffer, W. A. (1992). 
SLOSH: Sea, lake, and overland surges from hurricanes (Vol. 48). 
US Department of Commerce, National Oceanic and Atmospheric Administration, National Weather Service.

Lin, N., & Chavas, D. (2012). 
On hurricane parametric wind and applications in storm surge modeling. 
Journal of Geophysical Research: Atmospheres, 117(D9).

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
    case 9
        backwindCalcMethod='slosh'; Rmax=423e3; distCalcMethod='gc'; dispout='no';
    case 10
        Rmax=423e3; distCalcMethod='gc'; dispout='no';
    case 11
        distCalcMethod='gc'; dispout='no';
    case 12
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

if isrow(Vt)==1
    Vt=Vt';
end

if isrow(VtAzmdir)==1
    VtAzmdir=VtAzmdir';
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
    Vgridbackwind=zeros(Mgrid,Ngrid,length(xCenter(:,1))); %Pre-assigning array
    Vxgridbackwind=zeros(Mgrid,Ngrid,length(xCenter(:,1))); %Pre-assigning array
    Vygridbackwind=zeros(Mgrid,Ngrid,length(xCenter(:,1))); %Pre-assigning array
else
    Rgrid=zeros(Mgrid,Ngrid); %Pre-assigning array
    thetagrid=zeros(Mgrid,Ngrid); %Pre-assigning array
    Vgridbackwind=zeros(Mgrid,Ngrid); %Pre-assigning array
    Vxgridbackwind=zeros(Mgrid,Ngrid); %Pre-assigning array
    Vygridbackwind=zeros(Mgrid,Ngrid); %Pre-assigning array
end

VtR=zeros(Mgrid,Ngrid); %Pre-assigning array

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
%Adding background wind velocity due to hurricane motion to hurricane wind velocity

%Converting azimuth (bearing) to trigonometric direction
VtTrigdir=-VtAzmdir+90;

%Add 360 to all numbers to have them all positive
%Use mod(360) to take care of the ones larger than 360
VtTrigdir=mod((VtTrigdir+360),360); 

%Adding background wind velocity as constant value to all hurricane velocities
if strcmp(backwindCalcMethod,'constant')==1

    for i=1:length(xCenter(:,1))
        VtR(:)=Vt(i,1); %Translational wind velocity at radius R
        VtR(Rgrid(:,:,i)>=Rmax)=0; %e.g. Chavas & Emanuel (2010), Lin & Chavas (2012)

        Vxgridbackwind(:,:,i)=Vxgrid(:,:,i)+VtR.*cos(deg2rad(VtTrigdir(i,1))); %Hurricane velocity in x (East) direction
        Vygridbackwind(:,:,i)=Vygrid(:,:,i)+VtR.*sin(deg2rad(VtTrigdir(i,1))); %Hurricane velocity in y (North) direction
    end

    Vgridbackwind=sqrt(Vxgridbackwind.^2+Vygridbackwind.^2); %Assigning hurricane velocity at each track point

%Adding background wind velocity using SLOSH model by Jelesnianski et al. (1992)
elseif strcmp(backwindCalcMethod,'slosh')==1

    for i=1:length(xCenter(:,1))
        VtR=Vt(i,1).*(RVmax(i,1).*Rgrid(:,:,i))./((RVmax(i,1))^2+(Rgrid(:,:,i)).^2); %Translational wind velocity at radius R
        VtR(VtR<0)=0;
        VtR(Rgrid(:,:,i)>=Rmax)=0; %e.g. Chavas & Emanuel (2010), Lin & Chavas (2012)

        Vxgridbackwind(:,:,i)=Vxgrid(:,:,i)+VtR.*cos(deg2rad(VtTrigdir(i,1))); %Hurricane velocity in x (East) direction
        Vygridbackwind(:,:,i)=Vygrid(:,:,i)+VtR.*sin(deg2rad(VtTrigdir(i,1))); %Hurricane velocity in y (North) direction
    end

    Vgridbackwind=sqrt(Vxgridbackwind.^2+Vygridbackwind.^2); %Assigning hurricane velocity at each track point

%Adding background wind velocity using Lin & Chavas (2012)
elseif strcmp(backwindCalcMethod,'lin')==1

    for i=1:length(xCenter(:,1))
        %Hurricane is in northern hemisphere
        if mean(yCenter(:,1))>=0
            %In northern hemisphere, velocity vector rotate counter clockwise with right angle respect to radius

            VtR(:)=0.55*Vt(i,1); %Translational wind velocity at radius R
            VtR(Rgrid(:,:,i)>=Rmax)=0; %e.g. Chavas & Emanuel (2010), Lin & Chavas (2012)

            Vxgridbackwind(:,:,i)=Vxgrid(:,:,i)+VtR.*cos(deg2rad(VtTrigdir(i,1)+20)); %Hurricane velocity in x (East) direction
            Vygridbackwind(:,:,i)=Vygrid(:,:,i)+VtR.*sin(deg2rad(VtTrigdir(i,1)+20)); %Hurricane velocity in y (North) direction

        %Hurricane is in southern hemisphere
        elseif mean(yCenter(:,1))<0
            %In southern hemisphere, velocity vector rotate clockwise with right angle respect to radius

            VtR(:)=0.55*Vt(i,1); %Translational wind velocity at radius R
            VtR(Rgrid(:,:,i)>=Rmax)=0; %e.g. Chavas & Emanuel (2010), Lin & Chavas (2012)

            Vxgridbackwind(:,:,i)=Vxgrid(:,:,i)+VtR.*cos(deg2rad(VtTrigdir(i,1)-20)); %Hurricane velocity in x (East) direction
            Vygridbackwind(:,:,i)=Vygrid(:,:,i)+VtR.*sin(deg2rad(VtTrigdir(i,1)-20)); %Hurricane velocity in y (North) direction

        end
    end

    Vgridbackwind=sqrt(Vxgridbackwind.^2+Vygridbackwind.^2); %Assigning hurricane velocity at each track point

end

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'imagesc')==1 %2 dimensional plot using imagesc or imshow
    
    size(Vgridbackwind)
    imagesc(xgrid,ygrid,Vgridbackwind(:,:,end))
    set(gca,'YDir','normal'); %Correct y axis direction

elseif strcmp(dispout,'pcolor')==1 %2 dimensional plot using pcolor

    pcol=pcolor(xgrid,ygrid,Vgridbackwind(:,:,end));
    set(pcol,'edgecolor','none')

elseif strcmp(dispout,'contour')==1 %2 dimensional contour plot

    [Cont,hCont]=contourf(xgrid,ygrid,Vgridbackwind(:,:,end));
    set(hCont,'LineColor','none');

elseif strcmp(dispout,'quiver')==1 %Surface plot
    
    [Cont,hCont]=contour(xgrid,ygrid,Vgridbackwind(:,:,end));
    hold on
    qvr1=quiver(xgrid,ygrid,Vxgridbackwind(:,:,end),Vygridbackwind(:,:,end));
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
