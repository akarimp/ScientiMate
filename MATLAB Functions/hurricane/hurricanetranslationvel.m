function [Vt, VtAzmdir, VtTrigdir, distxy] = hurricanetranslationvel(xCenter, yCenter, dt, distCalcMethod, CalcMethod, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-10-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

hurricanetranslationvel
=======================

.. code:: MATLAB

    [Vt, VtAzmdir, VtTrigdir, distxy] = hurricanetranslationvel(xCenter, yCenter, dt, distCalcMethod, CalcMethod, dispout)

Description
-----------

Calculate hurricane center translational (forward motion) velocity

Inputs
------

xCenter
    x (longitude) of hurricane center (track)
yCenter
    y (latitude) of hurricane center (track)
dt=6*3600;
    | Time interval between pressure data points in (s)
    | National Hurricane Center reports data every 6 hours 
distCalcMethod='gc';
    | Distance calculation method 
    | 'cart': Distances are calculated on cartesian coordinate
    | 'gc': Distances are calculated on Great Circle based on Vincenty formula, Vincenty (1975)
    | Earth radius coonsidered as mean earth radius=6371000 m
CalcMethod='backward';
    | Calculation method 
    | 'forward': Calculate hurricane central pressure intensity change over time using forward difference method
    |     If CalcMethod='forward'; then last element is zero
    | 'backward': Calculate hurricane central pressure intensity change over time using backward difference method
    |     If CalcMethod='backward'; then first element is zero
    | 'central': Calculate hurricane central pressure intensity change over time using central difference method
    |     If CalcMethod='central'; then first and last elements are zero
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

Vt
    Hurricane central translational velocity in (m/s)
VtAzmdir
    | Hurricane center velocity azimuth (bearing) direction in (Degree)
    | azimuth (bearing) direction which is measured clockwise from the north:
    | 0 (degree): toward North, 90 (degree): toward East, 180 (degree): toward South, 270 (degree): toward West 
VtTrigdir
    Hurricane center velocity trigonometric direction in (Degree)
distxy
    Distance between hurricane locations of hurricane center in (m)

Examples
--------

.. code:: MATLAB

    %Longitude of Hurricane Katrine best track
    longtrack=[-75.1;-75.7;-76.2;-76.5;-76.9;-77.7;-78.4;-79.0;-79.6;-80.1;-80.3;-81.3;...
        -82.0;-82.6;-83.3;-84.0;-84.7;-85.3;-85.9;-86.7;-87.7;-88.6;-89.2;-89.6;...
        -89.6;-89.6;-89.6;-89.6;-89.1;-88.6;-88.0;-87.0;-85.3;-82.9];

    %Latitude of Hurricane Katrine best track
    lattrack=[23.1;23.4;23.8;24.5;25.4;26.0;26.1;26.2;26.2;26.0;25.9;25.4;...
        25.1;24.9;24.6;24.4;24.4;24.5;24.8;25.2;25.7;26.3;27.2;28.2;...
        29.3;29.5;30.2;31.1;32.6;34.1;35.6;37.0;38.6;40.1];

    [Vt,VtAzmdir,VtTrigdir,distxy]=hurricanetranslationvel(longtrack,lattrack,6*3600,'gc','backward','yes');

References
----------

Data

* www.nhc.noaa.gov/data/
* www.nhc.noaa.gov/data/hurdat/hurdat2-format-nencpac.pdf
* coast.noaa.gov/hurricanes
* www.aoml.noaa.gov/hrd/data_sub/re_anal.html

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
        dt=6*3600; distCalcMethod='gc'; CalcMethod='backward'; dispout='no';
    case 3
        distCalcMethod='gc'; CalcMethod='backward'; dispout='no';
    case 4
        CalcMethod='backward'; dispout='no';
    case 5
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

%--------------------------------------------------------------------------
%Pre-assigning array

distxy=zeros(length(xCenter(:,1)),1); %Pre-assigning array
thetarad=zeros(length(xCenter(:,1)),1); %Pre-assigning array
deltasigma=zeros(length(xCenter(:,1)),1); %Pre-assigning array
azimuthrad=zeros(length(xCenter(:,1)),1); %Pre-assigning array

%--------------------------------------------------------------------------
%Calculating distance between hurricane locations

%Calculating distance using cartesian formula
if strcmp(distCalcMethod,'cart')==1

    %Calculating distance between hurricane locations using forward difference method
    if strcmp(CalcMethod,'forward')==1

        for i=1:length(xCenter(:,1))-1

            %Calculation range: (1:end-1,1), last element is zero
            distxy(i,1)=sqrt((xCenter(i+1,1)-xCenter(i,1)).^2+(yCenter(i+1,1)-yCenter(i,1)).^2); %Calculating distance from (x,y) to (x(1),y(1))

            %Calculating angle of the line between start and end points
            thetarad(i,1)=atan2(yCenter(i+1,1)-yCenter(i,1),xCenter(i+1,1)-xCenter(i,1)); %Angle in radian

        end

    %Calculating distance between hurricane locations using backward difference method
    elseif strcmp(CalcMethod,'backward')==1

        for i=2:length(xCenter(:,1))

            %Calculation range: (2:end,1), first element is zero
            distxy(i,1)=sqrt((xCenter(i,1)-xCenter(i-1,1)).^2+(yCenter(i,1)-yCenter(i-1,1)).^2); %Calculating distance from (x,y) to (x(1),y(1))

            %Calculating angle of the line between start and end points
            thetarad(i,1)=atan2(yCenter(i,1)-yCenter(i-1,1),xCenter(i,1)-xCenter(i-1,1)); %Angle in radian
    
        end

    %Calculating distance between hurricane locations using centra difference method
    elseif strcmp(CalcMethod,'central')==1

        for i=2:length(xCenter(:,1))-1

            %Calculation range: (2:end-1,1), first and last elements are zero
            distxy(i,1)=sqrt((xCenter(i+1,1)-xCenter(i-1,1)).^2+(yCenter(i+1,1)-yCenter(i-1,1)).^2); %Calculating distance from (x,y) to (x(1),y(1))

            %Calculating angle of the line between start and end points
            thetarad(i,1)=atan2(yCenter(i+1,1)-yCenter(i-1,1),xCenter(i+1,1)-xCenter(i-1,1)); %Angle in radian

        end

    end
    
    theta=rad2deg(thetarad); %Angle in degree

    %Add 360 to all numbers to have them all positive
    %Use mod(360) to take care of the ones larger than 360
    theta=mod((theta+360),360); 


%Calculating distance using Vincenty formula
elseif strcmp(distCalcMethod,'gc')==1

    %Calculating distance between hurricane locations using forward difference method
    if strcmp(CalcMethod,'forward')==1

        for i=1:length(xCenter(:,1))-1

            %Calculation range: (1:end-1,1), last element is zero

            %Converting to radian
            lat1rad=deg2rad(yCenter(i,1));
            lon1rad=deg2rad(xCenter(i,1));
            lat2rad=deg2rad(yCenter(i+1,1));
            lon2rad=deg2rad(xCenter(i+1,1));

            deltalatrad21=lat2rad-lat1rad;
            deltalonrad21=lon2rad-lon1rad;

            deltasigma(i,1)=atan2(sqrt((cos(lat2rad).*sin(deltalonrad21)).^2+(cos(lat1rad).*sin(lat2rad)-sin(lat1rad).*cos(lat2rad).*cos(deltalonrad21)).^2),sin(lat1rad).*sin(lat2rad)+cos(lat1rad).*cos(lat2rad).*cos(deltalonrad21)); %Central angle

            %Calculating azimuth (bearing) between start and end of the line
            azimuthrad(i,1)=atan2(sin(deltalonrad21).*cos(lat2rad),cos(lat1rad).*sin(lat2rad)-sin(lat1rad).*cos(lat2rad).*cos(deltalonrad21)); %Azimuth (bearing) in radian

        end

    %Calculating distance between hurricane locations using backward difference method
    elseif strcmp(CalcMethod,'backward')==1

        for i=2:length(xCenter(:,1))

            %Calculation range: (2:end,1), first element is zero

            %Converting to radian
            lat1rad=deg2rad(yCenter(i-1,1));
            lon1rad=deg2rad(xCenter(i-1,1));
            lat2rad=deg2rad(yCenter(i,1));
            lon2rad=deg2rad(xCenter(i,1));

            deltalatrad21=lat2rad-lat1rad;
            deltalonrad21=lon2rad-lon1rad;

            deltasigma(i,1)=atan2(sqrt((cos(lat2rad).*sin(deltalonrad21)).^2+(cos(lat1rad).*sin(lat2rad)-sin(lat1rad).*cos(lat2rad).*cos(deltalonrad21)).^2),sin(lat1rad).*sin(lat2rad)+cos(lat1rad).*cos(lat2rad).*cos(deltalonrad21)); %Central angle

            %Calculating azimuth (bearing) between start and end of the line
            azimuthrad(i,1)=atan2(sin(deltalonrad21).*cos(lat2rad),cos(lat1rad).*sin(lat2rad)-sin(lat1rad).*cos(lat2rad).*cos(deltalonrad21)); %Azimuth (bearing) in radian

        end

    %Calculating distance between hurricane locations using centra difference method
    elseif strcmp(CalcMethod,'central')==1

        for i=2:length(xCenter(:,1))-1

            %Calculation range: (2:end-1,1), first and last elements are zero

            %Converting to radian
            lat1rad=deg2rad(yCenter(i-1,1));
            lon1rad=deg2rad(xCenter(i-1,1));
            lat2rad=deg2rad(yCenter(i+1,1));
            lon2rad=deg2rad(xCenter(i+1,1));

            deltalatrad21=lat2rad-lat1rad;
            deltalonrad21=lon2rad-lon1rad;

            deltasigma(i,1)=atan2(sqrt((cos(lat2rad).*sin(deltalonrad21)).^2+(cos(lat1rad).*sin(lat2rad)-sin(lat1rad).*cos(lat2rad).*cos(deltalonrad21)).^2),sin(lat1rad).*sin(lat2rad)+cos(lat1rad).*cos(lat2rad).*cos(deltalonrad21)); %Central angle

            %Calculating azimuth (bearing) between start and end of the line
            azimuthrad(i,1)=atan2(sin(deltalonrad21).*cos(lat2rad),cos(lat1rad).*sin(lat2rad)-sin(lat1rad).*cos(lat2rad).*cos(deltalonrad21)); %Azimuth (bearing) in radian

        end

    end

    REarth=6371000; %Earth radius in (m), mean earth radius=6371000 m
    arclen=REarth.*deltasigma; %Total distance of the line 
    distxy=arclen;

    azimuth=rad2deg(azimuthrad); %Azimuth (bearing) in degree

    %Add 360 to all numbers to have them all positive
    %Use mod(360) to take care of the ones larger than 360
    azimuth=mod((azimuth+360),360); 

end

%--------------------------------------------------------------------------
%Calculating hurricane central translational velocity

%Calculating hurricane central translational velocity using forward difference method
if strcmp(CalcMethod,'forward')==1

    %Calculation range: (1:end-1,1), last element is zero
    Vt=distxy./dt; %Central translational velocity

%Calculating hurricane central translational velocity using backward difference method
elseif strcmp(CalcMethod,'backward')==1

    %Calculation range: (2:end,1), first element is zero
    Vt=distxy./dt; %Central translational velocity

%Calculating hurricane central translational velocity using centra difference method
elseif strcmp(CalcMethod,'central')==1

    %Calculation range: (2:end-1,1), first and last elements are zero
    Vt=distxy./(2*dt); %Central translational velocity

end

%--------------------------------------------------------------------------
%Calculating hurricane central translational direction

%Calculating direction using cartesian formula
if strcmp(distCalcMethod,'cart')==1

    %Converting trigonometric direction to azimuth (bearing)
    VtAzmdir=90-theta; %Hurricane center velocity azimuth (bearing) direction in (Degree)

    VtTrigdir=theta; %Hurricane center velocity trigonometric direction in (Degree)

%Calculating direction using Vincenty formula
elseif strcmp(distCalcMethod,'gc')==1

    VtAzmdir=azimuth; %Hurricane center velocity azimuth (bearing) direction in (Degree)

    %Converting azimuth (bearing) to trigonometric direction
    VtTrigdir=-azimuth+90; %Hurricane center velocity trigonometric direction in (Degree)

end

%Add 360 to all numbers to have them all positive
%Use mod(360) to take care of the ones larger than 360
VtAzmdir=mod((VtAzmdir+360),360); 
VtTrigdir=mod((VtTrigdir+360),360); 

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1
    
    %Plotting data
    t(:,1)=[0:1:length(xCenter(:,1))-1].*dt./3600;

    subplot(2,1,1)
    plot(t,Vt)
    
    xlabel('Time (hr)')
    ylabel('Translational Velocity (m/s)')
    
    subplot(2,1,2)
    plot(t,VtAzmdir)
    
    xlabel('Time (hr)')
    ylabel('Velocity Azimuth (Degree)')

end

%--------------------------------------------------------------------------
