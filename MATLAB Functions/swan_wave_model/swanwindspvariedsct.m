function [swanwind, windvelpoint, winddirpoint] = swanwindspvariedsct(xgrid, ygrid, xpoint, ypoint, windvel, winddir, winddirtype, windvelmin, savedata, outfilename, outfilelocation, CalcMethod)
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

swanwindspvariedsct
===================

.. code:: MATLAB

    [swanwind, windvelpoint, winddirpoint] = swanwindspvariedsct(xgrid, ygrid, xpoint, ypoint, windvel, winddir, winddirtype, windvelmin, savedata, outfilename, outfilelocation, CalcMethod)

Description
-----------

Generate SWAN wind file for spatially varied wind from scattered input data

Inputs
------

xgrid
    x (longitude) of output grid points as a [M*N] array
ygrid
    y (latitude) of output grid points as a [M*N] array
xpoint
    | x (longitude) of the location (such as meteorological station) that wind is known in that location
    | as a [L] array
    | 1st element is 1st point, 2nd element is 2nd point, 3rd element is 3rd point, ...
ypoint
    | y (latitude) of the location (such as meteorological station) that wind is known in that location
    | as a [L] array
    | 1st element is 1st point, 2nd element is 2nd point, 3rd element is 3rd point, ...
windvel
    | Wind velocity as a [K*L] array
    | K is number of time steps for a time series
    |     1st row (i.e. [1,:]) is 1st time step, 2nd row (i.e. [2,:]) is 2nd time step, ...
    |     Example: windvel(1,:) is the wind velocity data for all points at the first time step
    | L is number of points (such as meteorological stations) that wind velocity is known in those locations
    |     L should be L>=3
    |     1st column (i.e. [:,1]) is 1st point, 2nd column (i.e. [:,2]) is 2nd point, 3rd column (i.e. [:,3]) is 3rd point, ...
    |     Example: windvel(:,1) is the wind velocity data for all time steps at the first point
winddir
    | Wind direction as a [K*L] array in (Degree)
    | K is number of time steps for a time series
    |     1st row (i.e. [1,:]) is 1st time step, 2nd row (i.e. [2,:]) is 2nd time step, ...
    |     Example: windvel(1,:) is the wind direction data for all points at the first time step
    | L is number of points (such as meteorological stations) that wind direction is known in those locations
    |     L should be L>=3
    |     1st column (i.e. [:,1]) is 1st point, 2nd column (i.e. [:,2]) is 2nd point, 3rd column (i.e. [:,3]) is 3rd point, ...
    |     Example: windvel(:,1) is the wind direction data for all time steps at the first point
winddirtype='mete';
    | Define wind direction type
    | 'mete': meteorological wind direction 
    |     Meteorological direction represents a direction wind comes from and is measured counter-clockwise from the North
    |     0 (degree): from North, 90 (degree): from East, 180 (degree): from South, 270 (degree): from West
    | 'trig': trigonometric wind direction
windvelmin=0;
    Minimum allowed wind velocity
savedata='no';
    | Define if save data in a file or not
    | 'no': does not save 
    | 'yes': save data as ascii file
outfilename='swanwind.wnd';
    | Name of output file between ' ' mark, example: 'swanwind.wnd'
    | outfilename should have '.wnd' extension
outfilelocation=pwd;
    Location of output file between ' ' mark, example: 'C:\' in MATLAB, or 'C:/' in Python
CalcMethod='linear';
    | Interpolation method 
    | 'linear': Use default or 'linear' method to interpolate
    | 'nearest': Use nearest neighbor method to interpolate

Outputs
-------

swanwind
    | Spatially varied wind velocity data formated for SWAN
    | Wind velocity data at each time step is assigned into the grid points
windvelpoint
    | Nearest interpolated wind velocity at (xpoint,ypoint) as a [K*L] array
    | K is number of time steps for a time series
    | L is number of points (such as meteorological stations) that wind velocity is known in those locations
winddirpoint
    | Nearest interpolated wind direction at (xpoint,ypoint) as a [K*L] array
    | K is number of time steps for a time series
    | L is number of points (such as meteorological stations) that wind direction is known in those locations

Examples
--------

.. code:: MATLAB

    [xgrid,ygrid]=meshgrid(linspace(-91,-90,100),linspace(28,30,100));
    windvel=[[10.5,10.55,10.6];[10.64,10.69,10.74];[10.7,10.75,10.8];[10.4,10.45,10.5]]; %Data for 4 time steps
    winddir=[[50,55,60];[64,69,74];[70,75,80];[40,45,50]]; %Data for 4 time steps
    xpoint=[-90.5;-90.3;-90.7]; %Data are known at 3 locations
    ypoint=[29.2;29;28.8]; %Data are known at 3 locations
    winddirtype='mete';
    windvelmin=2.5;
    savedata='no';
    outfilename='swanwind.wnd';
    outfilelocation=pwd;
    CalcMethod='linear';
    [swanwind,windvelpoint,winddirpoint]=swanwindspvariedsct(xgrid,ygrid,xpoint,ypoint,windvel,winddir,winddirtype,windvelmin,savedata,outfilename,outfilelocation,CalcMethod);


    [xgrid,ygrid]=meshgrid(linspace(-91,-90,100),linspace(28,30,100));
    windvel=[10.5,10.55,10.6]; %Data for 1 time step
    winddir=[50,55,60]; %Data for 1 time step
    xpoint=[-90.5;-90.3;-90.7]; %Data are known at 3 locations
    ypoint=[29.2;29;28.8]; %Data are known at 3 locations
    winddirtype='mete';
    windvelmin=2.5;
    savedata='no';
    outfilename='swanwind.wnd';
    outfilelocation=pwd;
    CalcMethod='linear';
    [swanwind,windvelpoint,winddirpoint]=swanwindspvariedsct(xgrid,ygrid,xpoint,ypoint,windvel,winddir,winddirtype,windvelmin,savedata,outfilename,outfilelocation,CalcMethod);

References
----------

Booij, N. R. R. C., Ris, R. C., & Holthuijsen, L. H. (1999). 
A third‐generation wave model for coastal regions: 1. Model description and validation. 
Journal of geophysical research: Oceans, 104(C4), 7649-7666.

SWAN Team. (2007). S
WAN user manual. 
Delft University of Technology. The Netherlands.

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
        winddirtype='mete'; windvelmin=0; savedata='no'; outfilename='swanwind.wnd'; outfilelocation=pwd; CalcMethod='linear';
    case 7
        windvelmin=0; savedata='no'; outfilename='swanwind.wnd'; outfilelocation=pwd; CalcMethod='linear';
    case 8
        savedata='no'; outfilename='swanwind.wnd'; outfilelocation=pwd; CalcMethod='linear';
    case 9
        outfilename='swanwind.wnd'; outfilelocation=pwd; CalcMethod='linear';
    case 10
        outfilelocation=pwd; CalcMethod='linear';
    case 11
        CalcMethod='linear';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(xpoint)==1
    xpoint=xpoint';
end

if isrow(ypoint)==1
    ypoint=ypoint';
end

%--------------------------------------------------------------------------
%Generating SWAN wind file (.wnd)

%If 
%xgrid has a format of:
% x1 x2 x3 x4
% x1 x2 x3 x4
% x1 x2 x3 x4
% x1 x2 x3 x4
%ygrid gas a format of:
% y1 y1 y1 y1
% y2 y2 y2 y2
% y3 y3 y3 y3
% y4 y4 y4 y4
%Then SWAN wind file (.wnd) for wind with three time steps has a format of (RowColumn,TimeStep):
% U(11,1) U(12,1) U(13,1) U(14,1)
% U(21,1) U(22,1) U(23,1) U(24,1)
% U(31,1) U(32,1) U(33,1) U(34,1)
% U(41,1) U(42,1) U(43,1) U(44,1)
% V(11,1) V(12,1) V(13,1) V(14,1)
% V(21,1) V(22,1) V(23,1) V(24,1)
% V(31,1) V(32,1) V(33,1) V(34,1)
% V(41,1) V(42,1) V(43,1) V(44,1)
% U(11,2) U(12,2) U(13,2) U(14,2)
% U(21,2) U(22,2) U(23,2) U(24,2)
% U(31,2) U(32,2) U(33,2) U(34,2)
% U(41,2) U(42,2) U(43,2) U(44,2)
% V(11,2) V(12,2) V(13,2) V(14,2)
% V(21,2) V(22,2) V(23,2) V(24,2)
% V(31,2) V(32,2) V(33,2) V(34,2)
% V(41,2) V(42,2) V(43,2) V(44,2)
% U(11,3) U(12,3) U(13,3) U(14,3)
% U(21,3) U(22,3) U(23,3) U(24,3)
% U(31,3) U(32,3) U(33,3) U(34,3)
% U(41,3) U(42,3) U(43,3) U(44,3)
% V(11,3) V(12,3) V(13,3) V(14,3)
% V(21,3) V(22,3) V(23,3) V(24,3)
% V(31,3) V(32,3) V(33,3) V(34,3)
% V(41,3) V(42,3) V(43,3) V(44,3)

%Checking minimum value of wind velocity
windvelmin(windvelmin<0)=0;

%Checking if windvel is 1d-row (1 time step)
%If windvel has 1 time step, then it should be 1d-row array
if iscolumn(windvel)==1, 
    windvel=windvel';
end

%Checking if winddir is 1d-row (1 time step)
%If winddir has 1 time step, then it should be 1d-row array
if iscolumn(winddir)==1, 
    winddir=winddir';
end

[M,N]=size(xgrid);
[K,L]=size(windvel);
windvelpointx=zeros(K,L); %Pre-assigning an array
windvelpointy=zeros(K,L); %Pre-assigning an array

%K is number of time steps for a time series
for i=1:K
    
    windveltimestep=(windvel(i,:))';
    winddirtimestep=(winddir(i,:))';
    
    %Checking for the minimum wind velocity
    windveltimestep(windveltimestep<windvelmin)=windvelmin;
    
    %Meteorological direction
    if strcmp(winddirtype,'mete')==1
        %Convert wind data from wind source to wind destination
        %Converting meteorological direction to trigonometric direction
        winddirtimesteptrig=-winddirtimestep+270;

        %Add 360 to all numbers to have them all positive
        %Use mod(360) to take care of the ones larger than 360
        winddirtimesteptrig=mod((winddirtimesteptrig+360),360); 

    %Trigonometric direction
    elseif strcmp(winddirtype,'trig')==1
        winddirtimesteptrig=winddirtimestep;
    
    end

    %Converting wind velocity into x and y components
    winddirtimesteptrigrad=deg2rad(winddirtimesteptrig); %Converting wind direction to Radian
    [windveltimestepx,windveltimestepy]=pol2cart(winddirtimesteptrigrad,windveltimestep); %Calculating x and y components of wind velocity
    
    %x direction
    %Interpolating data into grid using default or linear method
    if strcmp(CalcMethod,'linear')==1
        [xgrid,ygrid,windvelgridx]=griddata(xpoint,ypoint,windveltimestepx,xgrid,ygrid);

        %Replacing NaN data point resulted from default method with ones from nearest method
        if sum(isnan(windvelgridx(:)))>0
            [xgrid,ygrid,windvelgridxnearest]=griddata(xpoint,ypoint,windveltimestepx,xgrid,ygrid,'nearest');
            windvelgridx(isnan(windvelgridx)==1)=windvelgridxnearest(isnan(windvelgridx)==1); 
        end

    %Interpolating data into grid using nearest neighbor method
    elseif strcmp(CalcMethod,'nearest')==1
        [xgrid,ygrid,windvelgridx]=griddata(xpoint,ypoint,windveltimestepx,xgrid,ygrid,'nearest');

    end

    %y direction
    %Interpolating data into grid using default or linear method
    if strcmp(CalcMethod,'linear')==1
        [xgrid,ygrid,windvelgridy]=griddata(xpoint,ypoint,windveltimestepy,xgrid,ygrid);

        %Replacing NaN data point resulted from default method with ones from nearest method
        if sum(isnan(windvelgridy(:)))>0
            [xgrid,ygrid,windvelgridynearest]=griddata(xpoint,ypoint,windveltimestepy,xgrid,ygrid,'nearest');
            windvelgridy(isnan(windvelgridy)==1)=windvelgridynearest(isnan(windvelgridy)==1); 
        end

    %Interpolating data into grid using nearest neighbor method
    elseif strcmp(CalcMethod,'nearest')==1
        [xgrid,ygrid,windvelgridy]=griddata(xpoint,ypoint,windveltimestepy,xgrid,ygrid,'nearest');

    end

    %Creating the wind
    windgrid=[windvelgridx;windvelgridy];
    if i==1
        swanwind=windgrid;
    else
        swanwind=[swanwind;windgrid];
    end
    
    %Reshaping array into [:,1]
    %xgridrshp=xgrid(:); %Reshape xgrid into [:,1]
    %ygridrshp=ygrid(:); %Reshape ygrid into [:,1]
    %windvelgridxrshp=windvelgridx(:); %Reshape windvelgridx into [:,1]
    %windvelgridyrshp=windvelgridy(:); %Reshape windvelgridy into [:,1]
    xgridrshp=reshape(xgrid,[M*N,1]); %Reshape xgrid into [:,1]
    ygridrshp=reshape(ygrid,[M*N,1]); %Reshape ygrid into [:,1]
    windvelgridxrshp=reshape(windvelgridx,[M*N,1]); %Reshape windvelgridx into [:,1]
    windvelgridyrshp=reshape(windvelgridy,[M*N,1]); %Reshape windvelgridy into [:,1]

    %Export interpolated velocity at each station location
    for j=1:L
        distxy=sqrt((xgridrshp-xpoint(j,1)).^2+(ygridrshp-ypoint(j,1)).^2); %Calculating distance from (x,y) to (xpoint,ypoint)
        [distxymin,Indx]=min(distxy); %Finding minimum distance
        windvelpointx(i,j)=windvelgridxrshp(Indx,1); %Nearest interpolated wind velocity in x direction at each point of (xpoint,ypoint)
        windvelpointy(i,j)=windvelgridyrshp(Indx,1); %Nearest interpolated wind velocity in y direction at each point of (xpoint,ypoint)
    end
    
end

%--------------------------------------------------------------------------
%Converting wind data at point of (xpoint,ypoint) from cartesian to meteorological wind
%In meteorological wind, wind direction is a direction wind comes from and is measured counter-clockwise from the North
%0 (degree): from North, 90 (degree): from East, 180 (degree): from South, 270 (degree): from West

[winddirpointtrigrad,windvelpoint]=cart2pol(windvelpointx,windvelpointy); %Interpolated wind velocity at station locations
%winddirpointtrigrad=atan2(windvelpointx,windvelpointy); %Interpolated wind direction at station locations

winddirpointtrig=rad2deg(winddirpointtrigrad); %Converting wind direction to Degree

%Meteorological direction
if strcmp(winddirtype,'mete')==1
    %Converting trigonometric direction to meteorological direction
    winddirpoint=270-winddirpointtrig;

    %Add 360 to all numbers to have them all positive
    %Use mod(360) to take care of the ones larger than 360
    winddirpoint=mod((winddirpoint+360),360); 

%Trigonometric direction
elseif strcmp(winddirtype,'trig')==1
    winddirpoint=winddirpointtrig;

end

%--------------------------------------------------------------------------
%Saving data

if strcmp(savedata,'yes')==1

    %Changing directory to saving directory
    currentFolder=pwd;
    cd(outfilelocation)

    %Saving data
    %save(outfilename,'swanwind','-ASCII')
    dlmwrite(outfilename,swanwind,'delimiter',' ')

    %Changing directory to working directory
    cd(currentFolder)

end

%--------------------------------------------------------------------------
