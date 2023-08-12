function [swanvectorvariable, Vxpoint, Vypoint] = swanvectorvarspvariedsct(xgrid, ygrid, xpoint, ypoint, Vx, Vy, savedata, outfilename, outfilelocation, CalcMethod)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-12-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

swanvectorvarspvariedsct
========================

.. code:: MATLAB

    [swanvectorvariable, Vxpoint, Vypoint] = swanvectorvarspvariedsct(xgrid, ygrid, xpoint, ypoint, Vx, Vy, savedata, outfilename, outfilelocation, CalcMethod)

Description
-----------

Generate SWAN file for spatially varied vector variable from scattered input data

Inputs
------

xgrid
    x (longitude) of output grid points as a [M*N] array
ygrid
    y (latitude) of output grid points as a [M*N] array
xpoint
    | x (longitude) of the location (such as meteorological station) that vector variable is known in that location
    | as a [L] array
    | 1st element is 1st point, 2nd element is 2nd point, 3rd element is 3rd point, ...
ypoint
    | y (latitude) of the location (such as meteorological station) that vector variable is known in that location
    | as a [L] array
    | 1st element is 1st point, 2nd element is 2nd point, 3rd element is 3rd point, ...
Vx
    | Variable in x direction (x component of input variable) as a [K*L] array
    | K is number of time steps for a time series
    |     1st row (i.e. [1,:]) is 1st time step, 2nd row (i.e. [2,:]) is 2nd time step, ...
    |     Example: windvel(1,:) is the wind velocity data for all points at the first time step
    | L is number of points (such as meteorological stations) that wind velocity is known in those locations
    |     L should be L>=3
    |     1st column (i.e. [:,1]) is 1st point, 2nd column (i.e. [:,2]) is 2nd point, 3rd column (i.e. [:,3]) is 3rd point, ...
    |     Example: windvel(:,1) is the wind velocity data for all time steps at the first point
Vy
    | Variable in y direction (y component of input variable) as a [K*L] array
    | K is number of time steps for a time series
    |     1st row (i.e. [1,:]) is 1st time step, 2nd row (i.e. [2,:]) is 2nd time step, ...
    |     Example: windvel(1,:) is the wind velocity data for all points at the first time step
    | L is number of points (such as meteorological stations) that wind velocity is known in those locations
    |     L should be L>=3
    |     1st column (i.e. [:,1]) is 1st point, 2nd column (i.e. [:,2]) is 2nd point, 3rd column (i.e. [:,3]) is 3rd point, ...
    |     Example: windvel(:,1) is the wind velocity data for all time steps at the first point
savedata='no';
    | Define if save data in a file or not
    | 'no': does not save 
    | 'yes': save data as ascii file
outfilename='swanwind.wnd';
    | Name of output file between ' ' mark, example: 'swanwind.wnd'
    | outfilename should have proper name and extension
outfilelocation=pwd;
    Location of output file between ' ' mark, example: 'C:\' in MATLAB, or 'C:/' in Python
CalcMethod='linear';
    | Interpolation method 
    | 'linear': Use default or 'linear' method to interpolate
    | 'nearest': Use nearest neighbor method to interpolate

Outputs
-------

swanvectorvariable
    | Spatially varied vector variable data formated for SWAN
    | Vector variable data at each time step is assigned into the grid points
Vxpoint
    | Nearest interpolated vector variable in x direction at (xpoint,ypoint) as a [K*L] array
    | K is number of time steps for a time series
    | L is number of points (such as meteorological stations) that vector variable in x direction is known in those locations
Vypoint
    | Nearest interpolated vector variable in y direction at (xpoint,ypoint) as a [K*L] array
    | K is number of time steps for a time series
    | L is number of points (such as meteorological stations) that vector variable in y direction is known in those locations

Examples
--------

.. code:: MATLAB

    [xgrid,ygrid]=meshgrid(linspace(-91,-90,100),linspace(28,30,100));
    windvelx=[[10.5,10.55,10.6];[10.64,10.69,10.74];[10.7,10.75,10.8];[10.4,10.45,10.5]]; %Data for 4 time steps
    windvely=[[1.5,1.55,1.6];[1.64,1.69,1.74];[1.7,1.75,1.8];[1.4,1.45,1.5]]; %Data for 4 time steps
    xpoint=[-90.5;-90.3;-90.7]; %Data are known at 3 locations
    ypoint=[29.2;29;28.8]; %Data are known at 3 locations
    savedata='no';
    outfilename='swanwind.wnd';
    outfilelocation=pwd;
    CalcMethod='linear';
    [swanvectorvariable,windvelxpoint,windvelypoint]=swanvectorvarspvariedsct(xgrid,ygrid,xpoint,ypoint,windvelx,windvely,savedata,outfilename,outfilelocation,CalcMethod);


    [xgrid,ygrid]=meshgrid(linspace(-91,-90,100),linspace(28,30,100));
    windvelx=[10.5,10.55,10.6]; %Data for 1 time step
    windvely=[1.5,1.55,1.6]; %Data for 1 time step
    xpoint=[-90.5;-90.3;-90.7]; %Data are known at 3 locations
    ypoint=[29.2;29;28.8]; %Data are known at 3 locations
    savedata='no';
    outfilename='swanwind.wnd';
    outfilelocation=pwd;
    CalcMethod='linear';
    [swanvectorvariable,windvelxpoint,windvelypoint]=swanvectorvarspvariedsct(xgrid,ygrid,xpoint,ypoint,windvelx,windvely,savedata,outfilename,outfilelocation,CalcMethod);

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
        savedata='no'; outfilename='swanwind.wnd'; outfilelocation=pwd; CalcMethod='linear';
    case 7
        outfilename='swanwind.wnd'; outfilelocation=pwd; CalcMethod='linear';
    case 8
        outfilelocation=pwd; CalcMethod='linear';
    case 9
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
%Generating SWAN file

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

%Checking if Vx is 1d-row (1 time step)
%If Vx has 1 time step, then it should be 1d-row array
if iscolumn(Vx)==1, 
    Vx=Vx';
end

%Checking if Vy is 1d-row (1 time step)
%If Vy has 1 time step, then it should be 1d-row array
if iscolumn(Vy)==1, 
    Vy=Vy';
end

[M,N]=size(xgrid);
[K,L]=size(Vx);
Vxpoint=zeros(K,L); %Pre-assigning an array
Vypoint=zeros(K,L); %Pre-assigning an array

%K is number of time steps for a time series
for i=1:K
    
    Vtimestepx=(Vx(i,:))';
    Vtimestepy=(Vy(i,:))';
    
    %x direction
    %Interpolating data into grid using default or linear method
    if strcmp(CalcMethod,'linear')==1
        [xgrid,ygrid,Vgridx]=griddata(xpoint,ypoint,Vtimestepx,xgrid,ygrid);

        %Replacing NaN data point resulted from default method with ones from nearest method
        if sum(isnan(Vgridx(:)))>0
            [xgrid,ygrid,Vgridxnearest]=griddata(xpoint,ypoint,Vtimestepx,xgrid,ygrid,'nearest');
            Vgridx(isnan(Vgridx)==1)=Vgridxnearest(isnan(Vgridx)==1); 
        end

    %Interpolating data into grid using nearest neighbor method
    elseif strcmp(CalcMethod,'nearest')==1
        [xgrid,ygrid,Vgridx]=griddata(xpoint,ypoint,Vtimestepx,xgrid,ygrid,'nearest');

    end

    %y direction
    %Interpolating data into grid using default or linear method
    if strcmp(CalcMethod,'linear')==1
        [xgrid,ygrid,Vgridy]=griddata(xpoint,ypoint,Vtimestepy,xgrid,ygrid);

        %Replacing NaN data point resulted from default method with ones from nearest method
        if sum(isnan(Vgridy(:)))>0
            [xgrid,ygrid,Vgridynearest]=griddata(xpoint,ypoint,Vtimestepy,xgrid,ygrid,'nearest');
            Vgridy(isnan(Vgridy)==1)=Vgridynearest(isnan(Vgridy)==1); 
        end

    %Interpolating data into grid using nearest neighbor method
    elseif strcmp(CalcMethod,'nearest')==1
        [xgrid,ygrid,Vgridy]=griddata(xpoint,ypoint,Vtimestepy,xgrid,ygrid,'nearest');

    end

    %Creating the wind
    Vgrid=[Vgridx;Vgridy];
    if i==1
        swanvectorvariable=Vgrid;
    else
        swanvectorvariable=[swanvectorvariable;Vgrid];
    end
    
    %Reshaping array into [:,1]
    %xgridrshp=xgrid(:); %Reshape xgrid into [:,1]
    %ygridrshp=ygrid(:); %Reshape ygrid into [:,1]
    %Vgridxrshp=Vgridx(:); %Reshape Vgridx into [:,1]
    %Vgridyrshp=Vgridy(:); %Reshape Vgridy into [:,1]
    xgridrshp=reshape(xgrid,[M*N,1]); %Reshape xgrid into [:,1]
    ygridrshp=reshape(ygrid,[M*N,1]); %Reshape ygrid into [:,1]
    Vgridxrshp=reshape(Vgridx,[M*N,1]); %Reshape Vgridx into [:,1]
    Vgridyrshp=reshape(Vgridy,[M*N,1]); %Reshape Vgridy into [:,1]

    %Export interpolated velocity at each station location
    for j=1:L
        distxy=sqrt((xgridrshp-xpoint(j,1)).^2+(ygridrshp-ypoint(j,1)).^2); %Calculating distance from (x,y) to (xpoint,ypoint)
        [distxymin,Indx]=min(distxy); %Finding minimum distance
        Vxpoint(i,j)=Vgridxrshp(Indx,1); %Nearest interpolated wind velocity in x direction at each point of (xpoint,ypoint)
        Vypoint(i,j)=Vgridyrshp(Indx,1); %Nearest interpolated wind velocity in y direction at each point of (xpoint,ypoint)
    end
    
end

%--------------------------------------------------------------------------
%Saving data

if strcmp(savedata,'yes')==1

    %Changing directory to saving directory
    currentFolder=pwd;
    cd(outfilelocation)

    %Saving data
    %save(outfilename,'swanvectorvariable','-ASCII')
    dlmwrite(outfilename,swanvectorvariable,'delimiter',' ')

    %Changing directory to working directory
    cd(currentFolder)

end

%--------------------------------------------------------------------------
