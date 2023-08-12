function [swanwaterlevel] = swanwaterlevelspvariedgrid(xgrid, ygrid, xpointgrid, ypointgrid, waterlevelgrid, savedata, outfilename, outfilelocation, CalcMethod)
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

swanwaterlevelspvariedgrid
==========================

.. code:: MATLAB

    [swanwaterlevel] = swanwaterlevelspvariedgrid(xgrid, ygrid, xpointgrid, ypointgrid, waterlevelgrid, savedata, outfilename, outfilelocation, CalcMethod)

Description
-----------

| Generate SWAN water level file for spatially varied water level from gridded input data
| This function can be used for any other scalar variable as well

Inputs
------

xgrid
    x (longitude) of output grid points as a [M*N] array
ygrid
    y (latitude) of output grid points as a [M*N] array
xpointgrid
    x (longitude) of the locations that water level (or any scalar variable) is known in those locations as a [K*L] array
ypointgrid
    y (latitude) of the locations that water level (or any scalar variable) is known in those locations as a [K*L] array
waterlevelgrid
    | Water level (or any scalar variable) at (xpointgrid,ypointgrid) as a [K*L*P] array
    | P is number of time steps for a time series
savedata='no';
    | Define if save data in a file or not
    | 'no': does not save 
    | 'yes': save data as ascii file
outfilename='swanwaterlevel.wl';
    | Name of output file between ' ' mark, example: 'swanwaterlevel.wl'
    | outfilename should have '.wl' extension
    | In case of using other scalar varable than water level, use proper name and extension
outfilelocation=pwd;
    Location of output file between ' ' mark, example: 'C:\' in MATLAB, or 'C:/' in Python
CalcMethod='linear';
    | Interpolation method 
    | 'linear': Use default or 'linear' method to interpolate
    | 'nearest': Use nearest neighbor method to interpolate

Outputs
-------

swanwaterlevel
    | Spatially varied water level data (or any scalar variable) formated for SWAN
    | Water level data at each time step is assigned into the grid points

Examples
--------

.. code:: MATLAB

    [xgrid,ygrid]=meshgrid(linspace(-91,-90,100),linspace(28,30,100));
    [xpointgrid,ypointgrid]=meshgrid(linspace(-92,-89,100),linspace(27,31,100));
    waterlevelgrid=0.5+(0.6-0.5).*rand(100,100,4); %Data for 4 time steps
    savedata='no';
    outfilename='swanwaterlevel.wl';
    outfilelocation=pwd;
    CalcMethod='linear';
    [swanwaterlevel]=swanwaterlevelspvariedgrid(xgrid,ygrid,xpointgrid,ypointgrid,waterlevelgrid,savedata,outfilename,outfilelocation,CalcMethod);


    [xgrid,ygrid]=meshgrid(linspace(-91,-90,100),linspace(28,30,100));
    [xpointgrid,ypointgrid]=meshgrid(linspace(-92,-89,100),linspace(27,31,100));
    waterlevelgrid=0.5+(0.6-0.5).*rand(100,100); %Data for 1 time step
    savedata='no';
    outfilename='swanwaterlevel.wl';
    outfilelocation=pwd;
    CalcMethod='linear';
    [swanwaterlevel]=swanwaterlevelspvariedgrid(xgrid,ygrid,xpointgrid,ypointgrid,waterlevelgrid,savedata,outfilename,outfilelocation,CalcMethod);

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
    case 5
        savedata='no'; outfilename='swanwaterlevel.wl'; outfilelocation=pwd; CalcMethod='linear';
    case 6
        outfilename='swanwaterlevel.wl'; outfilelocation=pwd; CalcMethod='linear';
    case 7
        outfilelocation=pwd; CalcMethod='linear';
    case 8
        CalcMethod='linear';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

%if isrow(xpoint)==1
%    xpoint=xpoint';
%end

%--------------------------------------------------------------------------
%Generating SWAN water level file (.wl)

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
%Then SWAN water level file (.wl) for waterlevel with three time steps has a format of (RowColumn,TimeStep):
% wl(11,1) wl(12,1) wl(13,1) wl(14,1)
% wl(21,1) wl(22,1) wl(23,1) wl(24,1)
% wl(31,1) wl(32,1) wl(33,1) wl(34,1)
% wl(41,1) wl(42,1) wl(43,1) wl(44,1)
% wl(11,2) wl(12,2) wl(13,2) wl(14,2)
% wl(21,2) wl(22,2) wl(23,2) wl(24,2)
% wl(31,2) wl(32,2) wl(33,2) wl(34,2)
% wl(41,2) wl(42,2) wl(43,2) wl(44,2)
% wl(11,3) wl(12,3) wl(13,3) wl(14,3)
% wl(21,3) wl(22,3) wl(23,3) wl(24,3)
% wl(31,3) wl(32,3) wl(33,3) wl(34,3)
% wl(41,3) wl(42,3) wl(43,3) wl(44,3)

[M,N]=size(xgrid);
[K,L,P]=size(waterlevelgrid);

%Reshaping array into [K*L,1]
xpoint=reshape(xpointgrid,[K*L,1]);
ypoint=reshape(ypointgrid,[K*L,1]);

%P is number of time steps for a time series
for i=1:P
    
    %Reshaping array into [K*L,1]
    waterleveltimestep=reshape(waterlevelgrid(:,:,i),[K*L,1]);
    
    %Interpolating data into grid using default or linear method
    if strcmp(CalcMethod,'linear')==1
        [xgrid,ygrid,zgrid]=griddata(xpoint,ypoint,waterleveltimestep,xgrid,ygrid);

        %Replacing NaN data point resulted from default method with ones from nearest method
        if sum(isnan(zgrid(:)))>0
            [xgrid,ygrid,zgridnearest]=griddata(xpoint,ypoint,waterleveltimestep,xgrid,ygrid,'nearest');
            zgrid(isnan(zgrid)==1)=zgridnearest(isnan(zgrid)==1); 
        end

    %Interpolating data into grid using nearest neighbor method
    elseif strcmp(CalcMethod,'nearest')==1
        [xgrid,ygrid,zgrid]=griddata(xpoint,ypoint,waterleveltimestep,xgrid,ygrid,'nearest');

    end
    
    %Creating the water level
    if i==1
        swanwaterlevel=zgrid;
    else
        swanwaterlevel=[swanwaterlevel;zgrid];
    end
    
end

%--------------------------------------------------------------------------
%Saving data

if strcmp(savedata,'yes')==1

    %Changing directory to saving directory
    currentFolder=pwd;
    cd(outfilelocation)

    %Saving data
    %save(outfilename,'swanwaterlevel','-ASCII')
    dlmwrite(outfilename,swanwaterlevel,'delimiter',' ')

    %Changing directory to working directory
    cd(currentFolder)

end

%--------------------------------------------------------------------------
