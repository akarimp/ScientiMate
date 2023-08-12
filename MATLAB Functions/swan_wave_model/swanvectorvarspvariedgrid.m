function [swanvectorvariable] = swanvectorvarspvariedgrid(xgrid, ygrid, xpointgrid, ypointgrid, Vxgrid, Vygrid, savedata, outfilename, outfilelocation, CalcMethod)
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

swanvectorvarspvariedgrid
=========================

.. code:: MATLAB

    [swanvectorvariable] = swanvectorvarspvariedgrid(xgrid, ygrid, xpointgrid, ypointgrid, Vxgrid, Vygrid, savedata, outfilename, outfilelocation, CalcMethod)

Description
-----------

Generate SWAN file for spatially varied vector variable from gridded input data

Inputs
------

xgrid
    x (longitude) of output grid points as a [M*N] array
ygrid
    y (latitude) of output grid points as a [M*N] array
xpointgrid
    x (longitude) of the locations that vector variable is known in those locations as a [K*L] array
ypointgrid
    y (latitude) of the locations that vector variable is known in those locations as a [K*L] array
Vxgrid
    | Variable in x direction (x component of input variable) at (xpointgrid,ypointgrid) as a [K*L*P] array
    | P is number of time steps for a time series
Vygrid
    | Variable in y direction (y component of input variable) at (xpointgrid,ypointgrid) as a [K*L*P] array in (Degree)
    | P is number of time steps for a time series
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

Examples
--------

.. code:: MATLAB

    [xgrid,ygrid]=meshgrid(linspace(-91,-90,100),linspace(28,30,100));
    [xpointgrid,ypointgrid]=meshgrid(linspace(-92,-89,100),linspace(27,31,100));
    windvelxgrid=10+(12-10).*rand(100,100,4); %Data for 4 time steps
    windvelygrid=1+(2-1).*rand(100,100,4); %Data for 4 time steps
    savedata='no';
    outfilename='swanwind.wnd';
    outfilelocation=pwd;
    CalcMethod='linear';
    [swanvectorvariable]=swanvectorvarspvariedgrid(xgrid,ygrid,xpointgrid,ypointgrid,windvelxgrid,windvelygrid,savedata,outfilename,outfilelocation,CalcMethod);


    [xgrid,ygrid]=meshgrid(linspace(-91,-90,100),linspace(28,30,100));
    [xpointgrid,ypointgrid]=meshgrid(linspace(-92,-89,100),linspace(27,31,100));
    windvelxgrid=10+(12-10).*rand(100,100); %Data for 1 time step
    windvelygrid=1+(2-1).*rand(100,100); %Data for 1 time step
    savedata='no';
    outfilename='swanwind.wnd';
    outfilelocation=pwd;
    CalcMethod='linear';
    [swanvectorvariable]=swanvectorvarspvariedgrid(xgrid,ygrid,xpointgrid,ypointgrid,windvelxgrid,windvelygrid,savedata,outfilename,outfilelocation,CalcMethod);

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

%if isrow(xpoint)==1
%    xpoint=xpoint';
%end

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

[M,N]=size(xgrid);
[K,L,P]=size(Vxgrid);

%Reshaping array into [K*L,1]
xpoint=reshape(xpointgrid,[K*L,1]);
ypoint=reshape(ypointgrid,[K*L,1]);

%P is number of time steps for a time series
for i=1:P
    
    %Reshaping array into [K*L,1]
    Vxtimestep=reshape(Vxgrid(:,:,i),[K*L,1]);
    Vytimestep=reshape(Vygrid(:,:,i),[K*L,1]);
    
    %x direction
    %Interpolating data into grid using default or linear method
    if strcmp(CalcMethod,'linear')==1
        [xgrid,ygrid,Vgridx]=griddata(xpoint,ypoint,Vxtimestep,xgrid,ygrid);

        %Replacing NaN data point resulted from default method with ones from nearest method
        if sum(isnan(Vgridx(:)))>0
            [xgrid,ygrid,Vgridxnearest]=griddata(xpoint,ypoint,Vxtimestep,xgrid,ygrid,'nearest');
            Vgridx(isnan(Vgridx)==1)=Vgridxnearest(isnan(Vgridx)==1); 
        end

    %Interpolating data into grid using nearest neighbor method
    elseif strcmp(CalcMethod,'nearest')==1
        [xgrid,ygrid,Vgridx]=griddata(xpoint,ypoint,Vxtimestep,xgrid,ygrid,'nearest');

    end

    %y direction
    %Interpolating data into grid using default or linear method
    if strcmp(CalcMethod,'linear')==1
        [xgrid,ygrid,Vgridy]=griddata(xpoint,ypoint,Vytimestep,xgrid,ygrid);

        %Replacing NaN data point resulted from default method with ones from nearest method
        if sum(isnan(Vgridy(:)))>0
            [xgrid,ygrid,Vgridynearest]=griddata(xpoint,ypoint,Vytimestep,xgrid,ygrid,'nearest');
            Vgridy(isnan(Vgridy)==1)=Vgridynearest(isnan(Vgridy)==1); 
        end

    %Interpolating data into grid using nearest neighbor method
    elseif strcmp(CalcMethod,'nearest')==1
        [xgrid,ygrid,Vgridy]=griddata(xpoint,ypoint,Vytimestep,xgrid,ygrid,'nearest');

    end

    %Creating the output
    Vgrid=[Vgridx;Vgridy];
    if i==1
        swanvectorvariable=Vgrid;
    else
        swanvectorvariable=[swanvectorvariable;Vgrid];
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
