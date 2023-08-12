function [swandepth, swangrid, ncellx, ncelly, ngridx, ngridy] = swandepthgrid(xgrid, ygrid, zgrid, zmin, zmax, nanreplacement, savedata, xyoutfilename, zoutfilename, outfilelocation, dispout)
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

swandepthgrid
=============

.. code:: MATLAB

    [swandepth, swangrid, ncellx, ncelly, ngridx, ngridy] = swandepthgrid(xgrid, ygrid, zgrid, zmin, zmax, nanreplacement, savedata, xyoutfilename, zoutfilename, outfilelocation, dispout)

Description
-----------

Generate SWAN depth file and its associated x-y grid file

Inputs
------

xgrid
    x (longitude) of grid points as a [M*N] array
ygrid
    y (latitude) of grid points as a [M*N] array
zgrid
    | z (elevation) of grid points as a [M*N] array
    | z>0 is land, z<0 is water, z=0 is water surface
zmin=nanmin(zgrid);
    | Minimum z (elevation) to be considered
    | All z<zmin would be set to NaN
zmax=nanmax(zgrid);
    | Maximum z (elevation) to be considered
    | All z>zmax would be set to NaN
    | Example: to remove all land values from z data, set zmax=0
nanreplacement='no';
    | Replace NaN values with nanreplacement
    | By seting up an exception value equal to nanreplacement in SWAN input file, SWAN disregards these data points
    | If nanreplacement='no', then it does not replace NaN values
    | Example, nanreplacement=-999 will replace all NaN values with -999
    |     Then if -999 is set as an exception in SWAN input file, SWAN disregards all -999 values
savedata='no';
    | Define if save data in a file or not
    | 'no': does not save 
    | 'yes': save data as ascii file
xyoutfilename='swangrid.xy';
    | Name of output grid file between ' ' mark, example: 'swangrid.xy'
    | xyoutfilename should have '.xy' extension
zoutfilename='swandepth.dep';
    | Name of output depth file between ' ' mark, example: 'swandepth.dep'
    | zoutfilename should have '.dep' extension
outfilelocation=pwd;
    Location of output file between ' ' mark, example: 'C:\' in MATLAB, or 'C:/' in Python
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

swandepth
    | Depth data formated for SWAN
    | Note: SWAN (Delft3D) needs depth data not bathymetry data. 
    |     It means value below water surface should be positive and above surface negative.
    | Note: All NaN values are replaced with nanreplacement
    |     Set up an exception value equal to nanreplacement in SWAN to disregard these data points
swangrid
    Grid data formated for SWAN
ncellx
    | Number of cells in x direction (ncellx=ngridx-1)
    | In SWAN, ncellx is equal to a number of meshes in computational grid in x direction 
ncelly
    | Number of cells in y direction (ncelly=ngridy-1)
    | In SWAN, ncelly is equal to a number of meshes in computational grid in y direction 
ngridx
    Number of grid points in x direction
ngridy
    Number of grid points in y direction

Examples
--------

.. code:: MATLAB

    x(:,1)=10.*rand(1000,1);
    y(:,1)=10.*rand(1000,1);
    z=x.^2+y.^2-mean(x.^2+y.^2);
    [xgrid,ygrid]=meshgrid(linspace(min(x),max(x),100),linspace(min(y),max(y),100));
    zgrid=griddata(x,y,z,xgrid,ygrid);
    zmin=nanmin(zgrid(:));
    zmax=nanmax(zgrid(:));
    nanreplacement=-999;
    savedata='no';
    xyoutfilename='swangrid.xy';
    zoutfilename='swandepth.dep';
    outfilelocation=pwd;
    [swandepth,swangrid,ncellx,ncelly,ngridx,ngridy]=swandepthgrid(xgrid,ygrid,zgrid,zmin,zmax,nanreplacement,savedata,xyoutfilename,zoutfilename,outfilelocation,'yes');

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
    case 3
        zmin=nanmin(zgrid(:)); zmax=nanmax(zgrid(:)); nanreplacement='no'; savedata='no'; xyoutfilename='swangrid.xy'; zoutfilename='swandepth.dep'; outfilelocation=pwd; dispout='no';
    case 4
        zmax=nanmax(zgrid(:)); nanreplacement='no'; savedata='no'; xyoutfilename='swangrid.xy'; zoutfilename='swandepth.dep'; outfilelocation=pwd; dispout='no';
    case 5
        nanreplacement='no'; savedata='no'; xyoutfilename='swangrid.xy'; zoutfilename='swandepth.dep'; outfilelocation=pwd; dispout='no';
    case 6
        savedata='no'; xyoutfilename='swangrid.xy'; zoutfilename='swandepth.dep'; outfilelocation=pwd; dispout='no';
    case 7
        xyoutfilename='swangrid.xy'; zoutfilename='swandepth.dep'; outfilelocation=pwd; dispout='no';
    case 8
        zoutfilename='swandepth.dep'; outfilelocation=pwd; dispout='no';
    case 9
        outfilelocation=pwd; dispout='no';
    case 10
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

%if isrow(x)==1 
%    x=x';
%end

%--------------------------------------------------------------------------
%Checking input data

%Makeing sure all z data are in range of [zmin,zmax]
zgrid(zgrid<zmin)=NaN;
zgrid(zgrid>zmax)=NaN;

%--------------------------------------------------------------------------
%Generating SWAN grid file (.xy)

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
%Then SWAN x-y grid file (.xy) for this xgrid and ygrid has a format of:
% x1 x2 x3 x4
% x1 x2 x3 x4
% x1 x2 x3 x4
% x1 x2 x3 x4
% y1 y1 y1 y1
% y2 y2 y2 y2
% y3 y3 y3 y3
% y4 y4 y4 y4

%Generating grid file for SWAN, X=xgrid, Y=ygrid
swangrid=[xgrid;ygrid];

%--------------------------------------------------------------------------
%Generating SWAN depth file (.dep)

%SWAN (Delft3D) needs depth data not bathymetry data. 
%It means value below water surface should be positive and above surface negative.
swandepth=-zgrid;

%--------------------------------------------------------------------------
%Creating exception points by replacing NaN with nanreplacement 
%SWAN disregards data point with value equal to nanreplacement

if strcmp(nanreplacement,'no')~=1
    xgrid(isnan(xgrid)==1)=nanreplacement;
    ygrid(isnan(ygrid)==1)=nanreplacement;
    swandepth(isnan(swandepth)==1)=nanreplacement;
end

%--------------------------------------------------------------------------
%Number of grid points and cells 

%First method
%ngridx: number of grid points in x direction
%ngridy: number of grid points in y direction
[ngridy,ngridx]=size(xgrid);

ncellx=ngridx-1; %Number of cells in x direction
ncelly=ngridy-1; %Number of cells in y direction


%Second method
%ngridx=length(xgrid(1,:)); %Number of grid points in x direction (Number of cell in x direction is ngridx-1)
%ngridy=length(ygrid(:,1)); %Number of grid points in y direction (Number of cell in y direction is ngridy-1)

%ncellx=length(xgrid(1,:))-1; %Number of cells in x direction (Number of cell in x direction is ngridx-1)
%ncelly=length(ygrid(:,1))-1; %Number of cells in y direction (Number of cell in y direction is ngridy-1)

%--------------------------------------------------------------------------
%Saving data

if strcmp(savedata,'yes')==1

    %Changing directory to saving directory
    currentFolder=pwd;
    cd(outfilelocation)

    %Saving data
    %save(xyoutfilename,'swangrid','-ASCII')
    dlmwrite(xyoutfilename,swangrid,'delimiter',' ')
    %save(zoutfilename,'swandepth','-ASCII')
    dlmwrite(zoutfilename,swandepth,'delimiter',' ')

    %Changing directory to working directory
    cd(currentFolder)

end

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1
    
    subplot(1,2,1)
    %[Cont,hCont]=contourf(xgrid,ygrid,zgrid);
    %set(hCont,'LineColor','none');

    %pcol=pcolor(xgrid,ygrid,zgrid);
    %set(pcol,'edgecolor','none')

    imagesc(xgrid,ygrid,zgrid);
    set(gca,'YDir','normal'); %Correct y axis direction

    xlabel('x')
    ylabel('y')
    title('Bathymetry')
    colorbar
    

    subplot(1,2,2)
    %[Cont,hCont]=contourf(xgrid,ygrid,swandepth);
    %set(hCont,'LineColor','none');

    %pcol=pcolor(xgrid,ygrid,swandepth);
    %set(pcol,'edgecolor','none')

    imagesc(xgrid,ygrid,swandepth);
    set(gca,'YDir','normal'); %Correct y axis direction

    xlabel('x')
    ylabel('y')
    title('Depth')
    colorbar

end

%--------------------------------------------------------------------------
