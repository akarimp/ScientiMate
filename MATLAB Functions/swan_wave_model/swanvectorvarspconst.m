function [swanvectorvariable] = swanvectorvarspconst(Vx, Vy, savedata, outfilename, outfilelocation)
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

swanvectorvarspconst
====================

.. code:: MATLAB

    [swanvectorvariable] = swanvectorvarspconst(Vx, Vy, savedata, outfilename, outfilelocation)

Description
-----------

Generate SWAN file for spatially constant vector variable

Inputs
------

Vx
    | Variable in x direction (x component of input variable)
    | If size of Vx>1, then it is considered as a time series
    | 1st element is 1st time step, 2nd element is 2nd time step, ...
Vy
    | Variable in y direction (y component of input variable)
    | Should have a same size as Vx
savedata='no';
    | Define if save data in a file or not
    | 'no': does not save 
    | 'yes': save data as ascii file
outfilename='swanwind.wnd';
    | Name of output file between ' ' mark, example: 'swanwind.wnd'
    | outfilename should have proper name and extension
outfilelocation=pwd;
    Location of output file between ' ' mark, example: 'C:\' in MATLAB, or 'C:/' in Python

Outputs
-------

swanvectorvariable
    | Spatially constant vector variable formated for SWAN
    | Note: Vector variable at each time step is assigned into 4 points,
    |     assuming the vector variable domain is defined by 4 points, one at each corner

Examples
--------

.. code:: MATLAB

    windvelx=[10.5;10.6;10.55]; %Data for 3 time steps
    windvely=[2.5;2.6;2.55]; %Data for 3 time steps
    savedata='no';
    outfilename='swanwind.wnd';
    outfilelocation=pwd;
    [swanvectorvariable]=swanvectorvarspconst(windvelx,windvely,savedata,outfilename,outfilelocation);

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
    case 2
        savedata='no'; outfilename='swanwind.wnd'; outfilelocation=pwd;
    case 3
        outfilename='swanwind.wnd'; outfilelocation=pwd;
    case 4
        outfilelocation=pwd;
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(Vx)==1
    Vx=Vx';
end

if isrow(Vy)==1
    Vy=Vy';
end

%--------------------------------------------------------------------------
%Generating SWAN file

%If Vx=[U1;U2;U3]; and Vx=[V1;V2;V3];
%Then SWAN wind file (.wnd) for this windvel with three time steps has a format of:
% U1 U1
% U1 U1
% V1 V1
% V1 V1
% U2 U2
% U2 U2
% V2 V2
% V2 V2
% U3 U3
% U3 U3
% V3 V3
% V3 V3

%Creating the wind
M=length(Vx(:,1));
swanvectorvariable=zeros(4*M,2);
for i=1:length(Vx(:,1))

    %Wind velocity at time step i, first method
    swanvectorvariable(4*i-3:4*i-2,1:2)=Vx(i,1); %x wind velocity at time step i
    swanvectorvariable(4*i-1:4*i,1:2)=Vy(i,1); %y wind velocity at time step i

    %Wind velocity at time step i, second method
    %x wind velocity at time step i
    %swanvectorvariable(4*i-3,1)=Vx(i,1);
    %swanvectorvariable(4*i-3,2)=Vx(i,1);
    %swanvectorvariable(4*i-2,1)=Vx(i,1);
    %swanvectorvariable(4*i-2,2)=Vx(i,1);

    %y wind velocity at time step i
    %swanvectorvariable(4*i-1,1)=Vy(i,1);
    %swanvectorvariable(4*i-1,2)=Vy(i,1);
    %swanvectorvariable(4*i-0,1)=Vy(i,1);
    %swanvectorvariable(4*i-0,2)=Vy(i,1);

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