function [xgrid, ygrid] = gridgenerator(xmin, xmax, ymin, ymax, gridsize, gridsizetype, dispout)
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

gridgenerator
=============

.. code:: MATLAB

    [xgrid, ygrid] = gridgenerator(xmin, xmax, ymin, ymax, gridsize, gridsizetype, dispout)

Description
-----------

Generate 2d x-y grid

Inputs
------

xmin
    Minimum x of the domain to be generated
xmax
    Maximum x of the domain to be generated
ymin
    Minimum y of the domain to be generated
ymax
    Maximum y of the domain to be generated
gridsize=100;
    | Grid size in x (longitude) and y (latitude) directions to interpolate elevation data on them
    |     if gridsizetype='length' then gridsize is a distance between grid points
    |     if gridsizetype='points' then gridsize is number of grid points in each direction
gridsizetype='points';
    | Grid size type 
    |     'number': gridsize is considered as number of grid points in each direction
    |     'length': gridsize is considered as length between grid points
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

xgrid
    x of the defined mesh
ygrid
    y of the defined mesh

Examples
--------

.. code:: MATLAB

    [xgrid,ygrid]=gridgenerator(0,100,0,100,100,'points','yes');
    [xgrid,ygrid]=gridgenerator(0,100,0,100,20,'length','yes');

References
----------


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
    case 4
        gridsize=100; gridsizetype='points'; dispout='no';
    case 5
        gridsizetype='points'; dispout='no';
    case 6
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors


%--------------------------------------------------------------------------
%Generating grid

%Calculating grid size in each direction using grid length
if strcmp(gridsizetype,'length')==1

    %Checking x and y range to generate correct number of grid points
    if mod((xmax-xmin),gridsize)~=0
        xmax=xmax-mod((xmax-xmin),gridsize);
    end

    if mod((ymax-ymin),gridsize)~=0
        ymax=ymax-mod((ymax-ymin),gridsize);
    end

    gridsizex=gridsize; %Grid size in lat (x) direction
    gridsizey=gridsize; %Grid size in lon (y) direction

    ngridx=((xmax-xmin)/gridsizex)+1; %Number of grid points in each direction
    ngridy=((ymax-ymin)/gridsizey)+1; %Number of grid points in each direction

    %Generating grid
    [xgrid,ygrid]=meshgrid([xmin:gridsizex:xmax],[ymin:gridsizey:ymax]);

%Calculating grid size in each direction using number grid points 
elseif strcmp(gridsizetype,'points')==1

    ngridx=gridsize; %Number of grid points in each direction
    ngridy=gridsize; %Number of grid points in each direction

    gridsizex=(xmax-xmin)/(ngridx-1); %Grid size in x direction
    gridsizey=(ymax-ymin)/(ngridy-1); %Grid size in y direction

    %Generating grid
    [xgrid,ygrid]=meshgrid(linspace(xmin,xmax,ngridx),linspace(ymin,ymax,ngridy));

end

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1
    
    [M,N]=size(xgrid);

    %Plotting data
    for i=1:N
        plot([xgrid(1,i),xgrid(1,i)],[min(ygrid(:)),max(ygrid(:))])
        hold on
    end

    for i=1:M
        plot([min(xgrid(:)),max(xgrid(:))],[ygrid(i,1),ygrid(i,1)])
        hold on
    end
        
    xlabel('x')
    ylabel('y')
    
end

%--------------------------------------------------------------------------
