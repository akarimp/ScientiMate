function [z] = interpgrid2xyz(xgrid, ygrid, zgrid, x, y, dispout)
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

interpgrid2xyz
==============

.. code:: MATLAB

    [z] = interpgrid2xyz(xgrid, ygrid, zgrid, x, y, dispout)

Description
-----------

Interpolate 2d gridded data on given scatter point(s) using nearest neighbor method 

Inputs
------

xgrid
    x data as a [M*N] array
ygrid
    y data as a [M*N] array
zgrid
    z data as z(x,y) as a [M*N] array
x
    x of point that nearest point to that is desired to be found
y
    y of point that nearest point to that is desired to be found
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

z
    Value of interpolated data at (x,y) as z(x,y)

Examples
--------

.. code:: MATLAB

    [xgrid,ygrid]=meshgrid(linspace(0,10,100),linspace(0,10,100));
    zgrid=ygrid.*sin(xgrid)-xgrid.*cos(ygrid);
    x=10.*rand(100,1);
    y=10.*rand(100,1);
    [z]=interpgrid2xyz(xgrid,ygrid,zgrid,x,y,'yes');

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
    case 5
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(x)==1
    x=x';
end

if isrow(y)==1
    y=y';
end

%--------------------------------------------------------------------------
%Calculating z(x,y) values

%Interpolating data into grid using nearest neighbor method
%First method
z=interp2(xgrid,ygrid,zgrid,x,y,'nearest');

%Second method
%Mapping values of x and y to associated index
%[M,N]=size(xgrid);
%xmap(:,1)=interp1([nanmin(xgrid(:));nanmax(xgrid(:))],[1;N],x,'extrap');
%ymap(:,1)=interp1([nanmin(ygrid(:));nanmax(ygrid(:))],[1;M],y,'extrap');
%xmap(xmap<1)=1;
%xmap(xmap>N)=N;
%ymap(ymap<1)=1;
%ymap(ymap>M)=M;
%z=diag(zgrid(fix(ymap),fix(xmap)));

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1
    
    %Plotting data
    surf(xgrid,ygrid,zgrid)
    hold on
    scatter3(x,y,z,'r','filled')
    
    colorbar
    xlabel('x')
    ylabel('y')
    zlabel('z')
    legend('Input Data','Interpolated Data')
    
end

%--------------------------------------------------------------------------
