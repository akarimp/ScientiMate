function [xpoints, ypoints, distxy] = pointscart(x1, y1, x2, y2, nmidpoints, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-08-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

pointscart
==========

.. code:: MATLAB

    [xpoints, ypoints, distxy] = pointscart(x1, y1, x2, y2, nmidpoints, dispout)

Description
-----------

Generate points between point (x1,y1) and (x2,y2) on cartesian coordinate

Inputs
------

x1
    x of start point (first point)
y1
    y of start point (first point)
x2
    x of end point (last point) 
y2
    y of end point (last point) 
nmidpoints=1;
    | Number of middle points generated between first point and last point
    | nmidpoints=1;: 1 point between first point and last point is generated
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

xpoint
    | x of points from (x1,y1) and (x2,y2)
    | x1 and x2 are included, total number of points are nmidpoints+2
ypoint
    | y of points from (x1,y1) and (x2,y2)
    | y1 and y2 are included, total number of points are nmidpoints+2
distxy
    Distance at each points of (xpoint,ypoint) from the first point of (x1,y1)

Examples
--------

.. code:: MATLAB

    x1=-3;
    y1=-3;
    x2=3;
    y2=3;
    [xpoints,ypoints,distxy]=pointscart(x1,y1,x2,y2,10,'yes');

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
        nmidpoints=1; dispout='no';
    case 5
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

%if isrow(x)==1
%    x=x';
%end

%--------------------------------------------------------------------------
%Generating points from (x1,y1) to (x2,y2)

nmidpoints(nmidpoints<1)=1;

xpoints(:,1)=linspace(x1,x2,nmidpoints+2);
ypoints(:,1)=linspace(y1,y2,nmidpoints+2);

%--------------------------------------------------------------------------
%Calculating distance from (x1,y1) to (x2,y2)

%Distance at each points of (xpoint,ypoint) from the first point of (x1,y1)
distxy(:,1)=sqrt((xpoints-xpoints(1,1)).^2+(ypoints-ypoints(1,1)).^2); 

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1
    
    %Plotting data
    scatter(x1,y1,'filled')
    hold on
    scatter(x2,y2,'filled')
    scatter(xpoints(2:end-1,1),ypoints(2:end-1,1),'filled')
    
    xlabel('x')
    ylabel('y')
    legend('Start Point','End Point')
    
end

%--------------------------------------------------------------------------
