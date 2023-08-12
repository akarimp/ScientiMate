function [x2, y2, lineslope] = endpointcart(x1, y1, linelength, lineangle, dispout)
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

endpointcart
============

.. code:: MATLAB

    [x2, y2, lineslope] = endpointcart(x1, y1, linelength, lineangle, dispout)

Description
-----------

Find an end point of the straight line segment from its starting point (x1,y1) and its angle on cartesian coordinate

Inputs
------

x1
    x of start point (first point)
y1
    y of start point (first point)
linelength=1;
    Length of a line segment
lineangle=0;
    Angle of a line segment from start point (first point) toward end point (last point) in (Degree)
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

x2
    x of end point (last point) 
y2
    y of end point (last point) 
lineslope
    | Slope of a line segment from start point (first point) toward end point (last point) in (Degree)
    | lineslope=tan(deg2rad(lineangle))

Examples
--------

.. code:: MATLAB

    x1=3;
    y1=3;
    linelength=5;
    lineangle=45;
    [x2,y2,lineslope]=endpointcart(x1,y1,linelength,lineangle,'yes');

    x1=[0;0;0;0;0;0;0;0];
    y1=[0;0;0;0;0;0;0;0];
    linelength=5;
    lineangle=[0;45;90;135;180;225;270;315];
    [x2,y2,lineslope]=endpointcart(x1,y1,linelength,lineangle,'no');

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
    case 2
        linelength=1; lineangle=0; dispout='no';
    case 3
        lineangle=0; dispout='no';
    case 4
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(x1)==1
    x1=x1';
end

if isrow(y1)==1
    y1=y1';
end

if isrow(linelength)==1
    linelength=linelength';
end

if isrow(lineangle)==1
    lineangle=lineangle';
end

%--------------------------------------------------------------------------
%Calculating end point (x2,y2) from start point (x1,y1)

%x at an end of a line segment starts from (x1,y1) with length of linelength and angle of lineangle
x2=x1+linelength.*cos(deg2rad(lineangle)); 

%y at an end of a line segment starts from (x1,y1) with length of linelength and angle of lineangle
y2=y1+linelength.*sin(deg2rad(lineangle));

% Slope of a line segment from (x1,y1) toward (x2,y2) in (Degree)
lineslope=tan(deg2rad(lineangle));

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1
    
    %Plotting data
    scatter(x1,y1,'filled')
    hold on
    scatter(x2,y2,'filled')
    for i=1:length(x1(:,1))
        plot([x1(i,1),x2(i,1)],[y1(i,1),y2(i,1)])
    end
    
    xlabel('x')
    ylabel('y')
    legend('Start Point','End Point')
    
end

%--------------------------------------------------------------------------
