function [distxy, theta] = distancecart(x1, y1, x2, y2, CalcMethod, dispout)
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

distancecart
============

.. code:: MATLAB

    [distxy, theta] = distancecart(x1, y1, x2, y2, CalcMethod, dispout)

Description
-----------

Calculate distance from (x1,y1) to (x2,y2) on cartesian coordinate

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
CalcMethod='1d';
    | Distance calculation method 
    | '1d': use 1d array
    | 'pdist2': Use 2d distance function
    | 'vector': Use vectorized distance 
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

distxy
    | Distance from (x1,y1) to (x2,y2)
    | returns M*N array where M=length(x1) and N=length(x2)
    | mth row associated with mth point in (x,y)
    | nth column is associated with nth point in (x2,y2)
theta
    | Angle from start point to end point in (Degree)
    | returns M*N array where M=length(x1) and N=length(x2)
    | mth row associated with mth point in (x,y)
    | nth column is associated with nth point in (x2,y2)

Examples
--------

.. code:: MATLAB

    x1(:,1)=10.*rand(100,1);
    y1(:,1)=10.*rand(100,1);
    x2=[2.5;5;7.5];
    y2=[3;6;9];
    [distxy,theta]=distancecart(x1,y1,x2,y2,'1d','yes');

    x1(:,1)=10.*rand(100,1);
    y1(:,1)=10.*rand(100,1);
    x2(:,1)=100.*rand(10,1);
    y2(:,1)=100.*rand(10,1);
    [distxy,theta]=distancecart(x1,y1,x2,y2,'pdist2','yes');

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
        CalcMethod='1d'; dispout='no';
    case 5
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

if isrow(x2)==1
    x2=x2';
end

if isrow(y2)==1
    y2=y2';
end

%--------------------------------------------------------------------------
%Calculating distance from (x1,y1) to (x2,y2)

M=length(x1(:,1));
N=length(x2(:,1));

%First method, using 1d array
if strcmp(CalcMethod,'1d')==1

    distxy=zeros(M,N); %Pre-assigning array
    for i=1:N
        distxy(:,i)=sqrt((x1-x2(i,1)).^2+(y1-y2(i,1)).^2); %Calculating distance from (x1,y1) to (x2,y2)
    end

%Second method, using pdist2 function
elseif strcmp(CalcMethod,'pdist2')==1

    X=[x1,y1];
    Y=[x2,y2];
    distxy=pdist2(X,Y); %Calculating distance from (x1,y1) to (x2,y2)

%Third method, using vectorized distance 
elseif strcmp(CalcMethod,'vector')==1

    X=[x1,y1];
    Y=[x2,y2];
    
    X2=sum(X.^2,2);
    Y2=sum(Y.^2,2);
    distxy=sqrt(ones(M,1)*Y2'+X2*ones(1,N)-2*X*Y'); %Calculating distance from (x1,y1) to (x2,y2)
    
end

%--------------------------------------------------------------------------
%Calculating angle of the line between start and end points

thetarad=zeros(M,N); %Pre-assigning an array
for i=1:N
    thetarad(:,i)=atan2(y2(i,1)-y1,x2(i,1)-x1); %Angle in radian
end

theta=rad2deg(thetarad); %Angle in degree

%Add 360 to all numbers to have them all positive
%Use mod(360) to take care of the ones larger than 360
theta=mod((theta+360),360); 

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1
    
    %Plotting data
    scatter(x1,y1,'filled')
    hold on
    scatter(x2,y2,'filled')
    
    xlabel('x')
    ylabel('y')
    legend('Start Point','End Point')
    
end

%--------------------------------------------------------------------------
