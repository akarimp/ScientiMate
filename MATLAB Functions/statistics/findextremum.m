function [xmin, ymin, xmax, ymax] = findextremum(x, y, winlen, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2018-05-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

findextremum
============

.. code:: MATLAB

    [xmin, ymin, xmax, ymax] = findextremum(x, y, winlen, dispout)

Description
-----------

Find local extremum (minimum and maximum) in data

Inputs
------

x
    x data
y
    y data
winlen=3;
    | Window length, defines a number of points in sliding window used for defining maximum and minimum
    | Example: winlen=5 means two points on each side of each data point is used in calculation
    | Using a larger value for winlen makes it less sensitive 
    | winlen should be an odd number equal or larger than 3
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

xmin
    x of minmum points
ymin
    y of minmum points
xmax
    x of maximum points
ymax
    y of maximum points

Examples
--------

.. code:: MATLAB

    x(:,1)=linspace(0,30,1000);
    y=2.*exp(-0.1*2*pi/5.*x).*sin(sqrt(1-0.1^2)*2*pi/5.*x);
    [xmin,ymin,xmax,ymax]=findextremum(x,y,3,'yes');

    x(:,1)=linspace(0,50,1000);
    y(:,1)=sin(x)+0.1.*rand(1000,1);
    [xmin,ymin,xmax,ymax]=findextremum(x,y,15,'yes');

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
        winlen=3; dispout='no';
    case 3
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
%Checking initial inputs

%Length of the input data
L=length(x(:,1));

%Checking winlen be smaller than the length of input data
if winlen>=L
    winlen=L-1; %Subtracting 1 guarantees winlen stays smaller than L if L is even number
end

%Checking winlen be larger than 3
if winlen<3
    winlen=3;
end

%Checking if winlen is an odd number
if mod(winlen,2)==0
    winlen=winlen+1;
end

%Defining a number of points on each side of each point for defining maximum and minimum
npoints=fix(winlen/2);

%--------------------------------------------------------------------------
%Locating local minimum and maximum points

%Pre-assigning variables
xmin1=zeros(L,1);
ymin1=zeros(L,1);
xmax1=zeros(L,1);
ymax1=zeros(L,1);

%Locating local minimum and maximum points for each sliding window length
for i=1:npoints
    xmin1(i,1)=x(i,1);
    ymin1(i,1)=min(y(1:winlen,1));
    xmax1(i,1)=x(i,1);
    ymax1(i,1)=max(y(1:winlen,1));
end

for i=npoints+1:L-npoints
    xmin1(i,1)=x(i,1);
    ymin1(i,1)=min(y(i-npoints:i+npoints,1));
    xmax1(i,1)=x(i,1);
    ymax1(i,1)=max(y(i-npoints:i+npoints,1));
end

for i=L-npoints+1:L
    xmin1(i,1)=x(i,1);
    ymin1(i,1)=min(y(L-winlen+1:end,1));
    xmax1(i,1)=x(i,1);
    ymax1(i,1)=max(y(L-winlen+1:end,1));
end

%Locating local minimum and maximum points
n=1;
m=1;
for i=1:L
    
    if ymin1(i,1)==y(i,1)
        xmin(m,1)=xmin1(i,1);
        ymin(m,1)=ymin1(i,1);
        m=m+1;
    end
    
    if ymax1(i,1)==y(i,1)
        xmax(n,1)=xmax1(i,1);
        ymax(n,1)=ymax1(i,1);
        n=n+1;
    end
    
end

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1

    %Plotting data
    plot(x,y)
    hold on
    scatter(xmax,ymax,'*r')
    scatter(xmin,ymin,'m')
    xlabel('x data')
    ylabel('y data')
    legend('Data','Maximum','Minimum')

end
    
%--------------------------------------------------------------------------
