function [xintersect, yintersect] = intersectlineedge(x1, y1, x2, y2, x3, y3, x4, y4, CalcMethod, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-07-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

intersectlineedge
=================

.. code:: MATLAB

    [xintersect, yintersect] = intersectlineedge(x1, y1, x2, y2, x3, y3, x4, y4, CalcMethod, dispout)

Description
-----------

Find intersection point between two line segments (line edges)

Inputs
------

x1
    x of start point of first segment as (x1,y1)
y1
    y of start point of first segment as (x1,y1)
x2
    x of end point of first segment as (x2,y2)
y2
    | y of end point of first segment as (x2,y2)
    | First segment: p1(x1,y1) to p2(x2,y2)
x3
    x of start point of second segment as (x3,y3)
y3
    y of start point of second segment as (x3,y3)
x4
    x of end point of second segment as (x4,y4)
y4
    | y of end point of second segment as (x4,y4)
    | Second segment: p3(x3,y3) to p4(x4,y4)
CalcMethod='vector';
    | Intersection point calculation method 
    | 'vector': using vector intersection method
    | 'line': using line intersection method
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

xintersect
    x of intersection point between two segments
yintersect
    y of intersection point between two segments

Examples
--------

.. code:: MATLAB

    %Segment 1:
    x1=1;
    y1=1;
    x2=5;
    y2=5;
    %Segment 2:
    x3=5;
    y3=1;
    x4=1;
    y4=5;
    [xintersect,yintersect]=intersectlineedge(x1,y1,x2,y2,x3,y3,x4,y4,'vector','yes');

    %Colinear
    %Segment 1:
    x1=1;
    y1=1;
    x2=5;
    y2=5;
    %Segment 2:
    x3=2;
    y3=2;
    x4=6;
    y4=6;
    [xintersect,yintersect]=intersectlineedge(x1,y1,x2,y2,x3,y3,x4,y4,'vector','yes');

    %Parallel
    %Segment 1:
    x1=1;
    y1=1;
    x2=5;
    y2=5;
    %Segment 2:
    x3=2;
    y3=3;
    x4=6;
    y4=7;
    [xintersect,yintersect]=intersectlineedge(x1,y1,x2,y2,x3,y3,x4,y4,'vector','yes');

    %Segment 1:
    x1=[1;3];
    y1=[1;4];
    x2=[5;7];
    y2=[5;8];
    %Segment 2:
    x3=[5;4];
    y3=[1;7];
    x4=[1;7];
    y4=[5;2];
    [xintersect,yintersect]=intersectlineedge(x1,y1,x2,y2,x3,y3,x4,y4);

References
----------

Goldman, R. (1990, August). 
Intersection of two lines in three-space. In Graphics Gems (p. 304). 
Academic Press Professional, Inc..

| http://www.cs.swan.ac.uk/~cssimon/line_intersection.html
| https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
| https://www.cs.hmc.edu/ACM/lectures/intersections.html
| https://www.codeproject.com/Tips/862988/Find-the-Intersection-Point-of-Two-Line-Segments
| https://stackoverflow.com/questions/3838329/how-can-i-check-if-two-segments-intersect
| http://bryceboe.com/2006/10/23/line-segment-intersection-algorithm/
| https://en.wikipedia.org/wiki/Line-line_intersection

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
    case 8
        CalcMethod='vector'; dispout='no';
    case 9
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

if isrow(x3)==1
    x3=x3';
end

if isrow(y3)==1
    y3=y3';
end

if isrow(x4)==1
    x4=x4';
end

if isrow(y4)==1
    y4=y4';
end

%--------------------------------------------------------------------------
%Vector intersection method description 
%http://www.cs.swan.ac.uk/~cssimon/line_intersection.html

%Segment 1: p1+t12*(p2-p1)
%Segment 2: p3+t34*(p4-p3)
%or:
%Segment 1: x1+t12*(x2-x1) and y1+t12*(y2-y1)
%Segment 2: x3+t34*(x4-x3) and y3+t34*(y4-y3)

%Intersection:
%p1+t12*(p2-p1)=p3+t34*(p4-p3)
%or:
%x1+t12*(x2-x1)=x3+t34*(x4-x3)
%y1+t12*(y2-y1)=y3+t34*(y4-y3)
%or:
%(x1-x3)=t34*(x4-x3)-t12*(x2-x1)
%(y1-y3)=t34*(y4-y3)-t12*(y2-y1)
%or:
%[(x4-x3) -(x2-x1)]*[t34]=[x1-x3]
%[(y4-y3) -(y2-y1)]*[t12]=[y1-y3]
%or:
%[t34]=_______________1_______________*[(y1-y2) (x1-x2)][x1-x3]
%[t12]=(x4-x3)*(y1-y2)-(x1-x2)*(y4-y3)*[(y3-y4) (x4-x3)][y1-y3]
%or:
%t12=((y3-y4)*(x1-x3)+(x4-x3)*(y1-y3))/((x4-x3)*(y1-y2)-(x1-x2)*(y4-y3))
%t34=((y1-y2)*(x1-x3)+(x2-x1)*(y1-y3))/((x4-x3)*(y1-y2)-(x1-x2)*(y4-y3)) 

%xintersect=x1+t12*(x2-x1);
%yintersect=y1+t12*(y2-y1);


%Line intersection method description 
%https://en.wikipedia.org/wiki/Line-line_intersection

%xintersect=((x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));
%yintersect=((x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));

%Check if intersection point is on both segments

%--------------------------------------------------------------------------
%Find location of intersection

%Calculating intersection point using vector approach
%http://www.cs.swan.ac.uk/~cssimon/line_intersection.html
if strcmp(CalcMethod,'vector')==1

    %Find t12 and t34 which shows a location of intersection
    t12=((y3-y4).*(x1-x3)+(x4-x3).*(y1-y3))./((x4-x3).*(y1-y2)-(x1-x2).*(y4-y3));
    t34=((y1-y2).*(x1-x3)+(x2-x1).*(y1-y3))./((x4-x3).*(y1-y2)-(x1-x2).*(y4-y3)); 

    %Checking if Two segments intersects
    %if (t12>=0 & t12<=1) & (t34>=0 & t34<=1) %Two segments intersects
    %if (t12<0 | t12>1) | (t34<0 & t34>1) %Two segments do not intersect but two lines intersects
    %if ((x4-x3)*(y1-y2)-(x1-x2)*(y4-y3))==0 %Colinear lines with 0, 1 or more intersections
    SegInsctChk=(t12>=0 & t12<=1) & (t34>=0 & t34<=1);

    %Find location of intersection
    %Two segments intersects
    xintersect=zeros(length(x1(:,1)),1); %Pre-assigning vector
    yintersect=zeros(length(x1(:,1)),1); %Pre-assigning vector
    xintersect(SegInsctChk==1)=x1(SegInsctChk==1)+t12(SegInsctChk==1).*(x2(SegInsctChk==1)-x1(SegInsctChk==1));
    yintersect(SegInsctChk==1)=y1(SegInsctChk==1)+t12(SegInsctChk==1).*(y2(SegInsctChk==1)-y1(SegInsctChk==1));

    %Two segments do not intersect but two lines intersects
    xintersect(SegInsctChk==0)=NaN;
    yintersect(SegInsctChk==0)=NaN;

%Calculating intersection point using simple approach
%https://en.wikipedia.org/wiki/Line-line_intersection
elseif strcmp(CalcMethod,'line')==1

    %Slope of segments
    m12=(y2-y1)./(x2-x1); %Slope of first segment
    m34=(y4-y3)./(x4-x3); %Slope of second segment
    
    %Calculating intersection points
    xintersectall=((x1.*y2-y1.*x2).*(x3-x4)-(x1-x2).*(x3.*y4-y3.*x4))./((x1-x2).*(y3-y4)-(y1-y2).*(x3-x4));
    yintersectall=((x1.*y2-y1.*x2).*(y3-y4)-(y1-y2).*(x3.*y4-y3.*x4))./((x1-x2).*(y3-y4)-(y1-y2).*(x3-x4));

    %Finding a range in each segment
    minx1x2=min(x1,x2);
    miny1y2=min(y1,y2);
    minx3x4=min(x3,x4);
    miny3y4=min(y3,y4);
    maxx1x2=max(x1,x2);
    maxy1y2=max(y1,y2);
    maxx3x4=max(x3,x4);
    maxy3y4=max(y3,y4);
    
    %Checking if Two segments intersects
    SegInsctChkx1x2=(xintersectall>=minx1x2 & xintersectall<=maxx1x2); %Checking if xintersect is in [x1,x2] range
    SegInsctChky1y2=(yintersectall>=miny1y2 & yintersectall<=maxy1y2); %Checking if yintersect is in [y1,y2] range
    SegInsctChkx3x4=(xintersectall>=minx3x4 & xintersectall<=maxx3x4); %Checking if xintersect is in [x3,x4] range
    SegInsctChky3y4=(yintersectall>=miny3y4 & yintersectall<=maxy3y4); %Checking if yintersect is in [y3,y4] range
    SegInsctChkm12m34=(m12-m34)>1e-10; %Check if two segments are parallel (Colinear lines)
    
    SegInsctChk=(SegInsctChkx1x2==1 & SegInsctChky1y2==1 & SegInsctChkx3x4==1 & SegInsctChky3y4==1 & SegInsctChkm12m34==1);
    
    %Two segments intersects
    xintersect=zeros(length(x1(:,1)),1); %Pre-assigning vector
    yintersect=zeros(length(x1(:,1)),1); %Pre-assigning vector
    xintersect(SegInsctChk==1)=xintersectall(SegInsctChk==1);
    yintersect(SegInsctChk==1)=yintersectall(SegInsctChk==1);

    %Two segments do not intersect but two lines intersects
    xintersect(SegInsctChk==0)=NaN;
    yintersect(SegInsctChk==0)=NaN;
    
end

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1

    for i=1:length(x1(:,1))
        plot([x1(i,1),x2(i,1)],[y1(i,1),y2(i,1)])
        hold on
        plot([x3(i,1),x4(i,1)],[y3(i,1),y4(i,1)])
    end
    
    scatter(xintersect,yintersect,'filled')

    xlabel('x')
    ylabel('y')
    %legend('First Segment','Second Segment','Intersection Point')

end

%--------------------------------------------------------------------------
