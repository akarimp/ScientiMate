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
