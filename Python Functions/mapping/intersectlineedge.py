def intersectlineedge(x1, y1, x2, y2, x3, y3, x4, y4, CalcMethod='vector', dispout='no'):
    """
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

    scientimate.intersectlineedge
    =============================

    .. code:: python

        xintersect, yintersect = scientimate.intersectlineedge(x1, y1, x2, y2, x3, y3, x4, y4, CalcMethod='vector', dispout='no')

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
    CalcMethod='vector'
        | Intersection point calculation method 
        | 'vector': using vector intersection method
        | 'line': using line intersection method
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    xintersect
        x of intersection point between two segments
    yintersect
        y of intersection point between two segments

    Examples
    --------

    .. code:: python

        import scientimate as sm

        #Segment 1:
        x1=1
        y1=1
        x2=5
        y2=5
        #Segment 2:
        x3=5
        y3=1
        x4=1
        y4=5
        xintersect,yintersect=sm.intersectlineedge(x1,y1,x2,y2,x3,y3,x4,y4,'vector','yes')

        #Colinear
        #Segment 1:
        x1=1
        y1=1
        x2=5
        y2=5
        #Segment 2:
        x3=2
        y3=2
        x4=6
        y4=6
        xintersect,yintersect=sm.intersectlineedge(x1,y1,x2,y2,x3,y3,x4,y4,'vector','yes')

        #Parallel
        #Segment 1:
        x1=1
        y1=1
        x2=5
        y2=5
        #Segment 2:
        x3=2
        y3=3
        x4=6
        y4=7
        xintersect,yintersect=sm.intersectlineedge(x1,y1,x2,y2,x3,y3,x4,y4,'vector','yes')

        #Segment 1:
        x1=[1,3]
        y1=[1,4]
        x2=[5,7]
        y2=[5,8]
        #Segment 2:
        x3=[5,4]
        y3=[1,7]
        x4=[1,7]
        y4=[5,2]
        xintersect,yintersect=sm.intersectlineedge(x1,y1,x2,y2,x3,y3,x4,y4)

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
    """

    #--------------------------------------------------------------------------
    #CODE
    #--------------------------------------------------------------------------
    #Import required packages

    import numpy as np
    import scipy as sp
    if dispout=='yes':
        import matplotlib.pyplot as plt 

    #--------------------------------------------------------------------------
    #Convert inputs to numpy array

    #Changing type to numpy array
    def type2numpy(variable):
        if type(variable) is not str:
            if np.size(variable)==1:
                if ((type(variable) is list) or (type(variable) is np.ndarray)):
                    variable=np.array(variable)
                else:
                    variable=np.array([variable])
            elif np.size(variable)>1:
                if (type(variable).__module__)!='numpy':
                    variable=np.array(variable) 
        return variable
    
    x1=type2numpy(x1)
    y1=type2numpy(y1)
    x2=type2numpy(x2)
    y2=type2numpy(y2)
    x3=type2numpy(x3)
    y3=type2numpy(y3)
    x4=type2numpy(x4)
    y4=type2numpy(y4)

    #--------------------------------------------------------------------------
    #Vector intersection method description 
    #http://www.cs.swan.ac.uk/~cssimon/line_intersection.html

    #Segment 1: p1+t12*(p2-p1)
    #Segment 2: p3+t34*(p4-p3)
    #or:
    #Segment 1: x1+t12*(x2-x1) and y1+t12*(y2-y1)
    #Segment 2: x3+t34*(x4-x3) and y3+t34*(y4-y3)

    #Intersection:
    #p1+t12*(p2-p1)=p3+t34*(p4-p3)
    #or:
    #x1+t12*(x2-x1)=x3+t34*(x4-x3)
    #y1+t12*(y2-y1)=y3+t34*(y4-y3)
    #or:
    #(x1-x3)=t34*(x4-x3)-t12*(x2-x1)
    #(y1-y3)=t34*(y4-y3)-t12*(y2-y1)
    #or:
    #[(x4-x3) -(x2-x1)]*[t34]=[x1-x3]
    #[(y4-y3) -(y2-y1)]*[t12]=[y1-y3]
    #or:
    #[t34]=_______________1_______________*[(y1-y2) (x1-x2)][x1-x3]
    #[t12]=(x4-x3)*(y1-y2)-(x1-x2)*(y4-y3)*[(y3-y4) (x4-x3)][y1-y3]
    #or:
    #t12=((y3-y4)*(x1-x3)+(x4-x3)*(y1-y3))/((x4-x3)*(y1-y2)-(x1-x2)*(y4-y3))
    #t34=((y1-y2)*(x1-x3)+(x2-x1)*(y1-y3))/((x4-x3)*(y1-y2)-(x1-x2)*(y4-y3)) 

    #xintersect=x1+t12*(x2-x1);
    #yintersect=y1+t12*(y2-y1);


    #Line intersection method description 
    #https://en.wikipedia.org/wiki/Line-line_intersection

    #xintersect=((x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));
    #yintersect=((x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));

    #Check if intersection point is on both segments

    #--------------------------------------------------------------------------
    #Find location of intersection

    #Calculating intersection point using vector approach
    #http://www.cs.swan.ac.uk/~cssimon/line_intersection.html
    if CalcMethod=='vector':
    
        #Find t12 and t34 which shows a location of intersection
        t12=((y3-y4)*(x1-x3)+(x4-x3)*(y1-y3))/((x4-x3)*(y1-y2)-(x1-x2)*(y4-y3))
        t34=((y1-y2)*(x1-x3)+(x2-x1)*(y1-y3))/((x4-x3)*(y1-y2)-(x1-x2)*(y4-y3)) 
    
        #Checking if Two segments intersects
        #if (t12>=0 & t12<=1) & (t34>=0 & t34<=1) #Two segments intersects
        #if (t12<0 | t12>1) | (t34<0 & t34>1) #Two segments do not intersect but two lines intersects
        #if ((x4-x3)*(y1-y2)-(x1-x2)*(y4-y3))==False #Colinear lines with 0, 1 or more intersections
        SegInsctChk=((t12>=0) & (t12<=1) & (t34>=0) & (t34<=1))
    
        #Find location of intersection
        #Two segments intersects
        xintersect=np.zeros(len(x1)) #Pre-assigning vector
        yintersect=np.zeros(len(x1)) #Pre-assigning vector
        xintersect[SegInsctChk==True]=x1[SegInsctChk==True]+t12[SegInsctChk==True]*(x2[SegInsctChk==True]-x1[SegInsctChk==True])
        yintersect[SegInsctChk==True]=y1[SegInsctChk==True]+t12[SegInsctChk==True]*(y2[SegInsctChk==True]-y1[SegInsctChk==True])
    
        #Two segments do not intersect but two lines intersects
        xintersect[SegInsctChk==False]=np.nan
        yintersect[SegInsctChk==False]=np.nan
    
    #Calculating intersection point using simple approach
    #https://en.wikipedia.org/wiki/Line-line_intersection
    elif CalcMethod=='line':
    
        #Slope of segments
        m12=(y2-y1)/(x2-x1) #Slope of first segment
        m34=(y4-y3)/(x4-x3) #Slope of second segment
        
        #Calculating intersection points
        xintersectall=((x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4))
        yintersectall=((x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4))
    
        #Finding a range in each segment
        minx1x2=np.minimum(x1,x2)
        miny1y2=np.minimum(y1,y2)
        minx3x4=np.minimum(x3,x4)
        miny3y4=np.minimum(y3,y4)
        maxx1x2=np.maximum(x1,x2)
        maxy1y2=np.maximum(y1,y2)
        maxx3x4=np.maximum(x3,x4)
        maxy3y4=np.maximum(y3,y4)
        
        #Checking if Two segments intersects
        SegInsctChkx1x2=((xintersectall>=minx1x2) & (xintersectall<=maxx1x2)) #Checking if xintersect is in [x1,x2] range
        SegInsctChky1y2=((yintersectall>=miny1y2) & (yintersectall<=maxy1y2)) #Checking if yintersect is in [y1,y2] range
        SegInsctChkx3x4=((xintersectall>=minx3x4) & (xintersectall<=maxx3x4)) #Checking if xintersect is in [x3,x4] range
        SegInsctChky3y4=((yintersectall>=miny3y4) & (yintersectall<=maxy3y4)) #Checking if yintersect is in [y3,y4] range
        SegInsctChkm12m34=(m12-m34)>1e-10 #Check if two segments are parallel (Colinear lines)
        
        SegInsctChk=((SegInsctChkx1x2==True) & (SegInsctChky1y2==True) & (SegInsctChkx3x4==True) & (SegInsctChky3y4==True) & (SegInsctChkm12m34==True))
        
        #Two segments intersects
        xintersect=np.zeros(len(x1)) #Pre-assigning vector
        yintersect=np.zeros(len(x1)) #Pre-assigning vector
        xintersect[SegInsctChk==True]=xintersectall[SegInsctChk==True]
        yintersect[SegInsctChk==True]=yintersectall[SegInsctChk==True]
    
        #Two segments do not intersect but two lines intersects
        xintersect[SegInsctChk==False]=np.nan
        yintersect[SegInsctChk==False]=np.nan
        

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':

        for i in range(0,len(x1),1):
            plt.plot([x1[i],x2[i]],[y1[i],y2[i]],label='First Segment')
            plt.plot([x3[i],x4[i]],[y3[i],y4[i]],label='Second Segment')
        
        plt.scatter(xintersect,yintersect,label='Intersection Point')
    
        plt.xlabel('x')
        plt.ylabel('y')
        #plt.legend()
    
    #--------------------------------------------------------------------------
    #Outputs
    return xintersect, yintersect

    #--------------------------------------------------------------------------
