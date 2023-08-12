def intersectgc(lat1, lon1, lat2, lon2, lat3, lon3, lat4, lon4, dispout='no'):
    """
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

    scientimate.intersectgc
    =======================

    .. code:: python

        latint1, lonint1, latint2, lonint2, latintn, lonintn, latints, lonints, latinte, loninte, latintw, lonintw = scientimate.intersectgc(lat1, lon1, lat2, lon2, lat3, lon3, lat4, lon4, dispout='no')

    Description
    -----------

    Find intersection point between two line segments (line edges) on Great Circle

    Inputs
    ------

lat1
                                :Latitude (y) of start point of first segment as (lon1,lat1) in (Degree)
lon1
                                :Longitude (x) of start point of first segment as (lon1,lat1) in (Degree)
lat2
                                :Latitude (y) of end point of first segment as (lon2,lat2) in (Degree)
lon2
                                :Longitude (x) of end point of first segment as (lon2,lat2) in (Degree)
                                :First segment: p1(lon1,lat1) to p2(lon2,lat2)
lat3
                                :Latitude (y) of start point of second segment as (lon3,lat3) in (Degree)
lon3
                                :Longitude (x) of start point of second segment as (lon3,lat3) in (Degree)
lat4
                                :Latitude (y) of end point of second segment as (lon4,lat4) in (Degree)
lon4
                                :Longitude (x) of end point of second segment as (lon4,lat4) in (Degree)
                                :Second segment: p3(lon3,lat3) to p4(lon4,lat4)
dispout='no'
                                :Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

latint1
                                :Latitude of first intersection point between two segments in (Degree)
lonint1
                                :Longitude of first intersection point between two segments in (Degree)
latint2
                                :Latitude of second intersection point between two segments in (Degree)
lonint2
                                :Longitude of second intersection point between two segments in (Degree)
                                :Two path on sphere have two intesections on opposite side of sphere 
                                :(latint1,lonint1) and (latint2,lonint2) are antipodes of each other
latintn
                                :Latitude of intersection point between two segments in Northern Hemisphere in (Degree)
lonintn
                                :Longitude of intersection point between two segments in Northern Hemisphere in (Degree)
latints
                                :Latitude of intersection point between two segments in Southern Hemisphere in (Degree)
lonints
                                :Longitude of intersection point between two segments in Southern Hemisphere in (Degree)
latinte
                                :Latitude of intersection point between two segments in Eastern Hemisphere in (Degree)
loninte
                                :Longitude of intersection point between two segments in Eastern Hemisphere in (Degree)
latintw
                                :Latitude of intersection point between two segments in Western Hemisphere in (Degree)
lonintw
                                :Longitude of intersection point between two segments in Western Hemisphere in (Degree)

    Examples
    --------

    .. code:: python

        import scientimate as sm

        #Segment 1:
        lat1=-90
        lon1=0
        lat2=90 
        lon2=0
        #Segment 2:
        lat3=0
        lon3=-90
        lat4=0
        lon4=90
        latint1,lonint1,latint2,lonint2,latintn,lonintn,latints,lonints,latinte,loninte,latintw,lonintw=sm.intersectgc(lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4,'yes')

        #Segment 1:
        lat1=29.079710
        lon1=-90.457758
        lat2=29.344989
        lon2=-90.476782
        #Segment 2:
        lat3=29.190111
        lon3=-90.598023
        lat4=29.206355
        lon4=-90.337782
        latint1,lonint1,latint2,lonint2,latintn,lonintn,latints,lonints,latinte,loninte,latintw,lonintw=sm.intersectgc(lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4,'yes')
        #latint1 = -29.198
        #lonint1 =  89.534
        #latint2 =  29.198
        #lonint2 = -90.466

        #Segment 1:
        lat1=47.94713
        lon1=-131.211073
        lat2=24.207076
        lon2=-83.815088
        #Segment 2:
        lat3=28.257645
        lon3=-95.964404
        lat4=28.343359
        lon4=-42.233815
        latint1,lonint1,latint2,lonint2,latintn,lonintn,latints,lonints,latinte,loninte,latintw,lonintw=sm.intersectgc(lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4,'yes')
        #latint1 =  29.396
        #lonint1 = -89.930
        #latint2 = -29.396
        #lonint2 =  90.070

        #Segment 1:
        lat1=[29.079710,47.94713]
        lon1=[-90.457758,-131.211073]
        lat2=[29.344989,24.207076]
        lon2=[-90.476782,-83.815088]
        #Segment 2:
        lat3=[29.190111,28.257645]
        lon3=[-90.598023,-95.964404]
        lat4=[29.206355,28.343359]
        lon4=[-90.337782,-42.233815]
        latint1,lonint1,latint2,lonint2,latintn,lonintn,latints,lonints,latinte,loninte,latintw,lonintw=sm.intersectgc(lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4,'yes')

    References
    ----------

    | https://en.wikipedia.org/wiki/Geographic_coordinate_conversion
    | https://en.wikipedia.org/wiki/Geographic_coordinate_system
    | https://en.wikipedia.org/wiki/Earth_radius
    | https://stackoverflow.com/questions/29465468/python-intersection-point-of-two-great-circles-lat-long
    | https://stackoverflow.com/questions/2954337/great-circle-rhumb-line-intersection
    | http://www.codeguru.com/cpp/cpp/algorithms/article.php/c5115/Geographic-Distance-and-Azimuth-Calculations.htm

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
    
    lat1=type2numpy(lat1)
    lon1=type2numpy(lon1)
    lat2=type2numpy(lat2)
    lon2=type2numpy(lon2)
    lat3=type2numpy(lat3)
    lon3=type2numpy(lon3)
    lat4=type2numpy(lat4)
    lon4=type2numpy(lon4)

    #--------------------------------------------------------------------------
    #Find location of intersection

    #Converting to radian
    lat1rad=np.deg2rad(lat1)
    lon1rad=np.deg2rad(lon1)
    lat2rad=np.deg2rad(lat2)
    lon2rad=np.deg2rad(lon2)
    lat3rad=np.deg2rad(lat3)
    lon3rad=np.deg2rad(lon3)
    lat4rad=np.deg2rad(lat4)
    lon4rad=np.deg2rad(lon4)
    
    #Earth radius
    a=63781370 #Equatorial radius (m)
    b=63567523 #Polar radius (m)
    
    #Distance from the surface to the z axis along the ellipsoid normal (transverse radius of curvature)
    N1=a**2/np.sqrt(a**2*(np.cos(lat1rad))+b**2*(np.sin(lat1rad))**2)
    N2=a**2/np.sqrt(a**2*(np.cos(lat2rad))+b**2*(np.sin(lat2rad))**2)
    N3=a**2/np.sqrt(a**2*(np.cos(lat3rad))+b**2*(np.sin(lat3rad))**2)
    N4=a**2/np.sqrt(a**2*(np.cos(lat4rad))+b**2*(np.sin(lat4rad))**2)
    
    #Height from earth surface
    h=0 
    
    #Converting geodetic (lat,lon,h) to ECEF (x,y,z) coordinates
    x1=(N1+h)*np.cos(lat1rad)*np.cos(lon1rad)
    y1=(N1+h)*np.cos(lat1rad)*np.sin(lon1rad)
    z1=(N1+h)*np.sin(lat1rad)
    x2=(N2+h)*np.cos(lat2rad)*np.cos(lon2rad)
    y2=(N2+h)*np.cos(lat2rad)*np.sin(lon2rad)
    z2=(N2+h)*np.sin(lat2rad)
    x3=(N3+h)*np.cos(lat3rad)*np.cos(lon3rad)
    y3=(N3+h)*np.cos(lat3rad)*np.sin(lon3rad)
    z3=(N3+h)*np.sin(lat3rad)
    x4=(N4+h)*np.cos(lat4rad)*np.cos(lon4rad)
    y4=(N4+h)*np.cos(lat4rad)*np.sin(lon4rad)
    z4=(N4+h)*np.sin(lat4rad)
    
    #Normal planes that contain great circles
    if len(lat1)==1:
        NP1=np.cross(np.concatenate((x1,y1,z1)),np.concatenate((x2,y2,z2)))
        NP2=np.cross(np.concatenate((x3,y3,z3)),np.concatenate((x4,y4,z4)))
    else:
        NP1=np.cross(np.column_stack((x1,y1,z1)),np.column_stack((x2,y2,z2)))
        NP2=np.cross(np.column_stack((x3,y3,z3)),np.column_stack((x4,y4,z4)))
    
    #Line of intersection between two planes
    L=np.cross(NP1,NP2)
    
    #Intersection points
    X1=np.zeros((len(lat1),3)) #Pre-assigning vector
    if len(lat1)==1:
        X1=L/np.sqrt(L[0]**2+L[1]**2+L[2]**2) #Normalizing
    else:
        for i in range(0,len(lat1),1):
            X1[i,:]=L[i,:]/np.sqrt(L[i,0]**2+L[i,1]**2+L[i,2]**2) #Normalizing

    X2=-X1
    
    #Two path on sphere have two intesections on opposite side of sphere which are antipodes of each other
    if len(lat1)==1:
        latint1=np.rad2deg(np.arcsin(X1[2]))
        lonint1=np.rad2deg(np.arctan2(X1[1],X1[0]))
        latint2=np.rad2deg(np.arcsin(X2[2]))
        lonint2=np.rad2deg(np.arctan2(X2[1],X2[0]))
    else:
        latint1=np.rad2deg(np.arcsin(X1[:,2]))
        lonint1=np.rad2deg(np.arctan2(X1[:,1],X1[:,0]))
        latint2=np.rad2deg(np.arcsin(X2[:,2]))
        lonint2=np.rad2deg(np.arctan2(X2[:,1],X2[:,0]))

    #--------------------------------------------------------------------------
    #Intersection points based on their hemisphere locations

    #Checking hemisphere location of intersections
    SegInsctChkN1=((latint1>=0) & (latint1<=90)) #Checking if latint1 is in Northern Hemisphere
    SegInsctChkS1=((latint1>=-90) & (latint1<0)) #Checking if latint1 is in Southern Hemisphere
    SegInsctChkE1=((lonint1>=0) & (lonint1<=180)) #Checking if lonint1 is in Eastern Hemisphere
    SegInsctChkW1=((lonint1>=-180) & (lonint1<0)) #Checking if lonint1 is in Western Hemisphere
    
    SegInsctChkN2=((latint2>=0) & (latint2<=90)) #Checking if latint2 is in Northern Hemisphere
    SegInsctChkS2=((latint2>=-90) & (latint2<0)) #Checking if latint2 is in Southern Hemisphere
    SegInsctChkE2=((lonint2>=0) & (lonint2<=180)) #Checking if lonint2 is in Eastern Hemisphere
    SegInsctChkW2=((lonint2>=-180) & (lonint2<0)) #Checking if lonint2 is in Western Hemisphere
    
    #Pre-assigning vector
    latintn=np.zeros(len(lat1)) #Pre-assigning vector
    lonintn=np.zeros(len(lat1)) #Pre-assigning vector
    latints=np.zeros(len(lat1)) #Pre-assigning vector
    lonints=np.zeros(len(lat1)) #Pre-assigning vector
    latinte=np.zeros(len(lat1)) #Pre-assigning vector
    loninte=np.zeros(len(lat1)) #Pre-assigning vector
    latintw=np.zeros(len(lat1)) #Pre-assigning vector
    lonintw=np.zeros(len(lat1)) #Pre-assigning vector
    
    #Intersection points based on their locations
    if len(lat1)==1:
        #Intersection points in Northern Hemisphere
        if (SegInsctChkN1==True): latintn=latint1
        if (SegInsctChkN2==True): latintn=latint2
        if (SegInsctChkN1==True): lonintn=lonint1
        if (SegInsctChkN2==True): lonintn=lonint2
        
        #Intersection points in Southern Hemisphere
        if (SegInsctChkS1==True): latints=latint1
        if (SegInsctChkS2==True): latints=latint2
        if (SegInsctChkS1==True): lonints=lonint1
        if (SegInsctChkS2==True): lonints=lonint2
        
        #Intersection points in Eastern Hemisphere
        if (SegInsctChkE1==True): latinte=latint1
        if (SegInsctChkE2==True): latinte=latint2
        if (SegInsctChkE1==True): loninte=lonint1
        if (SegInsctChkE2==True): loninte=lonint2
        
        #Intersection points in Western Hemisphere
        if (SegInsctChkW1==True): latintw=latint1
        if (SegInsctChkW2==True): latintw=latint2
        if (SegInsctChkW1==True): lonintw=lonint1
        if (SegInsctChkW2==True): lonintw=lonint2
        
    else:
        #Intersection points based on their locations
        #Intersection points in Northern Hemisphere
        latintn[SegInsctChkN1==True]=latint1[SegInsctChkN1==True]
        latintn[SegInsctChkN2==True]=latint2[SegInsctChkN2==True]
        lonintn[SegInsctChkN1==True]=lonint1[SegInsctChkN1==True]
        lonintn[SegInsctChkN2==True]=lonint2[SegInsctChkN2==True]
        
        #Intersection points in Southern Hemisphere
        latints[SegInsctChkS1==True]=latint1[SegInsctChkS1==True]
        latints[SegInsctChkS2==True]=latint2[SegInsctChkS2==True]
        lonints[SegInsctChkS1==True]=lonint1[SegInsctChkS1==True]
        lonints[SegInsctChkS2==True]=lonint2[SegInsctChkS2==True]
        
        #Intersection points in Eastern Hemisphere
        latinte[SegInsctChkE1==True]=latint1[SegInsctChkE1==True]
        latinte[SegInsctChkE2==True]=latint2[SegInsctChkE2==True]
        loninte[SegInsctChkE1==True]=lonint1[SegInsctChkE1==True]
        loninte[SegInsctChkE2==True]=lonint2[SegInsctChkE2==True]
        
        #Intersection points in Western Hemisphere
        latintw[SegInsctChkW1==True]=latint1[SegInsctChkW1==True]
        latintw[SegInsctChkW2==True]=latint2[SegInsctChkW2==True]
        lonintw[SegInsctChkW1==True]=lonint1[SegInsctChkW1==True]
        lonintw[SegInsctChkW2==True]=lonint2[SegInsctChkW2==True]

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':

#        if len(lat1)==1:
#            plt.plot([lon1,lon2],[lat1,lat2],label='First Segment')
#            plt.plot([lon3,lon4],[lat3,lat4],label='Second Segment')
#        else:
        for i in range(0,len(lat1),1):
            plt.plot([lon1[i],lon2[i]],[lat1[i],lat2[i]],label='First Segment')
            plt.plot([lon3[i],lon4[i]],[lat3[i],lat4[i]],label='Second Segment')
        
        plt.scatter(lonint1,latint1,label='Intersection Point')
        plt.scatter(lonint2,latint2,label='Intersection Point')
    
        plt.plot([-180,180],[0,0],'k')
        plt.plot([0,0],[-90,90],'k')
    
        plt.xlabel('Longitude (Degree)')
        plt.ylabel('Latitude (Degree)')
        plt.xlim(-180,180)
        plt.ylim(-90,90)
        #plt.legend()
    
    #--------------------------------------------------------------------------
    #Outputs
    return latint1, lonint1, latint2, lonint2, latintn, lonintn, latints, lonints, latinte, loninte, latintw, lonintw

    #--------------------------------------------------------------------------
