def waypointsgc(lat1, lon1, lat2, lon2, nmidwaypoints=0, CalcMethod='haversine', R=6371000, dispout='no'):
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

    scientimate.waypointsgc
    =======================

    .. code:: python

        latwaypoints, lonwaypoints, arclenwaypoints, azimuthwaypoints, arclen, azimuthdir, metedir = scientimate.waypointsgc(lat1, lon1, lat2, lon2, nmidwaypoints=0, CalcMethod='haversine', R=6371000, dispout='no')

    Description
    -----------

    Generate (Latitude,Longitude) points between two (Latitude,Longitude) points using Great Circle

    Inputs
    ------

    lat1
        Latitude (y) of start point (first point) in (Degree)
    lon1
        Longitude (x) of start point (first point) in (Degree)
    lat2
        Latitude (y) of end point (last point) in (Degree)
    lon2
        Longitude (x) of end point (last point) in (Degree)
    nmidwaypoints=1
        | Number of middle waypoints
        | nmidwaypoints=1 : 1 point between first point and last point is generated
    CalcMethod='haversine'
        | Distance calculation method 
        | 'cos': Spherical law of cosines
        | 'haversine': Haversine formula
        | 'vincenty': Vincenty formula, Vincenty (1975)
    R=6371000
        Earth radius in (m), mean earth radius=6371000 m
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

latwaypoints
    | Latitude of waypoints along a line from start point to end point in (Degree)
    | lat1 and lat2 are included, total number of points are nmidwaypoints+2
lonwaypoints
    | Longitude of waypoints along a line from start point to end point in (Degree)
    | lon1 and lon2 are included, total number of points are nmidwaypoints+2
arclenwaypoints
    Distance of waypoints from start point in (m)
azimuthwaypoints
    Azimuth (bearing) of waypoints from start point in (Degree)
arclen
    Total distance from start point to end point in (m)
azimuthdir
    | Azimuth (bearing or compass direction) from start point to end point in (Degree)
    | 0 (degree): toward North, 90 (degree): toward East, 180 (degree): toward South, 270 (degree): toward West 
metedir
    | Meteorological direction from start point to end point in (Degree)
    | 0 (degree): from North, 90 (degree): from East, 180 (degree): from South, 270 (degree): from West 

    Examples
    --------

    .. code:: python

        import scientimate as sm

        lat1=29.5 #First point 
        lon1=-89.4 #First point 
        lat2=29.7 #last point
        lon2=-89.4 #last point
        latwaypoints,lonwaypoints,arclenwaypoints,azimuthwaypoints,arclen,azimuthdir,metedir=sm.waypointsgc(lat1,lon1,lat2,lon2)

        lat1=29 #First point 
        lon1=-89 #First point 
        lat2=30 #Last point
        lon2=-90 #Last point
        latwaypoints,lonwaypoints,arclenwaypoints,azimuthwaypoints,arclen,azimuthdir,metedir=sm.waypointsgc(lat1,lon1,lat2,lon2,3,'haversine',6371000,'yes')

    References
    ----------

    Vincenty, T. (1975). 
    Direct and inverse solutions of geodesics on the ellipsoid with application of nested equations. 
    Survey review, 23(176), 88-93.

    | http://www.movable-type.co.uk/scripts/latlong.html
    | https://en.wikipedia.org/wiki/Great-circle_distance
    | https://en.wikipedia.org/wiki/Great-circle_navigation
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

    #--------------------------------------------------------------------------
    #Calculating distace between start and end of the line (total length)

    #Converting to radian
    lat1rad=np.deg2rad(lat1)
    lon1rad=np.deg2rad(lon1)
    lat2rad=np.deg2rad(lat2)
    lon2rad=np.deg2rad(lon2)
    
    deltalatrad21=lat2rad-lat1rad
    deltalonrad21=lon2rad-lon1rad
    
    #Calculating distance using spherical law of cosines
    if CalcMethod=='cos':
    
        deltasigma=np.arccos(np.sin(lat1rad)*np.sin(lat2rad)+np.cos(lat1rad)*np.cos(lat2rad)*np.cos(deltalonrad21)) #Central angle
        arclen=R*deltasigma #Total distance of the line 
    
    #Calculating distance using Haversine formula
    elif CalcMethod=='haversine':
    
        a=(np.sin(deltalatrad21/2))**2+np.cos(lat1rad)*np.cos(lat2rad)*(np.sin(deltalonrad21/2))**2 #Angular distance in radians
        #c=2*np.arctan2(np.sqrt(a),np.sqrt(1-a)) #Square of half the chord length between the points (central angle)
        deltasigma=2*np.arctan2(np.sqrt(a),np.sqrt(1-a)) #Central angle
        arclen=R*deltasigma #Total distance of the line 
    
    #Calculating distance using Vincenty formula
    elif CalcMethod=='vincenty':
    
        deltasigma=np.arctan2(np.sqrt((np.cos(lat2rad)*np.sin(deltalonrad21))**2+(np.cos(lat1rad)*np.sin(lat2rad)-np.sin(lat1rad)*np.cos(lat2rad)*np.cos(deltalonrad21))**2),np.sin(lat1rad)*np.sin(lat2rad)+np.cos(lat1rad)*np.cos(lat2rad)*np.cos(deltalonrad21)) #Central angle
        arclen=R*deltasigma #Total distance of the line 
    

    #--------------------------------------------------------------------------
    #Calculating azimuth (bearing) between start and end of the line

    azimuthdirrad=np.arctan2(np.sin(deltalonrad21)*np.cos(lat2rad),np.cos(lat1rad)*np.sin(lat2rad)-np.sin(lat1rad)*np.cos(lat2rad)*np.cos(deltalonrad21)) #Azimuth (bearing) in radian
    azimuthdir=np.rad2deg(azimuthdirrad) #Azimuth (bearing) in degree
    
    #Add 360 to all numbers to have them all positive
    #Use mod(360) to take care of the ones larger than 360
    azimuthdir=((azimuthdir+360)%360) 

    #--------------------------------------------------------------------------
    #Generating lat and lon between (lon1,lat1) and (lon2,lat2)

    #Number of segments are nmidwaypoints+1 and number of ponits are nmidwaypoints+2
    #nmidwaypoints=0 means there is no point between start and end points, latpoints=[lat1;lat2]
    #latpoints=np.linspace(lat1,lat2,nmidwaypoints+2) #Latitude of points from start to end
    #latpointsrad=np.deg2rad(latpoints) #Converting to radian
    
    #Generating lat and lon between (lon1,lat1) and (lon2,lat2), first method
    #di(:,1)=linspace(0,d,nmidwaypoints+2); #Distance of each point from first point
    #deltasigmai=di/R #Central angle for each point from first point, deltasigma=d/R
    deltasigmai=np.linspace(0,deltasigma,nmidwaypoints+2) #Central angle for each point from first point, deltasigma=d/R
    latirad=np.arcsin(np.sin(lat1rad)*np.cos(deltasigmai)+np.cos(lat1rad)*np.sin(deltasigmai)*np.cos(azimuthdirrad))
    lonirad=lon1rad+np.arctan2(np.sin(azimuthdirrad)*np.sin(deltasigmai)*np.cos(lat1rad),np.cos(deltasigmai)-np.sin(lat1rad)*np.sin(latirad))
    
    #Generating lat and lon between (lon1,lat1) and (lon2,lat2), second method
    #frac=np.linspace(0,1,nmidwaypoints+2) #Fraction along great circle route, frac=0 for the first point 1 and frac=1 for the last point
    #A=np.sin((1-frac)*deltasigma)/np.sin(deltasigma)
    #B=np.sin(frac*deltasigma)/np.sin(deltasigma)
    #x=A*np.cos(lat1rad)*np.cos(lon1rad)+B*np.cos(lat2rad)*np.cos(lon2rad)
    #y=A*np.cos(lat1rad)*np.sin(lon1rad)+B*np.cos(lat2rad)*np.sin(lon2rad)
    #z=A*np.sin(lat1rad)+B*np.sin(lat2rad)
    #latirad=np.arctan2(z,np.sqrt(x**2+y**2))
    #lonirad=np.arctan2(y,x)
    
    #Make sure first and last points are correct
    #latirad[0]=lat1rad
    #lonirad[0]=lon1rad
    #latirad[-1]=lat1rad
    #lonirad[-1]=lon1rad
    
    #Converting to degree
    latideg=np.rad2deg(latirad)
    lonideg=np.rad2deg(lonirad)

    #--------------------------------------------------------------------------
    #Calculating distance and azimuth (bearing) for each waypoints respect to the first point

    #Calculating distance using spherical law of cosines
    if CalcMethod=='cos':
    
        deltalatirad=latirad-latirad[0]
        deltalonirad=lonirad-lonirad[0]
        deltasigmai=np.arccos(np.sin(latirad[0])*np.sin(latirad)+np.cos(latirad[0])*np.cos(latirad)*np.cos(deltalonirad)) #Central angle
        arcleni=R*deltasigmai #Distance of each point from first point 
        azimuthdirirad=np.arctan2(np.sin(deltalonirad)*np.cos(latirad),np.cos(latirad[0])*np.sin(latirad)-np.sin(latirad[0])*np.cos(latirad)*np.cos(deltalonirad)) #Azimuth (bearing) in radian
    
    #Calculating distance using Haversine formula
    elif CalcMethod=='haversine':
    
        deltalatirad=latirad-latirad[0]
        deltalonirad=lonirad-lonirad[0]
        ai=(np.sin(deltalatirad/2))**2+np.cos(latirad[0])*np.cos(latirad)*(np.sin(deltalonirad/2))**2 #Angular distance in radians
        deltasigmai=2*np.arctan2(np.sqrt(ai),np.sqrt(1-ai)) #Central angle
        arcleni=R*deltasigmai #Distance of each point from first point
        azimuthdirirad=np.arctan2(np.sin(deltalonirad)*np.cos(latirad),np.cos(latirad[0])*np.sin(latirad)-np.sin(latirad[0])*np.cos(latirad)*np.cos(deltalonirad)) #Azimuth (bearing) in radian
    
    #Calculating distance using Vincenty formula
    elif CalcMethod=='vincenty':
    
        deltalatirad=latirad-latirad[0]
        deltalonirad=lonirad-lonirad[0]
        deltasigmai=np.arctan2(np.sqrt((np.cos(latirad)*np.sin(deltalonirad))**2+(np.cos(latirad[0])*np.sin(latirad)-np.sin(latirad[0])*np.cos(latirad)*np.cos(deltalonirad))**2),np.sin(latirad[0])*np.sin(latirad)+np.cos(latirad[0])*np.cos(latirad)*np.cos(deltalonirad)) #Central angle
        arcleni=R*deltasigmai #Distance of each point from first point 
        azimuthdirirad=np.arctan2(np.sin(deltalonirad)*np.cos(latirad),np.cos(latirad[0])*np.sin(latirad)-np.sin(latirad[0])*np.cos(latirad)*np.cos(deltalonirad)) #Azimuth (bearing) in radian
    
    
    azimuthdiri=np.rad2deg(azimuthdirirad) #Azimuth (bearing) in degree
    
    #Add 360 to all numbers to have them all positive
    #Use mod(360) to take care of the ones larger than 360
    azimuthdiri=((azimuthdiri+360)%360) 

    #--------------------------------------------------------------------------
    #Calculating meteorological direction
    metedir=azimuthdir+180 #Convert Azimuth (bearing or compass) angle to meteorological direction
    
    #Add 360 to all numbers to have them all positive
    #Use mod(360) to take care of the ones larger than 360
    metedir=((metedir+360)%360) 

    #--------------------------------------------------------------------------
    #Assigning outputs

    latwaypoints=latideg.copy()
    lonwaypoints=lonideg.copy()
    arclenwaypoints=arcleni.copy()
    azimuthwaypoints=azimuthdiri.copy()

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':

        plt.scatter(lonideg,latideg,label='Intermediate Points')
        plt.scatter(lon1,lat1,label='Start Point')
        plt.scatter(lon2,lat2,label='End Point')
    
        plt.xlabel('Longitude (Degree)')
        plt.ylabel('Latitude (Degree)')
        plt.legend()
    
    #--------------------------------------------------------------------------
    #Outputs
    return latwaypoints, lonwaypoints, arclenwaypoints, azimuthwaypoints, arclen, azimuthdir, metedir

    #--------------------------------------------------------------------------
