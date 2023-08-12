def windfetch(xgrid, ygrid, zgrid, x_point, y_point, winddir=0, waterlevel=0, shorelinelevel=0, n_midpoints=100, distCalcMethod='gc', interpMethod='nearest', dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 202-12-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.windfetch
    =====================

    .. code:: python

        wind_fetch, z_mean, z_point, bed_slope_mean, x_fetch, y_fetch, z_fetch, dist_along_fetch = scientimate.windfetch(xgrid, ygrid, zgrid, x_point, y_point, winddir=0, waterlevel=0, shorelinelevel=0, n_midpoints=100, distCalcMethod='gc', interpMethod='nearest', dispout='no')

    Description
    -----------

    Calculate a wind fecth and z (elevation) profile along a path over water for a given 2d x-y domain (map, image, ...)

    Inputs
    ------

    xgrid
        x (longitude) data as a [M*N] array
    ygrid
        y (latitude) data as a [M*N] array
    zgrid
        | z (elevation) data as a [M*N] array
        | z>=0 is considered land, z<0 is considered water
    x_point
        | x (longitude) of the points to find fetches toward those points
        | If x data are longitude, they should be in Degree
    y_point
        | y (latitude) of the points to find fetches toward those points
        | If y data are latitude, they should be in Degree
    winddir=0
        | Meteorological wind direction in (Degree)
        | It represents a direction wind comes from and is measured counter-clockwise from the North
        | 0 (degree): from North, 90 (degree): from East, 180 (degree): from South, 270 (degree): from West
    waterlevel=0
        | Water surface level (Should have the same datum as zgrid)
        | Size of waterlevel should be either 1 or should be equal to size of winddir 
        | If water level is zero (waterlevel=0), it means that the water level is at a level of surface elevation (z=0) in bathymetry data  
        | If water surface level is positive (waterlevel>0), it means that the water level is above zero elevation (z=0) in bathymetry data  
        |     ,therefore water level will be subtracted from bathymetry elevation 
        | If water surface level is negative (waterlevel<0), it means that the water level is below zero elevation (z=0) in bathymetry data  
        |     ,therefore water level will be added to bathymetry elevation 
    shorelinelevel=0
        | Shoreline elevation, 
        | z>=shorelinelevel is considered land, z<shorelinelevel is considered water
    n_midpoints=100
        | Number of middle points generated between first point and last point
        | if n_midpoints=1 then 1 point between first point and last point is generated
    distCalcMethod='gc'
        | Distance calculation method 
        | 'cart': Distances are calculated on cartesian coordinate
        | 'gc': Distances are calculated on Great Circle based on Vincenty formula, Vincenty (1975)
        | Earth radius considered as mean earth radius=6371000 m
    interpMethod='nearest'
        | Interpolation method 
        | 'linear': Use default or 'linear' method to interpolate
        | 'nearest': Use nearest neighbor method to interpolate
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    wind_fetch
        | Total wind fecth distance along a path from shoreline to point (x_point,y_point) with angle of winddir
        | If input data are latitude and longitude in Degree, wind_fetch is in m
    z_mean
        Mean of z (elevation) along a wind fetch
    z_point
        Value of z (elevation) at (x_point,y_point) location
    bed_slope_mean
        Mean of bed slope along a wind fetch
    x_fetch
        x (longitude) along a path at given points (x_fetch,y_fetch)
    y_fetch
        x (latitude) along a path at given points (x_fetch,y_fetch)
    z_fetch
        z (elevation, ...) data along a path at given points (x_fetch,y_fetch)
    dist_along_fetch
        | Distance at each points of (x_fetch,y_fetch) from the domain boundary
        | If input data are latitude and longitude in Degree, dist_along_fetch is in m

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np
        import os

        xgrid,ygrid=np.meshgrid(np.linspace(-5,5,100),np.linspace(-5,5,100))
        R=np.sqrt(xgrid**2+ygrid**2)
        zgrid=np.sin(R)/R
        x_point=-4
        y_point=-2
        winddir=90
        waterlevel=0
        shorelinelevel=0
        n_midpoints=100
        wind_fetch,z_mean,z_point,bed_slope_mean,x_fetch,y_fetch,z_fetch,dist_along_fetch=sm.windfetch(xgrid,ygrid,zgrid,x_point,y_point,winddir,waterlevel,shorelinelevel,n_midpoints,'cart','nearest','yes')

        xgrid,ygrid=np.meshgrid(np.linspace(-5,5,100),np.linspace(-5,5,100))
        zgrid=ygrid*np.sin(xgrid)-xgrid*np.cos(ygrid)-3
        x_point=0
        y_point=0
        winddir=np.arange(0,360+30,30)
        waterlevel=np.ones(np.size(winddir))*0.1
        shorelinelevel=0
        n_midpoints=100
        wind_fetch,z_mean,z_point,bed_slope_mean,x_fetch,y_fetch,z_fetch,dist_along_fetch=sm.windfetch(xgrid,ygrid,zgrid,x_point,y_point,winddir,waterlevel,shorelinelevel,n_midpoints,'cart','nearest','yes')

        #Download Persina Gulf and Gulf of Oman with coordinate exteneds xmin=47, xmax=63, ymin=19, ymax=31 from https://maps.ngdc.noaa.gov/viewers/grid-extract/index.html
        xyzfilename='xyz.csv'; #e.g. xyzfilename='PersianGulf.csv'
        xyzfilelocation=os.getcwd(); #e.g. xyzfilelocation=r'C:/datafolder'
        x,y,z=sm.readxyz(xyzfilename,xyzfilelocation,1,'all');
        xgrid,ygrid,zgrid=sm.interpxyz2grid(x,y,z,100,'points',np.nanmin(x),np.nanmax(x),np.nanmin(y),np.nanmax(y),np.nanmin(z),np.nanmax(z),'all','nearest','no');
        x_point=58.0 #Or x_point=52.0
        y_point=24.5 #Or y_point=26.0
        winddir=[0,40,135,280]
        waterlevel=0
        shorelinelevel=0
        n_midpoints=100
        wind_fetch,z_mean,z_point,bed_slope_mean,x_fetch,y_fetch,z_fetch,dist_along_fetch=sm.windfetch(xgrid,ygrid,zgrid,x_point,y_point,winddir,waterlevel,shorelinelevel,n_midpoints,'gc','nearest','yes');

    References
    ----------

    Vincenty, T. (1975). 
    Direct and inverse solutions of geodesics on the ellipsoid with application of nested equations. 
    Survey review, 23(176), 88-93.

    * http://www.movable-type.co.uk/scripts/latlong.html
    * http://edwilliams.org/gccalc.htm
    * http://edwilliams.org/avform.htm
    * https://www.nhc.noaa.gov/gccalc.shtml

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
    from scipy import interpolate
    import matplotlib as mpl
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
    
    xgrid=type2numpy(xgrid)
    ygrid=type2numpy(ygrid)
    zgrid=type2numpy(zgrid)
    x_point=type2numpy(x_point)
    y_point=type2numpy(y_point)
    winddir=type2numpy(winddir)
    waterlevel=type2numpy(waterlevel)
    shorelinelevel=type2numpy(shorelinelevel)

    #--------------------------------------------------------------------------
    #Generate points on boundary of domain

    #Domain extend
    xmin=np.min(xgrid) #Minimum of x value
    xmax=np.max(xgrid) #Maximum of x value
    ymin=np.min(ygrid) #Minimum of y value
    ymax=np.max(ygrid) #Maximum of y value

    #Points on north, east, south, and west boundary rotates clockwise
    y_boundary_north = np.ones(15)*ymax
    x_boundary_north = np.linspace(xmin,xmax,15)

    y_boundary_east = np.linspace(ymax,ymin,15)
    x_boundary_east = np.ones(15)*xmax

    y_boundary_south = np.ones(15)*ymin
    x_boundary_south = np.linspace(xmax,xmin,15)

    y_boundary_west = np.linspace(ymin,ymax,15)
    x_boundary_west = np.ones(15)*xmin

    #Points on domain boundary, starts from east and rotate counter-clockwise
    #y_boundary = np.concatenate((y_boundary_east[y_boundary_east>=y_point],y_boundary_north[1:],y_boundary_west[1:],y_boundary_south[1:-1],y_boundary_east[y_boundary_east<y_point]))
    #x_boundary = np.concatenate((x_boundary_east[y_boundary_east>=y_point],x_boundary_north[1:],x_boundary_west[1:],x_boundary_south[1:-1],x_boundary_east[y_boundary_east<y_point]))

    #Points on domain boundary, starts from north and rotate clockwise
    y_boundary = np.concatenate((y_boundary_north[x_boundary_north>=x_point],y_boundary_east[1:],y_boundary_south[1:],y_boundary_west[1:-1],y_boundary_north[x_boundary_north<x_point]))
    x_boundary = np.concatenate((x_boundary_north[x_boundary_north>=x_point],x_boundary_east[1:],x_boundary_south[1:],x_boundary_west[1:-1],x_boundary_north[x_boundary_north<x_point]))

    #Add the the last point before the first point and the first point after the last point to have a closed loop
    y_boundary = np.insert(y_boundary, 0, y_boundary[-1])
    x_boundary = np.insert(x_boundary, 0, x_boundary[-1])

    y_boundary = np.append(y_boundary,y_boundary[0])
    x_boundary = np.append(x_boundary,x_boundary[0])

    #Change to column vector
    #y_boundary = y_boundary.T;
    #x_boundary = x_boundary.T;


    #--------------------------------------------------------------------------
    #Generating a line segment that start from (x,y) with a length of diagonallength and a same direction as a wind direction
    #This line starts from a point of interest inside a domain and extend toward outside a domain while intersecting with a domain boundary 

    #Converte meteorological direction to trigonometric direction
    trigdir=-winddir+270
    
    #Converte meteorological direction to azimuth (bearing)
    azimuthdir=winddir-180
    
    #Add 360 to all numbers to have them all positive
    #Use mod(360) to take care of the ones larger than 360
    trigdir=((trigdir+360)%360) 
    azimuthdir=((azimuthdir+360)%360) 
    
    
    #--------------------------------------------------------------------------
    #Generating points on transect line starts from point of interest toward points on boundary 

    #Calculate using cartesian formula
    if distCalcMethod=='cart':
        
        #Generating points from (x1,y1) to (x2,y2)
        if (n_midpoints<1): n_midpoints=1
        x_transect_360=np.linspace(x_point,x_boundary,n_midpoints+2)
        y_transect_360=np.linspace(y_point,y_boundary,n_midpoints+2)
    
    #Calculate using great circle formula
    elif distCalcMethod=='gc':
    
        #Converte to radian
        lat1_rad=np.deg2rad(y_point)
        lon1_rad=np.deg2rad(x_point)
        lat2_rad=np.deg2rad(y_boundary)
        lon2_rad=np.deg2rad(x_boundary)
    
        delta_lat_21_rad=lat2_rad-lat1_rad
        delta_lon_21_rad=lon2_rad-lon1_rad
    
        #Calculating azimuth (bearing) between start and end of the line
        azimuthdir_rad=np.arctan2(np.sin(delta_lon_21_rad)*np.cos(lat2_rad),np.cos(lat1_rad)*np.sin(lat2_rad)-np.sin(lat1_rad)*np.cos(lat2_rad)*np.cos(delta_lon_21_rad)) #Azimuth (bearing) in radian

        R=6371000 #Earth radius in (m), mean earth radius=6371000 m
        deltasigma=np.arctan2(np.sqrt((np.cos(lat2_rad)*np.sin(delta_lon_21_rad))**2+(np.cos(lat1_rad)*np.sin(lat2_rad)-np.sin(lat1_rad)*np.cos(lat2_rad)*np.cos(delta_lon_21_rad))**2),np.sin(lat1_rad)*np.sin(lat2_rad)+np.cos(lat1_rad)*np.cos(lat2_rad)*np.cos(delta_lon_21_rad)) #Central angle
    
        #Generating lat and lon between (lon1,lat1) and (lon2,lat2), first method
        if (n_midpoints<1): n_midpoints=1
        deltasigmai=np.linspace(0,deltasigma,n_midpoints+2) #Central angle for each point from first point, deltasigma=d/R
        lati_rad=np.arcsin(np.sin(lat1_rad)*np.cos(deltasigmai)+np.cos(lat1_rad)*np.sin(deltasigmai)*np.cos(azimuthdir_rad))
        loni_rad=lon1_rad+np.arctan2(np.sin(azimuthdir_rad)*np.sin(deltasigmai)*np.cos(lat1_rad),np.cos(deltasigmai)-np.sin(lat1_rad)*np.sin(lati_rad))
    
        #Converte to degree
        y_transect_360=np.rad2deg(lati_rad)
        x_transect_360=np.rad2deg(loni_rad)
    
    #Transpose array to have shape of (n_point_on_boundary,n_midpoints+2)
    y_transect_360=y_transect_360.T
    x_transect_360=x_transect_360.T

    #--------------------------------------------------------------------------
    #Calculate z (elevation) profile along transect path from boundary domain toward point of interest

    #Interpolating the z (elevation) values for given (x,y) from input data
    if interpMethod=='linear':
        fun=sp.interpolate.LinearNDInterpolator((np.reshape(xgrid,np.size(xgrid)),np.reshape(ygrid,np.size(ygrid))),np.reshape(zgrid,np.size(zgrid)))
        z_transect_360=fun(x_transect_360,y_transect_360)
    
        #Replacing NaN data point resulted from default method with ones from nearest method
        if np.sum(np.isnan(z_transect_360))>0:
            fun=sp.interpolate.NearestNDInterpolator((np.reshape(xgrid,np.size(xgrid)),np.reshape(ygrid,np.size(ygrid))),np.reshape(zgrid,np.size(zgrid)))
            z_nearest=fun(x_transect_360,y_transect_360)
            z_transect_360[np.isnan(z_transect_360)==True]=z_nearest[np.isnan(z_transect_360)==True] 

    
    #Interpolating data into grid using nearest neighbor method
    elif interpMethod=='nearest':
    
        #First method
        fun=sp.interpolate.NearestNDInterpolator((np.reshape(xgrid,np.size(xgrid)),np.reshape(ygrid,np.size(ygrid))),np.reshape(zgrid,np.size(zgrid)))
        z_transect_360=fun(x_transect_360,y_transect_360)
    
    #--------------------------------------------------------------------------
    #Find intersection points between shoreline and transect lines

    #Create land_water array
    z_transect_360_land_water=z_transect_360.copy()
    z_transect_360_land_water[z_transect_360>=shorelinelevel]=1 #All points on land are changed to 1
    z_transect_360_land_water[z_transect_360<shorelinelevel]=-1 #All points in water are changed to -1

    #Pre-assign array
    M, N=np.shape(x_transect_360)
    indx_beforeSHL=np.ones(M, dtype=np.int_)
    indx_afterSHL=np.ones(M, dtype=np.int_)
    x_shorelineintersect_360=np.zeros(M)
    y_shorelineintersect_360=np.zeros(M)

    #Find intersection points between shoreline and transect lines
    for i in range(0,M,1):

        #Check if point of interest (start of transect) is in water
        if z_transect_360_land_water[i,0]==-1:

            #Check if transect pass through both water and land
            if np.any(z_transect_360_land_water[i,:]==1)==1:
                
                #Find two points one before and one after shoreline
                indx_afterSHL[i]=np.nonzero(z_transect_360_land_water[i,:]==1)[0][0]
                indx_beforeSHL[i]=int(indx_afterSHL[i]-1)

                xpoint_beforeSHL=x_transect_360[i,indx_beforeSHL[i]] #Assign x of point before shoreline
                ypoint_beforeSHL=y_transect_360[i,indx_beforeSHL[i]] #Assign y of point before shoreline
                zpoint_beforeSHL=z_transect_360[i,indx_beforeSHL[i]] #Assign z of point before shoreline
                xpoint_afterSHL=x_transect_360[i,indx_afterSHL[i]] #Assign x of point after shoreline
                ypoint_afterSHL=y_transect_360[i,indx_afterSHL[i]] #Assign y of point after shoreline
                zpoint_afterSHL=x_transect_360[i,indx_afterSHL[i]] #Assign z of point after shoreline

                #Calculate intesction point between shoreline and trensect line
                m_Line_y=(zpoint_afterSHL-zpoint_beforeSHL)/(ypoint_afterSHL-ypoint_beforeSHL) #Slope of line between point before and point afer shoreline in latitudinal direction
                m_Line_x=(zpoint_afterSHL-zpoint_beforeSHL)/(xpoint_afterSHL-xpoint_beforeSHL) #Slope of line between point before and point after shoreline in longitudinal direction
                
                x_shorelineintersect_360[i]=(-zpoint_beforeSHL/m_Line_x)+xpoint_beforeSHL #x of intesction point between shoreline and trensect line
                y_shorelineintersect_360[i]=(-zpoint_beforeSHL/m_Line_y)+ypoint_beforeSHL #y of intesction point between shoreline and trensect line
                
                if ((np.isnan(y_shorelineintersect_360[i])==True) or (np.isnan(x_shorelineintersect_360[i])==True) or (np.isinf(y_shorelineintersect_360[i])==True) or (np.isinf(x_shorelineintersect_360[i])==True)):
                    #Assign mean velues to intesction point between shoreline and trensect line
                    x_shorelineintersect_360[i]=(xpoint_beforeSHL+xpoint_afterSHL)/2 #x of intesction point between shoreline and trensect line
                    y_shorelineintersect_360[i]=(ypoint_beforeSHL+ypoint_afterSHL)/2 #y of intesction point between shoreline and trensect line


            #Transect line passes completely over water and wind fetch extend to the boundary
            else:
                indx_afterSHL[i]=M
                indx_beforeSHL[i]=M
                x_shorelineintersect_360[i]=x_transect_360[i,-1] #Assign the farest point (on boundary) x value as x value of intersection points
                y_shorelineintersect_360[i]=y_transect_360[i,-1] #Assign the farest point (on boundary) y value as y value of intersection points


        #Check if point of interest (start of transect) is on land
        elif z_transect_360_land_water[i,0]==1:
            indx_afterSHL[i]=np.nan #Or float('nan')
            indx_beforeSHL[i]=np.nan #Or float('nan')
            x_shorelineintersect_360[i]=np.nan #Assign NaN as x value of intersection points
            y_shorelineintersect_360[i]=np.nan #Assign NaN as y value of intersection points


    #--------------------------------------------------------------------------
    #Generate points on transect line starts from shoreline toward the point of interest

    #Calculate using cartesian formula
    if distCalcMethod=='cart':
        
        #Generating points from (x1,y1) to (x2,y2)
        if (n_midpoints<1): n_midpoints=1
        x_fetch_360=np.linspace(x_shorelineintersect_360,x_point,n_midpoints+2)
        y_fetch_360=np.linspace(y_shorelineintersect_360,y_point,n_midpoints+2)
    
    #Calculate using great circle formula
    elif distCalcMethod=='gc':
    
        #Converte to radian
        lat1_rad=np.deg2rad(y_shorelineintersect_360)
        lon1_rad=np.deg2rad(x_shorelineintersect_360)
        lat2_rad=np.deg2rad(y_point)
        lon2_rad=np.deg2rad(x_point)
    
        delta_lat_21_rad=lat2_rad-lat1_rad
        delta_lon_21_rad=lon2_rad-lon1_rad
    
        #Calculating azimuth (bearing) between start and end of the line
        azimuthdir_rad=np.arctan2(np.sin(delta_lon_21_rad)*np.cos(lat2_rad),np.cos(lat1_rad)*np.sin(lat2_rad)-np.sin(lat1_rad)*np.cos(lat2_rad)*np.cos(delta_lon_21_rad)) #Azimuth (bearing) in radian

        R=6371000 #Earth radius in (m), mean earth radius=6371000 m
        deltasigma=np.arctan2(np.sqrt((np.cos(lat2_rad)*np.sin(delta_lon_21_rad))**2+(np.cos(lat1_rad)*np.sin(lat2_rad)-np.sin(lat1_rad)*np.cos(lat2_rad)*np.cos(delta_lon_21_rad))**2),np.sin(lat1_rad)*np.sin(lat2_rad)+np.cos(lat1_rad)*np.cos(lat2_rad)*np.cos(delta_lon_21_rad)) #Central angle
    
        #Generating lat and lon between (lon1,lat1) and (lon2,lat2), first method
        if (n_midpoints<1): n_midpoints=1
        deltasigmai=np.linspace(0,deltasigma,n_midpoints+2) #Central angle for each point from first point, deltasigma=d/R
        lati_rad=np.arcsin(np.sin(lat1_rad)*np.cos(deltasigmai)+np.cos(lat1_rad)*np.sin(deltasigmai)*np.cos(azimuthdir_rad))
        loni_rad=lon1_rad+np.arctan2(np.sin(azimuthdir_rad)*np.sin(deltasigmai)*np.cos(lat1_rad),np.cos(deltasigmai)-np.sin(lat1_rad)*np.sin(lati_rad))
    
        #Converte to degree
        y_fetch_360=np.rad2deg(lati_rad)
        x_fetch_360=np.rad2deg(loni_rad)
    
    #Transpose array to have shape of (n_point_on_boundary,n_midpoints+2)
    y_fetch_360=y_fetch_360.T
    x_fetch_360=x_fetch_360.T

    #--------------------------------------------------------------------------
    #Calculate distace of wind fetch points starts from shoreline toward point of interest

    #Calculate distance using cartesian formula
    if distCalcMethod=='cart':
    
        dist_xy_360=np.sqrt((x_fetch_360-x_point)**2+(y_fetch_360-y_point)**2) #Calculate distance from (x,y) to (x(1),y(1))
    
    #Calculate distance using Vincenty formula
    elif distCalcMethod=='gc':
    
        #Converte to radian
        lat1_rad=np.deg2rad(y_fetch_360)
        lon1_rad=np.deg2rad(x_fetch_360)
        lat2_rad=np.deg2rad(y_point)
        lon2_rad=np.deg2rad(x_point)
    
        delta_lat_21_rad=lat2_rad-lat1_rad
        delta_lon_21_rad=lon2_rad-lon1_rad
    
        R=6371000 #Earth radius in (m), mean earth radius=6371000 m
        deltasigma=np.arctan2(np.sqrt((np.cos(lat2_rad)*np.sin(delta_lon_21_rad))**2+(np.cos(lat1_rad)*np.sin(lat2_rad)-np.sin(lat1_rad)*np.cos(lat2_rad)*np.cos(delta_lon_21_rad))**2),np.sin(lat1_rad)*np.sin(lat2_rad)+np.cos(lat1_rad)*np.cos(lat2_rad)*np.cos(delta_lon_21_rad)) #Central angle
        arclen=R*deltasigma #Total distance of the line 
        dist_xy_360=arclen.copy()
    
    
    #Total wind fetch distance from shoreline to the point of interest
    wind_fetch_360=dist_xy_360[:,0]
    
    #Wind fetch distance for each point from shoreline
    dist_along_fetch_360=wind_fetch_360[:,np.newaxis]-dist_xy_360

    #--------------------------------------------------------------------------
    #Calculate z (elevation) profile along transect path from shoreline toward point of interest

    #Interpolating the z (elevation) values for given (x,y) from input data
    if interpMethod=='linear':
        fun=sp.interpolate.LinearNDInterpolator((np.reshape(xgrid,np.size(xgrid)),np.reshape(ygrid,np.size(ygrid))),np.reshape(zgrid,np.size(zgrid)))
        z_fetch_360=fun(x_fetch_360,y_fetch_360)
    
        #Replacing NaN data point resulted from default method with ones from nearest method
        if np.sum(np.isnan(z_fetch_360))>0:
            #If point of interest is on land then all elements of x_fetch_360 and x_fetch_360 would be NaN
            if ((np.all(np.isnan(x_fetch_360))==True) or (np.all(np.isnan(y_fetch_360))==True)):
                dummy=1
            else:
                fun=sp.interpolate.NearestNDInterpolator((np.reshape(xgrid,np.size(xgrid)),np.reshape(ygrid,np.size(ygrid))),np.reshape(zgrid,np.size(zgrid)))
                z_nearest=fun(x_fetch_360,y_fetch_360)
                z_fetch_360[np.isnan(z_fetch_360)==True]=z_nearest[np.isnan(z_fetch_360)==True] 

    
    #Interpolating data into grid using nearest neighbor method
    elif interpMethod=='nearest':
    
        #First method
        #If point of interest is on land then all elements of x_fetch_360 and x_fetch_360 would be NaN
        if ((np.all(np.isnan(x_fetch_360))==True) or (np.all(np.isnan(y_fetch_360))==True)):
            z_fetch_360=np.zeros(x_fetch_360.shape)
            z_fetch_360[:]=np.nan
        else:
            fun=sp.interpolate.NearestNDInterpolator((np.reshape(xgrid,np.size(xgrid)),np.reshape(ygrid,np.size(ygrid))),np.reshape(zgrid,np.size(zgrid)))
            z_fetch_360=fun(x_fetch_360,y_fetch_360)
    
    
    z_point=z_fetch_360[0,-1] #Value of z at (x_point,y_point) location

    #--------------------------------------------------------------------------
    #Calculate mean z (elevation)

    #First method
    z_mean_360=np.nanmean(z_fetch_360,axis=1)

    #Second method
    #Calculate zdx=z(i)*dx(i)
    #zdx=(0.5*(z_fetch_360[:,0:-1]+z_fetch_360[:,1:]))*np.diff(dist_along_fetch_360,n=1,axis=1)

    #Calculate mean z (elevation) as z_mean = sum(z(i)*dx(i)) / total_distance
    #z_mean_360=np.sum(zdx,axis=1)/dist_along_fetch_360[:,-1]

    #for i in range(1,M,1):
    #    if dist_along_fetch_360[i,-1]==0:
    #        z_mean_360[i]=z_fetch_360[i,-1]

    #--------------------------------------------------------------------------
    #Calculate mean slope along wind fetch

    #Calculate slope_along_fetch_360=dz(i)/dx(i)
    slope_along_fetch_360=np.diff(z_fetch_360,n=1,axis=1)/np.diff(dist_along_fetch_360,n=1,axis=1)

    #Calculate mean slope
    bed_slope_mean_360=np.nanmean(slope_along_fetch_360,axis=1)

    #--------------------------------------------------------------------------
    #Calculate fetch direction from shoreline toward point of interest

    #Calculate using cartesian formula
    if distCalcMethod=='cart':
        
        #Generating direction from (x1,y1) to (x2,y2)
        delta_x_360 = x_point-x_shorelineintersect_360
        delta_y_360 = y_point-y_shorelineintersect_360
        winddir_360_trig_rad=np.arctan2(delta_y_360,delta_x_360)

        winddir_360_trig=np.rad2deg(winddir_360_trig_rad) #Wind direction in degree

        #Converting trigonometric direction to meteorological direction
        winddir_360=270-winddir_360_trig

        #Add 360 to all numbers to have them all positive
        #Use mod(360) to take care of the ones larger than 360
        winddir_360=((winddir_360+360)%360) 

    #Calculate using great circle formula
    elif distCalcMethod=='gc':

        #Converte to radian
        lat1_rad=np.deg2rad(y_shorelineintersect_360)
        lon1_rad=np.deg2rad(x_shorelineintersect_360)
        lat2_rad=np.deg2rad(y_point)
        lon2_rad=np.deg2rad(x_point)

        delta_lat_21_rad=lat2_rad-lat1_rad
        delta_lon_21_rad=lon2_rad-lon1_rad

        #Calculating azimuth (bearing) between start and end of the line
        azimuthdir_rad=np.arctan2(np.sin(delta_lon_21_rad)*np.cos(lat2_rad),np.cos(lat1_rad)*np.sin(lat2_rad)-np.sin(lat1_rad)*np.cos(lat2_rad)*np.cos(delta_lon_21_rad)) #Azimuth (bearing) in radian
        azimuthdir_360=np.rad2deg(azimuthdir_rad) #Azimuth (bearing) in degree

        #Converting azimuth (bearing) to meteorological direction
        winddir_360=azimuthdir_360+180
        
        #Add 360 to all numbers to have them all positive
        #Use mod(360) to take care of the ones larger than 360
        winddir_360=((winddir_360+360)%360) 
        

    #Subtract 360 degree from the first point and add 360 degree to the last point to have a closed loop
    winddir_360[0]=winddir_360[0]-360
    winddir_360[-1]=winddir_360[-1]+360

    #--------------------------------------------------------------------------
    #Interpolate values for given winddir

    f = sp.interpolate.interp1d(winddir_360,wind_fetch_360)
    wind_fetch = f(winddir)

    f = sp.interpolate.interp1d(winddir_360,z_mean_360)
    z_mean = f(winddir)

    f = sp.interpolate.interp1d(winddir_360,bed_slope_mean_360)
    bed_slope_mean = f(winddir)
    
    f = sp.interpolate.interp1d(winddir_360,x_fetch_360,axis=0)
    x_fetch = f(winddir)
    
    f = sp.interpolate.interp1d(winddir_360,y_fetch_360,axis=0)
    y_fetch = f(winddir)

    f = sp.interpolate.interp1d(winddir_360,z_fetch_360,axis=0)
    z_fetch = f(winddir)

    f = sp.interpolate.interp1d(winddir_360,dist_along_fetch_360,axis=0)
    dist_along_fetch = f(winddir)

    #--------------------------------------------------------------------------
    #Adjust z (elevation) data for water surface elevation

    #For size(waterlevel)=1
    if np.size(waterlevel)==1:
        z_mean = z_mean-waterlevel
        z_point = z_point-waterlevel
        z_fetch = z_fetch-waterlevel

    #For size(waterlevel)>1
    elif np.size(waterlevel)==np.size(winddir):
        z_mean = z_mean-waterlevel
        z_point = z_point-waterlevel
        for i in range(0,len(winddir),1):
            z_fetch[i,:] = z_fetch[i,:]-waterlevel[i]


    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
        
        plt.subplot(2,1,1)
        #Plotting domain
        levels=[np.nanmin(zgrid), 0, -np.nanmin(zgrid)]
        plt.contourf(xgrid,ygrid,zgrid,levels=levels,cmap=plt.get_cmap())

        #Plotting line 
        plt.scatter(x_point,y_point)
        sct3 = plt.scatter(x_shorelineintersect_360,y_shorelineintersect_360,marker='.',label='Shoreline')
        for i in range(0,len(winddir),1):
            #plt.plot([x_point,x_shorelineintersect_360[i]],[y_point,y_shorelineintersect_360[i]],'k')
            #plt.scatter(x_fetch_360[i,:],y_fetch_360[i,:],c='k',marker='*')
            sct1 = plt.scatter(x_fetch[i,-1],y_fetch[i,-1],c='k',marker='*',label='(x-point,y-point)')
            sct2 = plt.scatter(x_fetch[i,0],y_fetch[i,0],c='k',marker='x',label='Start of the Fetch')
            plt.plot([x_fetch[i,0],x_fetch[i,-1]],[y_fetch[i,0],y_fetch[i,-1]],'k')

        plt.xlabel('x')
        plt.ylabel('y')
        plt.legend(handles=[sct1, sct2, sct3])

        plt.subplot(2,1,2)
        #Plotting depth along a line from point 2 to point1
        for i in range(0,len(winddir),1):
            sct1 = plt.scatter(dist_along_fetch[i,-1],z_fetch[i,-1],c='k',marker='*',label='(x-point,y-point)')
            sct2 = plt.scatter(dist_along_fetch[i,0],z_fetch[i,0],c='k',marker='x',label='Start of the Fetch')
            plt.plot(dist_along_fetch[i,:],z_fetch[i,:])

        plt.xlabel('Distance from Shoreline')
        plt.ylabel('z')
        plt.legend(handles=[sct1, sct2])


    #--------------------------------------------------------------------------
    #Outputs
    return wind_fetch, z_mean, z_point, bed_slope_mean, x_fetch, y_fetch, z_fetch, dist_along_fetch

    #--------------------------------------------------------------------------
