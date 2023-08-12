def hurricanewindinflowangle(xgrid, ygrid, Vgrid, xCenter, yCenter, RVmax, inflowdirCalcMethod='bretschneider', distCalcMethod='gc', dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2017-11-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.hurricanewindinflowangle
    ====================================

    .. code:: python

        Vxgrid, Vygrid, thetaVgrid, thetaVgridtan, thetagrid, Rgrid = scientimate.hurricanewindinflowangle(xgrid, ygrid, Vgrid, xCenter, yCenter, RVmax, inflowdirCalcMethod='bretschneider', distCalcMethod='gc', dispout='no')

    Description
    -----------

    Calculate hurricane velocity tangential and inflow angle and inflow velocity in x (East) and y (North) directions

    Inputs
    ------

    xgrid
        | x (longitude) of points which outputs are calculated at as a [M*N] array 
        | xgrid can be a single point or 1d or 2d array 
    ygrid
        | y (latitude) of points which outputs are calculated at as a [M*N] array 
        | ygrid can be a single point or 1d or 2d array
    Vgrid
        | Resultant hurricane wind velocity (Vx^2+Vy^2)^0.5 on (xgrid,ygrid) as a [M*N*L] array in (m/s)
        | L is a number of time steps
        | If only angle values are required, then set Vgrid equal to an arbitary constant such as Vgrid=1
        | For demonstaration purpose, set Vgrid equal to an arbitary constant such as Vgrid=1
    xCenter
        x (longitude) of hurricane center (track) as a [L] array
    yCenter
        y (latitude) of hurricane center (track) as a [L] array
    RVmax
        Distance (radius) from hurricane center to a location of maximum hurricane wind velocity as a [L] array in (m)
    inflowdirCalcMethod='bretschneider'
        | Inflow angle calculation method 
        | 'no': Inflow angle are not calculated, thetaVgrid=thetaVgridtan
        | 'bretschneider': Inflow angle are calculated based on Bretschneider (1972)
        | 'sobey': Inflow angle are calculated based on Sobey et al. (1977)
    distCalcMethod='gc'
        | Distance calculation method 
        | 'cart': Distances are calculated on cartesian coordinate
        | 'gc': Distances are calculated on Great Circle based on Vincenty formula, Vincenty (1975)
        | Earth radius coonsidered as mean earth radius=6371000 m
    dispout='no'
        | Define to display outputs or not
        | 'imagesc': 2 dimensional plot using imagesc or imshow
        | 'pcolor': 2 dimensional plot using pcolor
        | 'contour': 2 dimensional contour plot, number of contour=ncolor
        | 'quiver': 2 dimensional vector plot 
        | 'no': not display 
        | Use dispout='no' if calculation mesh is not 2d array
        | if there is more than one time step, only the last one is plotted
        | if flattendata='yes'; then dispout is set as dispout='no';

    Outputs
    -------

    Vxgrid
        Hurricane wind velocity after applying inflow angle in x (East) direction on defined mesh in (m/s)
    Vygrid
        Hurricane wind velocity after applying inflow angle in y (North) direction on defined mesh in (m/s)
    thetaVgrid
        | Inflow angle (trigonometric direction) of hurricane velocity at each grid point in (Degree)
        | Inflow angle: angle between the inwardly spiraling surface wind 
        |               and the circular isobars around the hurricane center (Boose et al., 2004)
    thetaVgridtan
        | Angle (trigonometric direction) of hurricane velocity at each grid point in (Degree)
        | thetaVgridtan is tangential angle respect to radius. 
        | Note: Outputs has dimension of [M,N,L] where [M,N] is size of the x-y grid and [L] is number of time steps
    thetagrid
        Angle from hurricane center to each point on the grid in (Degree)
    Rgrid
        Distance (radius) from hurricane center to each point on the grid

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        #Creating calculation mesh
        xgrid,ygrid=np.meshgrid(np.linspace(-98,-68,100),np.linspace(16,44,100))

        #Longitude of Hurricane Katrine center at max velocity
        longCenter=-88.6

        #Latitude of Hurricane Katrine center at max velocity
        latCenter=26.3

        #Hurricane Katrina translational velocity (m/s) at max velocity
        Vt=5.18467

        #Hurricane Katrina velocity azimuth (bearing) in (Degree) at max velocity
        VtAzmdir=306.76219

        #Hurricane Katrina 1-min sustained maximum velocity (m/s) at max velocity
        Vmax=76.5
        Vmax=Vmax-Vt #Removing hurricane translation velocity from Vmax
        Vgmax=Vmax/0.8 #Converting surface velocity to gradient velocity

        #Calculating distance using spherical law of cosines
        Rgrid=(np.arccos(np.sin(np.deg2rad(latCenter))*np.sin(np.deg2rad(ygrid))+np.cos(np.deg2rad(latCenter))*np.cos(np.deg2rad(ygrid))*np.cos(np.deg2rad(xgrid)-np.deg2rad(longCenter))))*6371000 #Radius

        #Calculating hurricane velocity at each radius using SLOSH model
        RVmax=32197 #Radius from hurricane center to a location of maximum hurricane wind
        Vgrid=Vgmax*(2*RVmax*Rgrid)/((RVmax)**2+(Rgrid)**2) #Hurricane wind velocity at radius R

        Vxgrid,Vygrid,thetaVgrid,thetaVgridtan,thetagrid,Rgrid=sm.hurricanewindinflowangle(xgrid,ygrid,Vgrid,longCenter,latCenter,RVmax,'bretschneider','gc','quiver')

    References
    ----------

    Data

    * www.nhc.noaa.gov/data/
    * www.nhc.noaa.gov/data/hurdat/hurdat2-format-nencpac.pdf
    * coast.noaa.gov/hurricanes
    * www.aoml.noaa.gov/hrd/data_sub/re_anal.html

    Boose, E. R., Serrano, M. I., & Foster, D. R. (2004). 
    Landscape and regional impacts of hurricanes in Puerto Rico. 
    Ecological Monographs, 74(2), 335-352.

    Bretschneider, C. L. (1972, January). 
    A non-dimensional stationary hurricane wave model. 
    In Offshore Technology Conference. Offshore Technology Conference.

    Phadke, A. C., Martino, C. D., Cheung, K. F., & Houston, S. H. (2003). 
    Modeling of tropical cyclone winds and waves for emergency management. 
    Ocean Engineering, 30(4), 553-578.

    Sobey, R. J., Harper, B. A., & Stark, K. P. (1977). 
    Numerical simulation of tropical cyclone storm surge. 
    James Cook University of North Queensland, Department of Civil & Systems Engineering.

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
    if dispout!='no':
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
    
    xCenter=type2numpy(xCenter)
    yCenter=type2numpy(yCenter)
    RVmax=type2numpy(RVmax)

    #--------------------------------------------------------------------------
    #Pre-assigning array

    #xgrid and ygrid are 2d array
    if np.ndim(xgrid)==2: 
        Mgrid,Ngrid=np.shape(xgrid)

    #xgrid and ygrid are 1d array
    else: 
        Mgrid=(np.shape(xgrid))[0]
        Ngrid=1 #Ngrid=1 when xgrid and ygrid are 1d array

    if len(xCenter)>1:
        Rgrid=np.zeros((Mgrid,Ngrid,len(xCenter))) #Pre-assigning array
        thetagrid=np.zeros((Mgrid,Ngrid,len(xCenter))) #Pre-assigning array
        thetaVgridtan=np.zeros((Mgrid,Ngrid,len(xCenter))) #Pre-assigning array
        thetaVgrid=np.zeros((Mgrid,Ngrid,len(xCenter))) #Pre-assigning array
        Beta=np.zeros((Mgrid,Ngrid,len(xCenter))) #Pre-assigning array
    else:
        Rgrid=np.zeros((Mgrid,Ngrid,1)) #Pre-assigning array
        thetagrid=np.zeros((Mgrid,Ngrid,1)) #Pre-assigning array
        thetaVgridtan=np.zeros((Mgrid,Ngrid,1)) #Pre-assigning array
        thetaVgrid=np.zeros((Mgrid,Ngrid,1)) #Pre-assigning array
        Beta=np.zeros((Mgrid,Ngrid,1)) #Pre-assigning array

    if len(xCenter)==1:
        Vgrid=np.reshape(Vgrid,(Mgrid,Ngrid,1))

    if ((len(xCenter)==1) and (np.size(Vxgrid)!=1)):
        Vgrid=np.reshape(Vgrid,(Mgrid,Ngrid,1))
    elif ((len(xCenter)==1) and (np.size(Vxgrid)==1)):
        Vgrid=np.reshape(Vgrid,(1,1,1))

    #--------------------------------------------------------------------------
    #Calculating distance (radius) from hurricane center to each point

    for i in range(0,len(xCenter),1):
    
        #Calculating distance using cartesian formula
        if distCalcMethod=='cart':
    
            distxy=np.sqrt((xgrid-xCenter[i])**2+(ygrid-yCenter[i])**2) #Calculating distance from (x,y) to (x(1),y(1))
    
            #Calculating angle of the line between start and end points
    
            thetarad=np.arctan2(ygrid-yCenter[i],xgrid-xCenter[i]) #Angle in radian
            theta=np.rad2deg(thetarad) #Angle in degree
    
            #Add 360 to all numbers to have them all positive
            #Use mod(360) to take care of the ones larger than 360
            theta=((theta+360)%360) 
    
            #xgrid and ygrid are 2d array
            if np.ndim(xgrid)==2:
                Rgrid[:,:,i]=distxy.copy() #Distance from hurricane center to each point in (m)
                thetagrid[:,:,i]=theta.copy() #Angle from hurricane center to each point in (degree)

            #xgrid and ygrid are 1d array
            else:
                #Using np.newaxis (or None) to convert shape (Mgrid) into shape (Mgrid,1) when xgrid and ygrid are 1d array
                Rgrid[:,:,i]=distxy[:,np.newaxis].copy() #Distance from hurricane center to each point in (m)
                thetagrid[:,:,i]=theta[:,np.newaxis].copy() #Angle from hurricane center to each point in (degree)
    
        #Calculating distance using Vincenty formula
        elif distCalcMethod=='gc':
    
            #Converting to radian
            lat1rad=np.deg2rad(yCenter[i])
            lon1rad=np.deg2rad(xCenter[i])
            lat2rad=np.deg2rad(ygrid)
            lon2rad=np.deg2rad(xgrid)
    
            deltalatrad21=lat2rad-lat1rad
            deltalonrad21=lon2rad-lon1rad
    
            REarth=6371000 #Earth radius in (m), mean earth radius=6371000 m
            deltasigma=np.arctan2(np.sqrt((np.cos(lat2rad)*np.sin(deltalonrad21))**2+(np.cos(lat1rad)*np.sin(lat2rad)-np.sin(lat1rad)*np.cos(lat2rad)*np.cos(deltalonrad21))**2),np.sin(lat1rad)*np.sin(lat2rad)+np.cos(lat1rad)*np.cos(lat2rad)*np.cos(deltalonrad21)) #Central angle
            arclen=REarth*deltasigma #Total distance of the line 
            distxy=arclen.copy()
    
            #Calculating azimuth (bearing) between start and end of the line
    
            azimuthrad=np.arctan2(np.sin(deltalonrad21)*np.cos(lat2rad),np.cos(lat1rad)*np.sin(lat2rad)-np.sin(lat1rad)*np.cos(lat2rad)*np.cos(deltalonrad21)) #Azimuth (bearing) in radian
            azimuth=np.rad2deg(azimuthrad) #Azimuth (bearing) in degree

            #Add 360 to all numbers to have them all positive
            #Use mod(360) to take care of the ones larger than 360
            azimuth=((azimuth+360)%360) 
    
            #xgrid and ygrid are 2d array
            if np.ndim(xgrid)==2:
                Rgrid[:,:,i]=distxy.copy() #Distance from hurricane center to each point in (m)
                thetagrid[:,:,i]=azimuth.copy() #Angle from hurricane center to each point in (degree)

            #xgrid and ygrid are 1d array
            else:
                #Using np.newaxis (or None) to convert shape (Mgrid) into shape (Mgrid,1) when xgrid and ygrid are 1d array
                Rgrid[:,:,i]=distxy[:,np.newaxis].copy() #Distance from hurricane center to each point in (m)
                thetagrid[:,:,i]=azimuth[:,np.newaxis].copy() #Angle from hurricane center to each point in (degree)
    

    #--------------------------------------------------------------------------
    #Calculating hurricane velocity tangential and inflow angle and inflow velocity in x (East) and y (North) directions
    #Tangential velocity has a right angle respect to radius
    #thetaVgridtan is tangential angle respect to radius 
    #thetaVgrid is inflow angle calculated based on Bretschneider (1972)
    #Inflow angle: angle between the inwardly spiraling surface wind 
    #and the circular isobars around the hurricane center (Boose et al., 2004)

    for i in range(0,len(xCenter),1):
    
        #Calculating inflow angle based on Bretschneider (1972), e.g Phadke et al.,(2003), 
        if inflowdirCalcMethod=='bretschneider':
            #Jeong, C. K., Panchang, V., & Demirbilek, Z. (2012). 
            #Parametric adjustments to the rankine vortex wind model for Gulf of Mexico hurricanes. 
            #Journal of Offshore Mechanics and Arctic Engineering, 134(4), 041102.
            
            #Beta is an inflow angle inward from the tangential flow (Martino et al., 2001)
            Beta[Rgrid<RVmax[i]]=10+(1+Rgrid[Rgrid<RVmax[i]]/RVmax[i])
            Beta[((Rgrid>=RVmax[i]) & (Rgrid<1.2*RVmax[i]))]=20+25*(Rgrid[((Rgrid>=RVmax[i]) & (Rgrid<1.2*RVmax[i]))]/RVmax[i]-1)
            Beta[Rgrid>=1.2*RVmax[i]]=25
    
        #Calculating inflow angle based on Sobey et al. (1977), e.g. Martino et al. (2001)
        elif inflowdirCalcMethod=='sobey':
            #Sobey, R. J., Harper, B. A., & Stark, K. P. (1977). 
            #Numerical simulation of tropical cyclone storm surge. 
            #James Cook University of North Queensland, Department of Civil & Systems Engineering.
    
            #Martino, C. D., Cheung, K. F., Phadke, A. C., & Houston, S. H. (2001, January). 
            #Modeling of hurricane waves in Hawaiian waters. 
            #In The Eleventh International Offshore and Polar Engineering Conference. International Society of Offshore and Polar Engineers.    
            
            #Beta is an inflow angle inward from the tangential flow (Martino et al., 2001)
            Beta[Rgrid<RVmax[i]]=10*(Rgrid[Rgrid<RVmax[i]]/RVmax[i])
            Beta[((Rgrid>=RVmax[i]) & (Rgrid<1.2*RVmax[i]))]=10+75*(Rgrid[((Rgrid>=RVmax[i]) & (Rgrid<1.2*RVmax[i]))]/RVmax[i]-1)
            Beta[Rgrid>=1.2*RVmax[i]]=25
    
        else:
            Beta=0
        
    
    #Calculating velocity vector direction using cartesian formula
    if distCalcMethod=='cart':
    
        #Hurricane is in northern hemisphere
        if np.mean(yCenter)>=0:
            #In northern hemisphere, velocity vector rotate counter clockwise with right angle respect to radius
            thetaVgridtan=thetagrid+90 #Angle of velocity vector in degree
            thetaVgrid=thetaVgridtan+Beta #Inflow angle of velocity vector in degree
    
        #Hurricane is in southern hemisphere
        elif np.mean(yCenter)<0:
            #In southern hemisphere, velocity vector rotate clockwise with right angle respect to radius
            thetaVgridtan=thetagrid-90 #Angle of velocity vector in degree
            thetaVgrid=thetaVgridtan-Beta #Inflow angle of velocity vector in degree
    
    
        #Add 360 to all numbers to have them all positive
        #Use mod(360) to take care of the ones larger than 360
        thetaVgridtan=((thetaVgridtan+360)%360) 
        thetaVgrid=((thetaVgrid+360)%360) 
    
    #Calculating velocity vector direction using great circle formula
    elif distCalcMethod=='gc':
    
        #Converting azimuth (bearing) to trigonometric direction
        thetagrid=-thetagrid+90
    
        #Hurricane is in northern hemisphere
        if np.mean(yCenter)>=0:
            #In northern hemisphere, velocity vector rotate counter clockwise with right angle respect to radius
            thetaVgridtan=thetagrid+90 #Angle of velocity vector in degree
            thetaVgrid=thetaVgridtan+Beta #Inflow angle of velocity vector in degree
    
        #Hurricane is in southern hemisphere
        elif mean(yCenter)<0:
            #In southern hemisphere, velocity vector rotate clockwise with right angle respect to radius
            thetaVgridtan=thetagrid-90 #Angle of velocity vector in degree
            thetaVgrid=thetaVgridtan-Beta #Inflow angle of velocity vector in degree
    
    
        #Add 360 to all numbers to have them all positive
        #Use mod(360) to take care of the ones larger than 360
        thetaVgridtan=((thetaVgridtan+360)%360) 
        thetaVgrid=((thetaVgrid+360)%360) 
    
    
    inflowangle=Beta.copy() #Inflow angle of velocity vector in degree
    Vxgrid=Vgrid*np.cos(np.deg2rad(thetaVgrid)) #Hurricane velocity in x (East) direction
    Vygrid=Vgrid*np.sin(np.deg2rad(thetaVgrid)) #Hurricane velocity in y (North) direction

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='imagesc': #2 dimensional plot using imagesc or imshow
        
        plt.imshow(Vgrid[:,:,-1],extent=[np.nanmin(xgrid),np.nanmax(xgrid),np.nanmin(ygrid),np.nanmax(ygrid)],interpolation=None,cmap=plt.get_cmap(),aspect='auto',origin='lower')
    
    
    elif dispout=='pcolor': #2 dimensional plot using pcolor
    
        plt.pcolormesh(xgrid,ygrid,Vgrid[:,:,-1],cmap=plt.get_cmap())    
    

    elif dispout=='contour': #2 dimensional contour plot
    
        plt.contourf(xgrid,ygrid,Vgrid[:,:,-1],cmap=plt.get_cmap())
    
    
    elif dispout=='quiver': #Surface plot
        
        plt.contour(xgrid,ygrid,Vgrid[:,:,-1],cmap=plt.get_cmap())
        plt.colorbar()
        plt.quiver(xgrid,ygrid,Vxgrid[:,:,-1],Vygrid[:,:,-1])
    
    
    #Setting plot properties
    if ((dispout=='imagesc') or (dispout=='pcolor') or (dispout=='contour') or (dispout=='quiver')):
    
        #Setting axis limits
        #plt.xlim(np.nanmin(xgrid),np.nanmax(xgrid))
        #plt.ylim(np.nanmin(ygrid),np.nanmax(ygrid))
    
        #Plotting colorbar
        if dispout!='quiver':
            plt.colorbar()
    
        #Setting label
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Velocity (m/s)')
        

    #--------------------------------------------------------------------------
    #Converting output array with size of (Mgrid,Ngrid,1) to array with size of (Mgrid,Ngrid)

    if ((len(xCenter)==1) and (np.ndim(xgrid)>1)):
        #xgrid and ygrid are 2d array
        #Converting array with size of (Mgrid,Ngrid,1) to array with size of (Mgrid,Ngrid)
        Vxgrid=np.reshape(Vxgrid,(Mgrid,Ngrid)) #Reshaping array
        Vygrid=np.reshape(Vygrid,(Mgrid,Ngrid)) #Reshaping array
        Vgrid=np.reshape(Vgrid,(Mgrid,Ngrid)) #Reshaping array
        Rgrid=np.reshape(Rgrid,(Mgrid,Ngrid)) #Reshaping array
        thetagrid=np.reshape(thetagrid,(Mgrid,Ngrid)) #Reshaping array
        thetaVgrid=np.reshape(thetaVgrid,(Mgrid,Ngrid)) #Reshaping array
        thetaVgridtan=np.reshape(thetaVgridtan,(Mgrid,Ngrid)) #Reshaping array

    elif ((len(xCenter)==1) and (np.ndim(xgrid)==1)):
        #xgrid and ygrid are 1d array
        #Converting output array with size of (Mgrid,1,1) to array with size of (Mgrid,)
        Vxgrid=np.reshape(Vxgrid,(Mgrid)) #Reshaping array
        Vygrid=np.reshape(Vygrid,(Mgrid)) #Reshaping array
        Vgrid=np.reshape(Vgrid,(Mgrid)) #Reshaping array
        Rgrid=np.reshape(Rgrid,(Mgrid)) #Reshaping array
        thetagrid=np.reshape(thetagrid,(Mgrid)) #Reshaping array
        thetaVgrid=np.reshape(thetaVgrid,(Mgrid)) #Reshaping array
        thetaVgridtan=np.reshape(thetaVgridtan,(Mgrid)) #Reshaping array
    
    #--------------------------------------------------------------------------
    #Output
    return Vxgrid, Vygrid, thetaVgrid, thetaVgridtan, thetagrid, Rgrid

    #--------------------------------------------------------------------------
