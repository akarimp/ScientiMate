def hurricanebackgroundwind(xgrid, ygrid, Vxgrid, Vygrid, xCenter, yCenter, Vt, VtAzmdir, RVmax, backwindCalcMethod='slosh', Rmax=423e3, distCalcMethod='gc', dispout='no'):
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

    scientimate.hurricanebackgroundwind
    ===================================

    .. code:: python

        Vxgridbackwind, Vygridbackwind, Vgridbackwind, Rgrid = scientimate.hurricanebackgroundwind(xgrid, ygrid, Vxgrid, Vygrid, xCenter, yCenter, Vt, VtAzmdir, RVmax, backwindCalcMethod='slosh', Rmax=423e3, distCalcMethod='gc', dispout='no')

    Description
    -----------

    Calculate and add background wind velocity due to hurricane front motion to hurricane rotational wind velocity

    Inputs
    ------

    xgrid
        | x (longitude) of points which outputs are calculated at as a [M*N] array 
        | xgrid can be a single point or 1d or 2d array 
    ygrid
        | y (latitude) of points which outputs are calculated at as a [M*N] array 
        | ygrid can be a single point or 1d or 2d array
    Vxgrid
        | Hurricane wind velocity in x (East) direction as a [M*N*L] array in (m/s)
        | L is a number of time steps
        | If only background wind is required, set Vxgrid=0 and Vygrid=0
    Vygrid
        | Hurricane wind velocity in y (North) direction as a [M*N*L] array in (m/s)
        | L is a number of time steps
        | If only background wind is required, set Vxgrid=0 and Vygrid=0
    xCenter
        x (longitude) of hurricane center (track) as a [L] array
    yCenter
        y (latitude) of hurricane center (track) as a [L] array
    Vt
        Hurricane central translational velocity as a [L] array in (m/s)
    VtAzmdir
        | Hurricane center velocity azimuth (bearing) direction as a [L] array in (Degree)
        | azimuth (bearing) direction which is measured clockwise from the north:
        | 0 (degree): toward North, 90 (degree): toward East, 180 (degree): toward South, 270 (degree): toward West 
    RVmax
        Distance (radius) from hurricane center to a location of maximum hurricane wind velocity as a [L] array in (m)
    backwindCalcMethod='slosh'
        | Calculation method for adding background wind velocity due to hurricane motion
        | background wind velocity is added to points with R<=Rmax, e.g. Chavas & Emanuel (2010), Lin & Chavas (2012)
        | 'no': background wind velocity is not added to hurricane velocities
        | 'constant': background wind velocity is added as constant value to hurricane velocities
        | 'slosh': Background wind velocity is calculated and added using SLOSH model by Jelesnianski et al. (1992)
        | 'lin': Background wind velocity is calculated and added based on Lin & Chavas (2012)
    Rmax=423e3
        | Maximum radius of hurricane from hurricane center in (m)
        | Outputs for points with R>Rmax is set to zero
        | Median values of Rmax is 423 km (e.g. Chavas & Emanuel 2010; Lin & Chavas, 2012)
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

    Vxgridbackwind
        Hurricane wind velocity with background wind in x (East) direction on defined mesh in (m/s)
    Vygridbackwind
        Hurricane wind velocity with background wind in y (North) direction on defined mesh in (m/s)
    Vgridbackwind
        | Resultant hurricane wind velocity (Vx^2+Vy^2)^0.5 with background wind on defined mesh in (m/s)
        | Gradient wind velocity converted to wind velocity at 10 m above surface by V=VgToVCoeff*Vg
    Rgrid
        | Distance (radius) from hurricane center to each point on the grid
        | Note: Outputs has dimension of [M,N,L] where [M,N] is size of the x-y grid and [L] is number of time steps

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np


        #EXAMPLE 1

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

        #Calculating distance and angle using spherical law of cosines
        Rgrid=(np.arccos(np.sin(np.deg2rad(latCenter))*np.sin(np.deg2rad(ygrid))+np.cos(np.deg2rad(latCenter))*np.cos(np.deg2rad(ygrid))*np.cos(np.deg2rad(xgrid)-np.deg2rad(longCenter))))*6371000 #Radius
        thetagrid=np.arctan2(np.sin(np.deg2rad(xgrid)-np.deg2rad(longCenter))*np.cos(np.deg2rad(ygrid)),np.cos(np.deg2rad(latCenter))*np.sin(np.deg2rad(ygrid))-np.sin(np.deg2rad(latCenter))*np.cos(np.deg2rad(ygrid))*np.cos(np.deg2rad(xgrid)-np.deg2rad(longCenter))) #Azimuth in radian
        thetagrid=-thetagrid+np.pi/2 #Converting azimuth to trigonometric direction
        thetagrid=thetagrid+np.pi/2 #Angle of velocity vector in degree (right angle respect to radius)

        #Calculating hurricane velocity at each radius using SLOSH model
        RVmax=32197 #Radius from hurricane center to a location of maximum hurricane wind
        Vgrid=Vgmax*(2*RVmax*Rgrid)/((RVmax)**2+(Rgrid)**2) #Hurricane wind velocity at radius R
        Vxgrid=Vgrid*np.cos(thetagrid) #Hurricane velocity in x (East) direction
        Vygrid=Vgrid*np.sin(thetagrid) #Hurricane velocity in y (North) direction

        Vxgridbackwind,Vygridbackwind,Vgridbackwind,Rgrid=sm.hurricanebackgroundwind(xgrid,ygrid,Vxgrid,Vygrid,longCenter,latCenter,Vt,VtAzmdir,RVmax,'slosh',423e3,'gc','quiver')


        #EXAMPLE 2

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

        RVmax=32197 #Radius from hurricane center to a location of maximum hurricane wind

        Vxgrid=0 #Hurricane velocity in x (East) direction
        Vygrid=0 #Hurricane velocity in y (North) direction

        Vxgridbackwind,Vygridbackwind,Vgridbackwind,Rgrid=sm.hurricanebackgroundwind(xgrid,ygrid,Vxgrid,Vygrid,longCenter,latCenter,Vt,VtAzmdir,RVmax,'slosh',423e3,'gc','quiver')

    
    References
    ----------

    Data

    * www.nhc.noaa.gov/data/
    * www.nhc.noaa.gov/data/hurdat/hurdat2-format-nencpac.pdf
    * coast.noaa.gov/hurricanes
    * www.aoml.noaa.gov/hrd/data_sub/re_anal.html

    Chavas, D. R., & Emanuel, K. A. (2010). 
    A QuikSCAT climatology of tropical cyclone size. 
    Geophysical Research Letters, 37(18).

    Jelesnianski, C. P., Chen, J., & Shaffer, W. A. (1992). 
    SLOSH: Sea, lake, and overland surges from hurricanes (Vol. 48). 
    US Department of Commerce, National Oceanic and Atmospheric Administration, National Weather Service.

    Lin, N., & Chavas, D. (2012). 
    On hurricane parametric wind and applications in storm surge modeling. 
    Journal of Geophysical Research: Atmospheres, 117(D9).

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
    Vt=type2numpy(Vt)
    VtAzmdir=type2numpy(VtAzmdir)
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
        Vgridbackwind=np.zeros((Mgrid,Ngrid,len(xCenter))) #Pre-assigning array
        Vxgridbackwind=np.zeros((Mgrid,Ngrid,len(xCenter))) #Pre-assigning array
        Vygridbackwind=np.zeros((Mgrid,Ngrid,len(xCenter))) #Pre-assigning array
    else:
        Rgrid=np.zeros((Mgrid,Ngrid,1)) #Pre-assigning array
        thetagrid=np.zeros((Mgrid,Ngrid,1)) #Pre-assigning array
        Vgridbackwind=np.zeros((Mgrid,Ngrid,1)) #Pre-assigning array
        Vxgridbackwind=np.zeros((Mgrid,Ngrid,1)) #Pre-assigning array
        Vygridbackwind=np.zeros((Mgrid,Ngrid,1)) #Pre-assigning array

    
    VtR=np.zeros((Mgrid,Ngrid)) #Pre-assigning array

    if ((len(xCenter)==1) and (np.size(Vxgrid)!=1)):
        Vxgrid=np.reshape(Vxgrid,(Mgrid,Ngrid,1))
        Vygrid=np.reshape(Vygrid,(Mgrid,Ngrid,1))
    elif ((len(xCenter)==1) and (np.size(Vxgrid)==1)):
        Vxgrid=np.reshape(Vxgrid,(1,1,1))
        Vygrid=np.reshape(Vygrid,(1,1,1))

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
    #Adding background wind velocity due to hurricane motion to hurricane wind velocity

    #Converting azimuth (bearing) to trigonometric direction
    VtTrigdir=-VtAzmdir+90
    
    #Add 360 to all numbers to have them all positive
    #Use mod(360) to take care of the ones larger than 360
    VtTrigdir=((VtTrigdir+360)%360) 
    
    #Adding background wind velocity as constant value to all hurricane velocities
    if backwindCalcMethod=='constant':
    
        for i in range(0,len(xCenter),1):
            VtR[:]=Vt[i] #Translational wind velocity at radius R
            VtR[Rgrid[:,:,i]>=Rmax]=0 #e.g. Chavas & Emanuel (2010), Lin & Chavas (2012)
    
            Vxgridbackwind[:,:,i]=Vxgrid[:,:,i]+VtR*np.cos(np.deg2rad(VtTrigdir[i])) #Hurricane velocity in x (East) direction
            Vygridbackwind[:,:,i]=Vygrid[:,:,i]+VtR*np.sin(np.deg2rad(VtTrigdir[i])) #Hurricane velocity in y (North) direction
    

        Vgridbackwind=np.sqrt(Vxgridbackwind**2+Vygridbackwind**2) #Assigning hurricane velocity at each track point
    
    #Adding background wind velocity using SLOSH model by Jelesnianski et al. (1992)
    elif backwindCalcMethod=='slosh':
    
        for i in range(0,len(xCenter),1):
            VtR=Vt[i]*(RVmax[i]*Rgrid[:,:,i])/((RVmax[i])**2+(Rgrid[:,:,i])**2) #Translational wind velocity at radius R
            VtR[VtR<0]=0
            VtR[Rgrid[:,:,i]>=Rmax]=0 #e.g. Chavas & Emanuel (2010), Lin & Chavas (2012)
    
            Vxgridbackwind[:,:,i]=Vxgrid[:,:,i]+VtR*np.cos(np.deg2rad(VtTrigdir[i])) #Hurricane velocity in x (East) direction
            Vygridbackwind[:,:,i]=Vygrid[:,:,i]+VtR*np.sin(np.deg2rad(VtTrigdir[i])) #Hurricane velocity in y (North) direction

    
        Vgridbackwind=np.sqrt(Vxgridbackwind**2+Vygridbackwind**2) #Assigning hurricane velocity at each track point
    
    #Adding background wind velocity using Lin & Chavas (2012)
    elif backwindCalcMethod=='lin':
    
        for i in range(0,len(xCenter),1):
            #Hurricane is in northern hemisphere
            if np.mean(yCenter)>=0:
                #In northern hemisphere, velocity vector rotate counter clockwise with right angle respect to radius
    
                VtR[:]=0.55*Vt[i] #Translational wind velocity at radius R
                VtR[Rgrid[:,:,i]>=Rmax]=0 #e.g. Chavas & Emanuel (2010), Lin & Chavas (2012)
    
                Vxgridbackwind[:,:,i]=Vxgrid[:,:,i]+VtR*np.cos(np.deg2rad(VtTrigdir[i]+20)) #Hurricane velocity in x (East) direction
                Vygridbackwind[:,:,i]=Vygrid[:,:,i]+VtR*np.sin(np.deg2rad(VtTrigdir[i]+20)) #Hurricane velocity in y (North) direction
    
            #Hurricane is in southern hemisphere
            elif np.mean(yCenter)<0:
                #In southern hemisphere, velocity vector rotate clockwise with right angle respect to radius
    
                VtR[:]=0.55*Vt[i] #Translational wind velocity at radius R
                VtR[Rgrid[:,:,i]>=Rmax]=0 #e.g. Chavas & Emanuel (2010), Lin & Chavas (2012)
    
                Vxgridbackwind[:,:,i]=Vxgrid[:,:,i]+VtR*np.cos(np.deg2rad(VtTrigdir[i]-20)) #Hurricane velocity in x (East) direction
                Vygridbackwind[:,:,i]=Vygrid[:,:,i]+VtR*np.sin(np.deg2rad(VtTrigdir[i]-20)) #Hurricane velocity in y (North) direction
    

        Vgridbackwind=np.sqrt(Vxgridbackwind**2+Vygridbackwind**2) #Assigning hurricane velocity at each track point
    

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='imagesc': #2 dimensional plot using imagesc or imshow
        
        plt.imshow(Vgridbackwind[:,:,-1],extent=[np.nanmin(xgrid),np.nanmax(xgrid),np.nanmin(ygrid),np.nanmax(ygrid)],interpolation=None,cmap=plt.get_cmap(),aspect='auto',origin='lower')
    
    
    elif dispout=='pcolor': #2 dimensional plot using pcolor
    
        plt.pcolormesh(xgrid,ygrid,Vgridbackwind[:,:,-1],cmap=plt.get_cmap())    
    

    elif dispout=='contour': #2 dimensional contour plot
    
        plt.contourf(xgrid,ygrid,Vgridbackwind[:,:,-1],cmap=plt.get_cmap())
    
    
    elif dispout=='quiver': #Surface plot
        
        plt.contour(xgrid,ygrid,Vgridbackwind[:,:,-1],cmap=plt.get_cmap())
        plt.colorbar()
        plt.quiver(xgrid,ygrid,Vxgridbackwind[:,:,-1],Vygridbackwind[:,:,-1])
    
    
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
        Vxgridbackwind=np.reshape(Vxgridbackwind,(Mgrid,Ngrid)) #Reshaping array
        Vygridbackwind=np.reshape(Vygridbackwind,(Mgrid,Ngrid)) #Reshaping array
        Vgridbackwind=np.reshape(Vgridbackwind,(Mgrid,Ngrid)) #Reshaping array
        Rgrid=np.reshape(Rgrid,(Mgrid,Ngrid)) #Reshaping array

    elif ((len(xCenter)==1) and (np.ndim(xgrid)==1)):
        #xgrid and ygrid are 1d array
        #Converting output array with size of (Mgrid,1,1) to array with size of (Mgrid,)
        Vxgridbackwind=np.reshape(Vxgridbackwind,(Mgrid)) #Reshaping array
        Vygridbackwind=np.reshape(Vygridbackwind,(Mgrid)) #Reshaping array
        Vgridbackwind=np.reshape(Vgridbackwind,(Mgrid)) #Reshaping array
        Rgrid=np.reshape(Rgrid,(Mgrid)) #Reshaping array
    
    #--------------------------------------------------------------------------
    #Outputs
    return Vxgridbackwind, Vygridbackwind, Vgridbackwind, Rgrid

    #--------------------------------------------------------------------------
