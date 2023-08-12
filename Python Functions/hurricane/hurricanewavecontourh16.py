def hurricanewavecontourh16(xgrid, ygrid, Vgrid, xCenter, yCenter, VtAzmdir=0, G=0.88, fetchCalcMethod='constant', distCalcMethod='gc', dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2017-10-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.hurricanewavecontourh16
    ===================================

    .. code:: python

        Hsgrid, Tpgrid, Hsmax, Tpmax = scientimate.hurricanewavecontourh16(xgrid, ygrid, Vgrid, xCenter, yCenter, VtAzmdir=0, G=0.88, fetchCalcMethod='constant', distCalcMethod='gc', dispout='no')

    Description
    -----------

    | Calculate hurricane wave height and wave period field (contours) on given mesh using 
    | method from Hwang (2016) and Hwang & Walsh (2016) 
    | For method from Young (1988) use hurricanewavecontoury88
    | For method from U.S. Army Corps of Engineers use hurricanewavecontourhcem

    Inputs
    ------

    xgrid
        | x (longitude) of points which outputs are calculated at
        | xgrid can be a single point or 1d or 2d array 
    ygrid
        | y (latitude) of points which outputs are calculated at
        | ygrid can be a single point or 1d or 2d array 
    Vgrid
        Resultant hurricane 1-min averaged wind velocity at 10 m above surface (Vx^2+Vy^2)^0.5 on defined mesh in (m/s)
    xCenter
        x (longitude) of hurricane center (track)
    yCenter
        y (latitude) of hurricane center (track)
    VtAzmdir=0
        | Hurricane center velocity azimuth (bearing) direction in (Degree)
        | azimuth (bearing) direction which is measured clockwise from the north:
        | 0 (degree): toward North, 90 (degree): toward East, 180 (degree): toward South, 270 (degree): toward West 
    G=0.88
        | Wind gust factor to convert 1-min averaged wind to 10-min averaged wind
        | e.g. Young (2017); Liu et al. (2017)
        | G=U(600s)/U(60s), therefore U(600s)=G*U(60s)
        | G=1/1.1; based on Powell and Houston (1996)
        | G=1/1.08; based on Harper (2013)
        | G=0.88; based on World Meteorological Organization (2015)
        | G=1/1.08 to 1/1.16; based on Liu et al. (2017)
    fetchCalcMethod='constant'
        | Effective wind fetch calculation method 
        | 'constant': Use constant coefficients needed for effective wind fetch calculation
        | 'interp': Use interpolateed coefficients needed for effective wind fetch calculation
        | Earth radius coonsidered as mean earth radius=6371000 m
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
        | 'no': not display 
        | Use dispout='no' if calculation mesh is not 2d array

    Outputs
    -------

    Hsgrid
        Hurricane significant wave height on grid mesh in (m)
    Tpgrid
        Hurricane peak wave period on grid mesh in (s)
    Hsmax
        Hurricane maximum significant wave height in (m) 
    Tpmax
        | Hurricane maximum peak wave period in (s) 
        | Note: Maximum values of wave height and wave period should be limited to fully developed values

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np
        import matplotlib.pyplot as plt


        #EXAMPLE 1

        #Creating calculation mesh
        xgrid,ygrid=np.meshgrid(np.linspace(-98,-68,100),np.linspace(16,44,100))

        #Longitude of Hurricane Katrine center at max velocity
        longCenter=-88.6

        #Latitude of Hurricane Katrine center at max velocity
        latCenter=26.3

        #Hurricane Katrina centeral pressure (Pa) at max velocity
        Pc=90200

        #Hurricane Katrina translational velocity (m/s) at max velocity
        Vt=5.18467

        #Hurricane Katrina velocity azimuth (bearing) in (Degree) at max velocity
        VtAzmdir=306.76219

        #Hurricane Katrina 1-min sustained maximum velocity (m/s) at max velocity
        Vmax=76.5
        Vmax=Vmax-Vt #Removing hurricane translation velocity from Vmax
        Vgmax=Vmax/0.8 #Converting surface velocity to gradient velocity

        #34 kt (17.49 m/s) wind radii maximum extent in northeastern quadrant in (m) for Hurricane Katrina at max velocity
        Rknown=370400
        VRknown=17.49
        VRknown=VRknown-Vt #Removing hurricane translation velocity from VRknown
        VgRknown=VRknown/0.8 #Converting surface velocity to gradient velocity

        Pn=101325 #Ambient surface pressure (external pressure) in (Pa)
        Rhoa=1.204 #Air density in (kg/m3)

        #Calculating distance (radius) from hurricane center to each point
        Rgrid=(np.arccos(np.sin(np.deg2rad(latCenter))*np.sin(np.deg2rad(ygrid))+np.cos(np.deg2rad(latCenter))*np.cos(np.deg2rad(ygrid))*np.cos(np.deg2rad(xgrid)-np.deg2rad(longCenter))))*6371000

        #Generting wind velocity for Hurricane Katrine at max velocity using SLOSH model
        Vggrid=Vgmax*(2*32197*Rgrid)/(32197**2+Rgrid**2) #Gradient wind velocity
        Vggrid[Rgrid>=423e3]=0 
        Vgrid=Vggrid*0.8 #Wind velocity at 10m height

        Hsgrid,Tpgrid,Hsmax,Tpmax=sm.hurricanewavecontourh16(xgrid,ygrid,Vgrid,longCenter,latCenter,VtAzmdir,0.88,'constant','gc','contour')


        #EXAMPLE 2

        xgrid=np.linspace(0,10,100) #(Degree)
        ygrid=np.ones(100)*20 #(Degree)
        longCenter=0 #(Degree)
        latCenter=20 #(Degree)
        Pc=90200 #(Pa)
        Vt=5.18467 #(m/s)
        VtAzmdir=306.76219 #(Degree) 
        Vmax=76.5 #(m/s)
        Vmax=Vmax-Vt
        Vgmax=Vmax/0.8 #(m/s)
        Rknown=370400 #(m)
        VRknown=17.49 #(m/s)
        VRknown=VRknown-Vt
        VgRknown=VRknown/0.8 #(m/s)
        Pn=101325 #Ambient surface pressure (external pressure) in (Pa)
        Rhoa=1.204 #Air density in (kg/m3)
        Rgrid=(np.arccos(np.sin(np.deg2rad(latCenter))*np.sin(np.deg2rad(ygrid))+np.cos(np.deg2rad(latCenter))*np.cos(np.deg2rad(ygrid))*np.cos(np.deg2rad(xgrid)-np.deg2rad(longCenter))))*6371000
        Vggrid=Vgmax*(2*32197*Rgrid)/(32197**2+Rgrid**2) #Gradient wind velocity
        Vggrid[Rgrid>=423e3]=0 
        Vgrid=Vggrid*0.8 #Wind velocity at 10m height

        Hsgrid,Tpgrid,Hsmax,Tpmax=sm.hurricanewavecontourh16(xgrid,ygrid,Vgrid,longCenter,latCenter,VtAzmdir,0.88,'constant','gc','no')
        plt.plot(Rgrid,Hsgrid)

    References
    ----------

    Data

    * www.nhc.noaa.gov/data/
    * www.nhc.noaa.gov/data/hurdat/hurdat2-format-nencpac.pdf
    * coast.noaa.gov/hurricanes
    * www.aoml.noaa.gov/hrd/data_sub/re_anal.html

    Harper, B.A. (2013)
    Best practice in tropical cyclone wind hazard modelling: In search of data and emptying the skeleton cupboard. 
    In Proceedings of the 16th Australasian Wind Engineering Society Workshop, Brisbane, Qld, Australia, 18–19 July 2013

    Hwang, P. A. (2016). 
    Fetch-and duration-limited nature of surface wave growth inside tropical cyclones: With applications to air–sea exchange and remote sensing. 
    Journal of Physical Oceanography, 46(1), 41-56.

    Hwang, P. A., & Walsh, E. J. (2016). 
    Azimuthal and radial variation of wind-generated surface waves inside tropical cyclones. 
    Journal of Physical Oceanography, 46(9), 2605-2621.

    Liu, Q., Babanin, A., Fan, Y., Zieger, S., Guan, C., & Moon, I. J. (2017). 
    Numerical simulations of ocean surface waves under hurricane conditions: Assessment of existing model performance. 
    Ocean Modelling, 118, 73-93.

    Powell, M. D., & Houston, S. H. (1996). 
    Hurricane Andrew's landfall in South Florida. Part II: Surface wind fields and potential real-time applications. 
    Weather and Forecasting, 11(3), 329-349.

    World Meteorological Organization. Tropical Cyclone Programme, & Holland, G. J. (2015). 
    Global guide to tropical cyclone forecasting. 
    Secretariat of the World Meteorological Organization.

    Young, I.R. (2017)
    A Review of Parametric Descriptions of Tropical Cyclone Wind-Wave Generation.
    Atmosphere 2017, 8, 194.

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
    
    xgrid=type2numpy(xgrid)
    ygrid=type2numpy(ygrid)
    Vgrid=type2numpy(Vgrid)
    xCenter=type2numpy(xCenter)
    yCenter=type2numpy(yCenter)
    VtAzmdir=type2numpy(VtAzmdir)

    #--------------------------------------------------------------------------
    #Pre-assigning array

    #xgrid and ygrid are 2d array
    if np.ndim(xgrid)==2: 
        Mgrid,Ngrid=np.shape(xgrid)

    #xgrid and ygrid are 1d array
    else: 
        Mgrid=(np.shape(xgrid))[0]
        Ngrid=1 #Ngrid=1 when xgrid and ygrid are 1d array

    Rgrid=np.zeros((Mgrid,Ngrid)) #Pre-assigning array
    thetagrid=np.zeros((Mgrid,Ngrid)) #Pre-assigning array
    Ietax=np.zeros((Mgrid,Ngrid)) #Pre-assigning array
    Setax=np.zeros((Mgrid,Ngrid)) #Pre-assigning array
    Iomegax=np.zeros((Mgrid,Ngrid)) #Pre-assigning array
    Somegax=np.zeros((Mgrid,Ngrid)) #Pre-assigning array
    Beta=np.zeros((Mgrid,Ngrid)) #Pre-assigning array

    
    #--------------------------------------------------------------------------
    #Checking input array size

    #xgrid and ygrid are 1d array
    VgridOriginalDim=np.ndim(Vgrid) #Original dimension of the array
    if ((np.ndim(xgrid)==1) and (np.ndim(Vgrid)==1)): 
        #Converting array with size of (Mgrid,) to array with size of (Mgrid,1)
        Vgrid=np.reshape(Vgrid,(Mgrid,Ngrid))

    
    #--------------------------------------------------------------------------
    #Calculating distance (radius) from hurricane center to each point

    #Calculating distance using cartesian formula
    if distCalcMethod=='cart':
    
        distxy=np.sqrt((xgrid-xCenter)**2+(ygrid-yCenter)**2) #Calculating distance from (x,y) to (x(1),y(1))
    
        #Calculating angle of the line between start and end points
    
        thetarad=np.arctan2(ygrid-yCenter,xgrid-xCenter) #Angle in radian
        theta=np.rad2deg(thetarad) #Angle in degree
    
        #Add 360 to all numbers to have them all positive
        #Use mod(360) to take care of the ones larger than 360
        theta=((theta+360)%360) 
    
        #xgrid and ygrid are 2d array
        if np.ndim(xgrid)==2:
            Rgrid=distxy.copy() #Distance from hurricane center to each point in (m)
            thetagrid=theta.copy() #Angle from hurricane center to each point in (degree)

        #xgrid and ygrid are 1d array
        else:
            #Using np.newaxis (or None) to convert shape (Mgrid) into shape (Mgrid,1) when xgrid and ygrid are 1d array
            Rgrid=distxy[:,np.newaxis].copy() #Distance from hurricane center to each point in (m)
            thetagrid=theta[:,np.newaxis].copy() #Angle from hurricane center to each point in (degree)

    #Calculating distance using Vincenty formula
    elif distCalcMethod=='gc':
    
        #Converting to radian
        lat1rad=np.deg2rad(yCenter)
        lon1rad=np.deg2rad(xCenter)
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
            Rgrid=distxy.copy() #Distance from hurricane center to each point in (m)
            thetagrid=azimuth.copy() #Angle from hurricane center to each point in (degree)

        #xgrid and ygrid are 1d array
        else:
            #Using np.newaxis (or None) to convert shape (Mgrid) into shape (Mgrid,1) when xgrid and ygrid are 1d array
            Rgrid=distxy[:,np.newaxis].copy() #Distance from hurricane center to each point in (m)
            thetagrid=azimuth[:,np.newaxis].copy() #Angle from hurricane center to each point in (degree)

    #--------------------------------------------------------------------------
    #Converting trigonometric direction to azimuth (bearing)

    if distCalcMethod=='cart':
    
        #Converting trigonometric direction to azimuth (bearing)
        thetagrid=90-thetagrid
    
        #Add 360 to all numbers to have them all positive
        #Use mod(360) to take care of the ones larger than 360
        thetagrid=((thetagrid+360)%360) 
    

    #--------------------------------------------------------------------------
    #Calculating effective wind fetch

    #Calculating effective wind fetch using costant coefficients
    if fetchCalcMethod=='constant':
    
        #Hurricane is in northern hemisphere
        if np.mean(yCenter)>=0:
            #In northern hemisphere, velocity vector rotate counter clockwise with right angle respect to radius
    
            #Azimuth counterclockwise angle of grid points respect to the hurricane direction
            Beta=VtAzmdir-thetagrid
    
        #Hurricane is in southern hemisphere
        elif np.mean(yCenter)<0:
            #In southern hemisphere, velocity vector rotate clockwise with right angle respect to radius
    
            #Azimuth counterclockwise angle of grid points respect to the hurricane direction
            Beta=-(VtAzmdir-thetagrid)
    
    
        #Add 360 to all numbers to have them all positive
        #Use mod(360) to take care of the ones larger than 360
        Beta=((Beta+360)%360) 
    
        #Calculating effective wind fetch
        RgridKm=Rgrid/1000 #Convert m to km
        xetaxKm=np.zeros((Mgrid,Ngrid)) #Pre-assigning array
        xomegaxKm=np.zeros((Mgrid,Ngrid)) #Pre-assigning array
    
        xetaxKm[((Beta>225) & (Beta<=360))]=-0.26*RgridKm[((Beta>225) & (Beta<=360))]+259.79 #Right sector
        xetaxKm[((Beta>=0) & (Beta<=135))]=1.25*RgridKm[((Beta>=0) & (Beta<=135))]+58.25 #Left sector
        xetaxKm[((Beta>135) & (Beta<=225))]=0.71*RgridKm[((Beta>135) & (Beta<=225))]+30.02 #Back sector
    
        xomegaxKm[((Beta>225) & (Beta<=360))]=0.21*RgridKm[((Beta>225) & (Beta<=360))]+170.00 #Right sector
        xomegaxKm[((Beta>=0) & (Beta<=135))]=2.25*RgridKm[((Beta>=0) & (Beta<=135))]+24.85 #Left sector
        xomegaxKm[((Beta>135) & (Beta<=225))]=0.50*RgridKm[((Beta>135) & (Beta<=225))]+14.16 #Back sector
    
        #Convert km to m
        xetax=xetaxKm*1000 #Convert km to m
        xomegax=xomegaxKm*1000 #Convert km to m
    
    #Calculating effective wind fetch using interpolated coefficients
    elif fetchCalcMethod=='interp':
    
        #Values of Phi, A_etax, a_etax, A_omegax and a_omegax as a function of Phi
        PhiAzmdir=np.array([-13,7,22,67,115,157,200,242,270,292,330,347,367]) #Azimuth angle referenced to the hurricane heading
        IetaxPhi=np.array([100.00,94.72,108.79,77.47,-33.58,37.75,43.45,134.73,107.65,149.25,109.42,100.00,94.72]) #Intercept, I_Eta_x(Phi) or A_Eta_x(Phi)
        SetaxPhi=np.array([0.50,0.14,-0.01,0.93,1.78,0.40,0.70,0.69,0.46,0.30,0.81,0.50,0.14]) #Slope, S_Eta_x(Phi) or a_Eta_x(Phi)
        IomegaxPhi=np.array([100.00,-104.50,79.43,181.65,-9.26,27.46,40.29,114.20,160.67,135.62,35.62,100.00,-104.50]) #Intercept, I_Omega_x(Phi), A_Omega_x(Phi)
        SomegaxPhi=np.array([0.50,3.60,1.12,0.71,2.60,0.29,0.37,0.37,-0.07,0.45,1.59,0.50,3.60]) #Slope, S_Omega_x(Phi) or a_Omega_x(Phi)
    
        #Hurricane is in northern hemisphere
        if np.mean(yCenter)>=0:
            #In northern hemisphere, velocity vector rotate counter clockwise with right angle respect to radius
    
            #Azimuth counterclockwise angle of grid points respect to the hurricane direction
            Beta=VtAzmdir-thetagrid
    
        #Hurricane is in southern hemisphere
        elif np.mean(yCenter)<0:
            #In southern hemisphere, velocity vector rotate clockwise with right angle respect to radius
    
            #Azimuth counterclockwise angle of grid points respect to the hurricane direction
            Beta=-(VtAzmdir-thetagrid)
    
    
        #Add 360 to all numbers to have them all positive
        #Use mod(360) to take care of the ones larger than 360
        Beta=((Beta+360)%360) 
    
#        #Calculating coefficient for given grid azimuth
#        for i in range(0,Mgrid,1):
#            for j in range(0,Ngrid,1):
#                Ietax[i,j]=np.interp(Beta[i,j],PhiAzmdir,IetaxPhi) #Intercept, I_Eta_x or A_Eta_x
#                Setax[i,j]=np.interp(Beta[i,j],PhiAzmdir,SetaxPhi) #Slope, S_Eta_x or a_Eta_x
#                Iomegax[i,j]=np.interp(Beta[i,j],PhiAzmdir,IomegaxPhi) #Intercept, I_Omega_x, A_Omega_x
#                Somegax[i,j]=np.interp(Beta[i,j],PhiAzmdir,SomegaxPhi) #Slope, S_Omega_x or a_Omega_x
        #Calculating coefficient for given grid azimuth
        for j in range(0,Ngrid,1):
            Ietax[:,j]=np.interp(Beta[:,j],PhiAzmdir,IetaxPhi) #Intercept, I_Eta_x or A_Eta_x
            Setax[:,j]=np.interp(Beta[:,j],PhiAzmdir,SetaxPhi) #Slope, S_Eta_x or a_Eta_x
            Iomegax[:,j]=np.interp(Beta[:,j],PhiAzmdir,IomegaxPhi) #Intercept, I_Omega_x, A_Omega_x
            Somegax[:,j]=np.interp(Beta[:,j],PhiAzmdir,SomegaxPhi) #Slope, S_Omega_x or a_Omega_x
    

        Ietax[np.isnan(Ietax)==1]=0
        Setax[np.isnan(Setax)==1]=0
        Iomegax[np.isnan(Iomegax)==1]=0
        Somegax[np.isnan(Somegax)==1]=0
    
        #Calculating effective wind fetch
        RgridKm=Rgrid/1000 #Convert m to km
        xetaxKm=Setax*RgridKm+Ietax #Effective wind fetch for wave height, x_Etax
        xomegaxKm=Somegax*RgridKm+Iomegax #Effective wind fetch for wave period, x_Omegax
    
        #Convert km to m
        xetax=xetaxKm*1000 #Convert km to m
        xomegax=xomegaxKm*1000 #Convert km to m
    

    #--------------------------------------------------------------------------
    #Calculating hurricane wave contour

    #Converting 1-min sustained wind to 10-min averaged wind using gust factor
    #e.g. World Meteorological Organization (2015)
    U10grid=Vgrid*G
    
    Hsgrid=8.10e-4*U10grid**1.19*xetax**0.405
    #Numpy does not allow fractional powers of negative numbers, 
    #Therefore values for xetax<0 is recalculated as follow:
    Hsgrid[xetax<0]=8.10e-4*(U10grid[xetax<0])**1.19*(-1*(np.abs(xetax[xetax<0]))**0.405)

    #Tpgrid=9.28e-2*U10**0.526*xomegax**0.237
    Tpgrid=9.28e-2*U10grid**0.526*xomegax**0.237
    #Numpy does not allow fractional powers of negative numbers, 
    #Therefore values for xetax<0 is recalculated as follow:
    Tpgrid[xomegax<0]=9.28e-2*(U10grid[xomegax<0])**0.526*(-1*(np.abs(xomegax[xomegax<0]))**0.237)
    
    Hsmax=np.max(Hsgrid) #Maximum significant wave height
    Tpmax=np.max(Tpgrid) #Maximum peak wave period

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='imagesc': #2 dimensional plot using imagesc or imshow
        
        plt.subplot(1,2,1)
        plt.imshow(Hsgrid,extent=[np.nanmin(xgrid),np.nanmax(xgrid),np.nanmin(ygrid),np.nanmax(ygrid)],interpolation=None,cmap=plt.get_cmap(),aspect='auto',origin='lower')
        plt.plot([xCenter,xCenter],[np.min(ygrid),np.max(ygrid)],color=[1,1,1])
        plt.plot([np.min(xgrid),np.max(xgrid)],[yCenter,yCenter],color=[1,1,1])
    
        plt.subplot(1,2,2)
        plt.imshow(Tpgrid,extent=[np.nanmin(xgrid),np.nanmax(xgrid),np.nanmin(ygrid),np.nanmax(ygrid)],interpolation=None,cmap=plt.get_cmap(),aspect='auto',origin='lower')
        plt.plot([xCenter,xCenter],[np.min(ygrid),np.max(ygrid)],color=[1,1,1])
        plt.plot([np.min(xgrid),np.max(xgrid)],[yCenter,yCenter],color=[1,1,1])
    
    elif dispout=='pcolor': #2 dimensional plot using pcolor
    
        plt.subplot(1,2,1)
        plt.pcolormesh(xgrid,ygrid,Hsgrid,cmap=plt.get_cmap())    
        plt.plot([xCenter,xCenter],[np.min(ygrid),np.max(ygrid)],color=[1,1,1])
        plt.plot([np.min(xgrid),np.max(xgrid)],[yCenter,yCenter],color=[1,1,1])
    
        plt.subplot(1,2,2)
        plt.pcolormesh(xgrid,ygrid,Tpgrid,cmap=plt.get_cmap())    
        plt.plot([xCenter,xCenter],[np.min(ygrid),np.max(ygrid)],color=[1,1,1])
        plt.plot([np.min(xgrid),np.max(xgrid)],[yCenter,yCenter],color=[1,1,1])
    
    elif dispout=='contour': #2 dimensional contour plot
    
        plt.subplot(1,2,1)
        plt.contourf(xgrid,ygrid,Hsgrid,cmap=plt.get_cmap())
        plt.plot([xCenter,xCenter],[np.min(ygrid),np.max(ygrid)],color=[1,1,1])
        plt.plot([np.min(xgrid),np.max(xgrid)],[yCenter,yCenter],color=[1,1,1])
    
        plt.subplot(1,2,2)
        plt.contourf(xgrid,ygrid,Tpgrid,cmap=plt.get_cmap())
        plt.plot([xCenter,xCenter],[np.min(ygrid),np.max(ygrid)],color=[1,1,1])
        plt.plot([np.min(xgrid),np.max(xgrid)],[yCenter,yCenter],color=[1,1,1])
    
    
    #Setting plot properties
    if ((dispout=='imagesc') or (dispout=='pcolor') or (dispout=='contour')):
    
        plt.subplot(1,2,1)
    
        #Setting axis limits
        #plt.xlim([np.min(xgrid),np.max(xgrid)])
        #plt.ylim([np.min(ygrid),np.max(ygrid)])
    
        #Plotting colorbar
        plt.colorbar()
    
        #Setting label
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Significant Wave Height (m)')
        
        plt.subplot(1,2,2)
    
        #Setting axis limits
        #plt.xlim([np.min(xgrid),np.max(xgrid)])
        #plt.ylim([np.min(ygrid),np.max(ygrid)])
    
        #Plotting colorbar
        plt.colorbar()
    
        #Setting label
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Peak Wave Period (s)')
        

    #--------------------------------------------------------------------------
    #Converting output array with size of (Mgrid,1) to array with size of (Mgrid,)

    #xgrid and ygrid are 1d array
    if ((np.ndim(xgrid)==1) and (VgridOriginalDim==1)): 
        Hsgrid=np.reshape(Hsgrid,(Mgrid)) #Reshaping array
        Tpgrid=np.reshape(Tpgrid,(Mgrid)) #Reshaping array

    
    #--------------------------------------------------------------------------
    #Outputs
    return Hsgrid, Tpgrid, Hsmax, Tpmax

    #--------------------------------------------------------------------------
