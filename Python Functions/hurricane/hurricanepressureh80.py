def hurricanepressureh80(xgrid, ygrid, xCenter, yCenter, Pc, Vgmax, Rknown, VgRknown, Pn=101325, Rhoa=1.204, \
    distCalcMethod='gc', flattendata='no', savedata='no', dispout='no'):
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

    scientimate.hurricanepressureh80
    ================================

    .. code:: python

        Pgrid, Rgrid, RVmax, thetagrid = scientimate.hurricanepressureh80(xgrid, ygrid, xCenter, yCenter, Pc, Vgmax, Rknown, VgRknown, Pn=101325, Rhoa=1.204, \
            distCalcMethod='gc', flattendata='no', savedata='no', dispout='no')

    Description
    -----------

    Generate hurricane pressure data on given (x,y) points using method from Holland (1980)

    Inputs
    ------

xgrid
                                :x (longitude) of points which outputs are calculated at
                                :xgrid can be a single point or 1d or 2d array 
ygrid
                                :y (latitude) of points which outputs are calculated at
                                :ygrid can be a single point or 1d or 2d array 
xCenter
                                :x (longitude) of hurricane center (track)
yCenter
                                :y (latitude) of hurricane center (track)
Pc
                                :Hurricane central surface pressure in (Pa)
Vgmax
                                :Maximum hurricane 1-min averaged wind velocity at the gradient level in (m/s)
                                :It can be estimated from Vmax=3.44*(1010-Pc*1e-2)^0.644; where Pc is in (Pa)
                                :and, Vmax is 1-min averaged wind at a 10-m elevation in (m/s) (Atkinson & Holliday, 1977)
                                :Hurricane translation velocity (forward velocity) needs to be removed from Vgmax before applying (e.g. Hu et al., 2012)
                                :Hurricane translation velocity (forward velocity) needs to be added after rotational velocity is calculated
Rknown
                                :Radius that hurricane wind velocity is known at that radius in (m)
                                :Rknown should be larger than radius associated to Vmax
VgRknown
                                :Hurricane wind velocity at the gradient level which is known at radius Rknown in (m/s)
Pn=101325
                                :Ambient surface pressure (external pressure) in (Pa)
                                :Standard atmosphere pressure is 101325 (Pa) 
                                :Typical value: Pn=101500 (Pa) for the western North Pacific, Pn= 101000 (Pa) for the North Atlantic
                                :(Batke et al., 2014)
Rhoa=1.204
                                :Air density at the gradient level in (kg/m3)
distCalcMethod='gc'
                                :Distance calculation method 
                                :'cart': Distances are calculated on cartesian coordinate
                                :'gc': Distances are calculated on Great Circle based on Vincenty formula, Vincenty (1975)
                                :Earth radius coonsidered as mean earth radius=6371000 m
flattendata='no'
                                :Define if flat data or not
                                :'no': does not flat the results, outputs are in 3d array
                                :'yes': flat the results, outputs are in 2d array
savedata='no'
                                :Define if save data in a file or not in working folder
                                :'no': does not save, 
                                :'yes': save data as ascii 'dat' file, data are flatten regrdless of flattendata value
dispout='no'
                                :Define to display outputs or not
                                :'imagesc': 2 dimensional plot using imagesc or imshow
                                :'pcolor': 2 dimensional plot using pcolor
                                :'contour': 2 dimensional contour plot, number of contour=ncolor
                                :'no': not display 
                                :Use dispout='no' if calculation mesh is not 2d array
                                :if there is more than one time step, only the last one is plotted
                                :if flattendata='yes' then dispout is set as dispout='no'

    Outputs
    -------

Pgrid                           
                                :Hurricane surface pressure on defined mesh in (Pa)
Rgrid                           
                                :Distance (radius) from hurricane center to each point on the grid
RVmax
                                :Distance (radius) from hurricane center to a location of maximum hurricane wind velocity (m)
thetagrid
                                :Angle from hurricane center to each point on the grid in (Degree)

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

        Pgrid,Rgrid,RVmax,thetagrid=sm.hurricanepressureh80(xgrid,ygrid,longCenter,latCenter,Pc,Vgmax,Rknown,VgRknown,Pn,Rhoa,\
            'gc','no','no','imagesc')


        #EXAMPLE 2

        #Creating calculation mesh
        xgrid,ygrid=np.meshgrid(np.linspace(-98,-68,100),np.linspace(16,44,100))

        #Longitude of Hurricane Katrine best track
        longtrack=[-75.1,-75.7,-76.2,-76.5,-76.9,-77.7,-78.4,-79.0,-79.6,-80.1,-80.3,-81.3,\
            -82.0,-82.6,-83.3,-84.0,-84.7,-85.3,-85.9,-86.7,-87.7,-88.6,-89.2,-89.6,\
            -89.6,-89.6,-89.6,-89.6,-89.1,-88.6,-88.0,-87.0,-85.3,-82.9]

        #Latitude of Hurricane Katrine best track
        lattrack=[23.1,23.4,23.8,24.5,25.4,26.0,26.1,26.2,26.2,26.0,25.9,25.4,\
            25.1,24.9,24.6,24.4,24.4,24.5,24.8,25.2,25.7,26.3,27.2,28.2,\
            29.3,29.5,30.2,31.1,32.6,34.1,35.6,37.0,38.6,40.1]

        #Hurricane Katrina centeral pressure (Pa)
        Pc=[100800,100700,100700,100600,100300,100000,99700,99400,98800,98400,98300,98700,\
            97900,96800,95900,95000,94200,94800,94100,93000,90900,90200,90500,91300,\
            92000,92300,92800,94800,96100,97800,98500,99000,99400,99600]

        #Hurricane Katrina translational velocity (m/s)
        Vt=np.array([0.00000,3.23091,3.13105,3.86928,4.99513,4.82816,3.27813,2.81998,2.77140,2.53041,\
            1.05928,5.30662,3.60661,2.98269,3.61863,3.43691,3.28168,2.85849,3.20404,4.26279,\
            5.31340,5.18467,5.39195,5.46121,5.66270,1.02958,3.60354,4.63312,8.02540,8.01558,\
            8.12721,8.31580,10.75406,12.28350])
            
        #Hurricane Katrina 1-min sustained maximum velocity (m/s)
        Vmax=np.array([15.3,15.3,15.3,17.850,20.4,22.950,25.5,28.050,30.6,35.7,35.7,33.150,\
            38.250,43.350,45.9,48.450,51.0,51.0,51.0,63.750,73.950,76.5,71.4,63.750,\
            56.1,56.1,53.550,40.8,25.5,20.4,15.3,15.3,15.3,12.750])

        Vmax=Vmax-Vt #Removing hurricane translation velocity from Vmax
        Vgmax=Vmax/0.8 #Converting surface velocity to gradient velocity

        #34 kt (17.49 m/s) wind radii maximum extent in northeastern quadrant in (m) for Hurricane Katrina
        RknownRaw=[0,0,0,111120,111120,111120,111120,111120,129640,np.nan,129640,138900,\
            138900,138900,166680,240760,240760,259280,259280,296320,333360,370400,370400,370400,\
            np.nan,370400,np.nan,185200,138900,138900,0,0,0,0]

        #34 kt (17.49 m/s) wind radii maximum extent in northeastern quadrant in (m) for Hurricane Katrina
        Rknown=[0,0,0,111120,111120,111120,111120,111120,129640,129640,129640,138900,\
            138900,138900,166680,240760,240760,259280,259280,296320,333360,370400,370400,370400,\
            370400,370400,277800,185200,138900,138900,0,0,0,0]
        VRknown=np.ones(34)*17.49
        VRknown=VRknown-Vt #Removing hurricane translation velocity from VRknown
        VgRknown=VRknown/0.8 #Converting surface velocity to gradient velocity

        Pn=101325 #Ambient surface pressure (external pressure) in (Pa)
        Rhoa=1.204 #Air density in (kg/m3)

        Pgrid,Rgrid,RVmax,thetagrid=sm.hurricanepressureh80(xgrid,ygrid,longtrack[3:27],lattrack[3:27],Pc[3:27],Vgmax[3:27],Rknown[3:27],VgRknown[3:27],Pn,Rhoa,\
            'gc','no','no','imagesc')


        #EXAMPLE 3

        xgrid=np.linspace(0,10,100) #(Degree)
        ygrid=np.ones(100)*20 #(Degree)
        longCenter=0 #(Degree)
        latCenter=20 #(Degree)
        Pc=90200 #(Pa)
        Vt=5.18467 #(m/s)
        Vmax=76.5 #(m/s)
        Vmax=Vmax-Vt
        Vgmax=Vmax/0.8 #(m/s)
        Rknown=370400 #(m)
        VRknown=17.49 #(m/s)
        VRknown=VRknown-Vt
        VgRknown=VRknown/0.8 #(m/s)
        Pn=101325 #Ambient surface pressure (external pressure) in (Pa)
        Rhoa=1.204 #Air density in (kg/m3)

        Pgrid,Rgrid,RVmax,thetagrid=sm.hurricanepressureh80(xgrid,ygrid,longCenter,latCenter,Pc,Vgmax,Rknown,VgRknown,Pn,Rhoa,\
        'gc','no','no','no')
        plt.plot(Rgrid,Pgrid)

    References
    ----------

    Data
    
    * www.nhc.noaa.gov/data/
    * www.nhc.noaa.gov/data/hurdat/hurdat2-format-nencpac.pdf
    * coast.noaa.gov/hurricanes
    * www.aoml.noaa.gov/hrd/data_sub/re_anal.html

    Atkinson, G. D., & Holliday, C. R. (1977). 
    Tropical cyclone minimum sea level pressure/maximum sustained wind relationship for the western north Pacific. 
    Monthly Weather Review, 105(4), 421-427.

    Batke, S. P., Jocque, M., & Kelly, D. L. (2014). 
    Modelling hurricane exposure and wind speed on a mesoclimate scale: a case study from Cusuco NP, Honduras. 
    PloS one, 9(3), e91306.

    Holland, G. J. (1980). 
    An analytic model of the wind and pressure profiles in hurricanes. 
    Monthly weather review, 108(8), 1212-1218.

    Holland, G. (2008). 
    A revised hurricane pressure–wind model. 
    Monthly Weather Review, 136(9), 3432-3445.

    Holland, G. J., Belanger, J. I., & Fritz, A. (2010). 
    A revised model for radial profiles of hurricane winds. 
    Monthly Weather Review, 138(12), 4393-4401.

    Phadke, A. C., Martino, C. D., Cheung, K. F., & Houston, S. H. (2003). 
    Modeling of tropical cyclone winds and waves for emergency management. 
    Ocean Engineering, 30(4), 553-578.

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
    from scipy import optimize
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
    Pc=type2numpy(Pc)
    Vgmax=type2numpy(Vgmax)
    Rknown=type2numpy(Rknown)
    VgRknown=type2numpy(VgRknown)
    Pn=type2numpy(Pn)

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
        Pgrid=np.zeros((Mgrid,Ngrid,len(xCenter))) #Pre-assigning array
    else:
        Rgrid=np.zeros((Mgrid,Ngrid,1)) #Pre-assigning array
        thetagrid=np.zeros((Mgrid,Ngrid,1)) #Pre-assigning array
        Pgrid=np.zeros((Mgrid,Ngrid,1)) #Pre-assigning array

    
    RVmax=np.zeros(len(xCenter)) #Pre-assigning array

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
    #Calculating hurricane wind and pressure using Holland (1980)

    deltaP=(Pn-Pc) #Pressure drop from a defined external pressure to the cyclone center
    #deltaP=deltaP*0.01 #Convert from (Pa) to (hPa)
    
    B=Vgmax**2*Rhoa*np.exp(1)/(deltaP) #Hurricane shape parameter
    B[B>2.5]=2.5 #B is typically between 1 and 2.5
    B[B<1]=1 #B is typically between 1 and 2.5
    
    #w=7.2921150e-5 #Rotational frequency of the earth (Radian/s)
    #f=2*w*np.sin(np.deg2rad(yCenter)) #Coriolis parameter (Coriolis frequency)
    
    for i in range(0,len(xCenter),1):
    
        #Initial guess for RVmax using Rankine vortex model
        #Rankine vortex model for R>=RVmax: V=Vgmax*(RVmax/R)^x, where 0.4<x<0.6 (e,g. Phadke et al., 2003; Holland et al., 2010)
        #Rankine vortex model for R<RVmax: V=Vgmax*(R/RVmax)^x, where 0.4<x<0.6 (e,g. Phadke et al., 2003; Holland et al., 2010)
        RVmaxini=Rknown[i]*(VgRknown[i]/Vgmax[i])**2 #Initial guess
    
        #Calculating radius of maximum wind from known velocity and its radius
        fun = lambda x : (VgRknown[i]-(B[i]/Rhoa*(x/Rknown[i])**B[i]*deltaP[i]*np.exp(-(x/Rknown[i])**B[i]))**0.5) #Hurricane radius for Vgmax without coriolis parameter
        #fun = lambda x : (VgRknown[i]-(B[i]/Rhoa*(x/Rknown[i])**B[i]*deltaP[i]*np.exp(-(x/Rknown[i])**B[i])+(Rknown[i]*f[i]/2)**2)**0.5-(Rknown[i]*f[i]/2)) #Hurricane radius for Vgmax with coriolis parameter
        RVmax[i]=sp.optimize.fsolve(fun,RVmaxini) #Radius of maximum wind
    
        #Surface pressure at radius R
        Ps=Pc[i]+deltaP[i]*np.exp(-(RVmax[i]/Rgrid[:,:,i])**B[i]) #Surface pressure at radius R from hurricane center
        
        Pgrid[:,:,i]=Ps.copy() #Assigning hurricane pressure at each track point
    

    #--------------------------------------------------------------------------
    #Flating results

    if flattendata=='yes':
    
        #Pgrid is 3d array
        if np.ndim(Pgrid)==3: 
            M,N,L=np.shape(Pgrid)
    
        #Pgrid is 2d array
        elif np.ndim(Pgrid)==2: 
            M,N=np.shape(Pgrid)
            L=1 #L=1 when Pgrid is 2d array
    
        #Pgrid is 1d array
        else: 
            M=(np.shape(Pgrid))[0]
            N=1 #N=1 when Pgrid is 1d array
            L=1 #L=1 when Pgrid is 1d array
        
        Pgrid=np.reshape(Pgrid,(M*L,N))
    
        dispout='no' #Results of flatten array cannot be plotted
    

    #--------------------------------------------------------------------------
    #Saving data

    if savedata=='yes':
    
        if flattendata=='yes':
    
            np.savetxt('P.dat',Pgrid,delimiter=' ')
    
        #Flating results if they were not flatten before
        elif flattendata=='no':
    
            #Pgrid is 3d array
            if np.ndim(Pgrid)==3: 
                M,N,L=np.shape(Pgrid)
        
            #Pgrid is 2d array
            elif np.ndim(Pgrid)==2: 
                M,N=np.shape(Pgrid)
                L=1 #L=1 when xgrid and ygrid are 2d array
        
            #Pgrid is 1d array
            else: 
                M=(np.shape(Pgrid))[0]
                N=1 #N=1 when xgrid and ygrid are 1d array
                L=1 #L=1 when xgrid and ygrid are 1d array

            Pgridrshp=np.reshape(Pgrid,(M*L,N))
     
            np.savetxt('P.dat',Pgridrshp,delimiter=' ')
    

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='imagesc': #2 dimensional plot using imagesc or imshow
        
        plt.imshow(Pgrid[:,:,-1],extent=[np.nanmin(xgrid),np.nanmax(xgrid),np.nanmin(ygrid),np.nanmax(ygrid)],interpolation=None,cmap=plt.get_cmap(),aspect='auto',origin='lower')
    
    elif dispout=='pcolor': #2 dimensional plot using pcolor
    
        plt.pcolormesh(xgrid,ygrid,Pgrid[:,:,-1],cmap=plt.get_cmap())    
    
    elif dispout=='contour': #2 dimensional contour plot
    
        plt.contourf(xgrid,ygrid,Pgrid[:,:,-1],cmap=plt.get_cmap())
    
    
    #Setting plot properties
    if ((dispout=='imagesc') or (dispout=='pcolor') or (dispout=='contour')):
    
        #Setting axis limits
        #plt.xlim(np.nanmin(xgrid),np.nanmax(xgrid))
        #plt.ylim(np.nanmin(ygrid),np.nanmax(ygrid))
    
        #Plotting colorbar
        plt.colorbar()
    
        #Setting label
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Pressure (Pa)')
        

    #--------------------------------------------------------------------------
    #Converting output array with size of (Mgrid,Ngrid,1) to array with size of (Mgrid,Ngrid)

    if ((len(xCenter)==1) and (np.ndim(xgrid)>1)):
        #xgrid and ygrid are 2d array
        #Converting array with size of (Mgrid,Ngrid,1) to array with size of (Mgrid,Ngrid)
        Pgrid=np.reshape(Pgrid,(Mgrid,Ngrid)) #Reshaping array
        Rgrid=np.reshape(Rgrid,(Mgrid,Ngrid)) #Reshaping array
        thetagrid=np.reshape(thetagrid,(Mgrid,Ngrid)) #Reshaping array

    elif ((len(xCenter)==1) and (np.ndim(xgrid)==1)):
        #xgrid and ygrid are 1d array
        #Converting output array with size of (Mgrid,1,1) to array with size of (Mgrid,)
        Pgrid=np.reshape(Pgrid,(Mgrid)) #Reshaping array
        Rgrid=np.reshape(Rgrid,(Mgrid)) #Reshaping array
        thetagrid=np.reshape(thetagrid,(Mgrid)) #Reshaping array
    
    #--------------------------------------------------------------------------
    #Outputs
    return Pgrid, Rgrid, RVmax, thetagrid

    #--------------------------------------------------------------------------
