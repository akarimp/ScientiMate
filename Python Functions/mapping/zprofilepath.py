def zprofilepath(xgrid, ygrid, zgrid, x, y, distCalcMethod='cart', CalcMethod='nearest', dispout='no'):
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

    scientimate.zprofilepath
    ========================

    .. code:: python

        z, zmean, distxy = scientimate.zprofilepath(xgrid, ygrid, zgrid, x, y, distCalcMethod='cart', CalcMethod='nearest', dispout='no')

    Description
    -----------

    Calculate z (elevation, ...) profile along a path over a given 2d x-y domain (map, image, ...)

    Inputs
    ------

    xgrid
        x (longitude, pixel, ...) data as a [M*N] array
    ygrid
        y (latitude, pixel, ...) data as a [M*N] array
    zgrid
        z (elevation, ...) data as a [M*N] array
    x
        x (longitude, pixel, ...) of the points along the line as (x,y)
    y
        | y (latitude, pixel, ...) of the points along the line as (x,y)
        | If input data are latitude and longitude, they should be in Degree
    distCalcMethod='cart'
        | Distance calculation method 
        | 'cart': Distances are calculated on cartesian coordinate
        | 'gc': Distances are calculated on Great Circle based on Vincenty formula, Vincenty (1975)
        | Earth radius coonsidered as mean earth radius=6371000 m
    CalcMethod='nearest'
        | Interpolation method 
        | 'linear': Use default or 'linear' method to interpolate
        | 'nearest': Use nearest neighbor method to interpolate
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    z
        z (elevation, ...) data along a path at given points (x,y)
    zmean
        Weighted mean of z (elevation) along a line calculated az zmean=1/(total_distance)*(sum(z(i)*dx(i)))
    distxy
        | Distance at each points of (x,y) from the first point of (x(1),y(1))
        | If input data are latitude and longitude in Degree, distxy is in m

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        xgrid,ygrid=np.meshgrid(np.linspace(-5,5,100),np.linspace(-5,5,100))
        zgrid=np.sin(xgrid**2+ygrid**2)/(xgrid**2+ygrid**2)
        x1=-3
        y1=-3
        x2=3
        y2=3
        x=np.linspace(x1,x2,1000)
        y=np.linspace(y1,y2,1000)
        z,zmean,distxy=sm.zprofilepath(xgrid,ygrid,zgrid,x,y,'cart','nearest','yes')

    References
    ----------

    Vincenty, T. (1975). 
    Direct and inverse solutions of geodesics on the ellipsoid with application of nested equations. 
    Survey review, 23(176), 88-93.

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
    x=type2numpy(x)
    y=type2numpy(y)

    #--------------------------------------------------------------------------
    #Calculating z (elevation) profile

    #Interpolating the z (elevation) values for given (x,y) from input data
    if CalcMethod=='linear':
        #fun=sp.interpolate.interp2d(xgrid[0,:],ygrid[:,0],zgrid)
        #z=np.zeros(len(x))
        #for i in range(0,len(x),1):
        #    z[i]=fun(x[i],y[i])
        ##z2d=fun(x,y) #interp2d returns 2d array
        ##z=np.diag(z2d) #Extract a diagonal
        fun=sp.interpolate.LinearNDInterpolator((np.reshape(xgrid,np.size(xgrid)),np.reshape(ygrid,np.size(ygrid))),np.reshape(zgrid,np.size(zgrid)))
        z=fun(x,y)
    
        #Replacing NaN data point resulted from default method with ones from nearest method
        if np.sum(np.isnan(z))>0:
            fun=sp.interpolate.NearestNDInterpolator((np.reshape(xgrid,np.size(xgrid)),np.reshape(ygrid,np.size(ygrid))),np.reshape(zgrid,np.size(zgrid)))
            znearest=fun(x,y)
            z[np.isnan(z)==True]=znearest[np.isnan(z)==True] 

            #M,N=np.shape(xgrid)
            #xmap=np.interp(x,np.array([np.nanmin(xgrid),np.nanmax(xgrid)]),np.array([0,N-1]))
            #ymap=np.interp(y,np.array([np.nanmin(ygrid),np.nanmax(ygrid)]),np.array([0,M-1]))
            #xmap[xmap<0]=0
            #xmap[xmap>N-1]=N-1
            #ymap[ymap<0]=0
            #ymap[ymap>M-1]=M-1
            #znearest=zgrid[np.int64(ymap),np.int64(xmap)]
            #z[np.isnan(z)==True]=znearest[np.isnan(z)==True] 
    
    #Interpolating data into grid using nearest neighbor method
    elif CalcMethod=='nearest':

        #First method
        fun=sp.interpolate.NearestNDInterpolator((np.reshape(xgrid,np.size(xgrid)),np.reshape(ygrid,np.size(ygrid))),np.reshape(zgrid,np.size(zgrid)))
        z=fun(x,y)

        #Second method
        #Mapping values of x and y to associated index
        #M,N=np.shape(xgrid)
        #xmap=np.interp(x,np.array([np.nanmin(xgrid),np.nanmax(xgrid)]),np.array([0,N-1]))
        #ymap=np.interp(y,np.array([np.nanmin(ygrid),np.nanmax(ygrid)]),np.array([0,M-1]))
        #xmap[xmap<0]=0
        #xmap[xmap>N-1]=N-1
        #ymap[ymap<0]=0
        #ymap[ymap>M-1]=M-1
        #z=zgrid[np.int64(ymap),np.int64(xmap)]
    
    #Calculating distance using cartesian formula
    if distCalcMethod=='cart':
    
        distxy=np.sqrt((x-x[0])**2+(y-y[0])**2) #Calculating distance from (x,y) to (x(1),y(1))
    
    #Calculating distance using Vincenty formula
    elif distCalcMethod=='gc':
    
        #Converting to radian
        lat1rad=np.deg2rad(y[0])
        lon1rad=np.deg2rad(x[0])
        lat2rad=np.deg2rad(y)
        lon2rad=np.deg2rad(x)
    
        deltalatrad21=lat2rad-lat1rad
        deltalonrad21=lon2rad-lon1rad
    
        R=6371000 #Earth radius in (m), mean earth radius=6371000 m
        deltasigma=np.arctan2(np.sqrt((np.cos(lat2rad)*np.sin(deltalonrad21))**2+(np.cos(lat1rad)*np.sin(lat2rad)-np.sin(lat1rad)*np.cos(lat2rad)*np.cos(deltalonrad21))**2),np.sin(lat1rad)*np.sin(lat2rad)+np.cos(lat1rad)*np.cos(lat2rad)*np.cos(deltalonrad21)) #Central angle
        arclen=R*deltasigma #Total distance of the line 
        distxy=arclen.copy()
    
    #--------------------------------------------------------------------------
    #Calculating mean z (elevation)

    #Calculating zdx=z(i)*dx(i)
    zdx=(0.5*(z[0:-1]+z[1:]))*np.diff(distxy)
    
    if distxy[-1]!=0:
        #Calculating mean z (elevation) as zmean=1/(total_distance)*(sum(z(i)*dx(i)))
        zmean=np.sum(zdx)/distxy[-1]
    else:
        zmean=z[-1]

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
        
        plt.subplot(2,1,1)
        #Plotting domain
        plt.imshow(zgrid,extent=[np.min(xgrid),np.max(xgrid),np.min(ygrid),np.max(ygrid)],cmap=plt.cm.get_cmap(),aspect='auto',origin='lower')
        plt.colorbar()
    
        #Plotting line 
        plt.plot(x,y)
    
        plt.xlabel('x')
        plt.ylabel('y')
        
        plt.subplot(2,1,2)
        #Plotting depth along a line from point 2 to point1
        plt.plot(distxy,z)
        plt.xlabel('distance')
        plt.ylabel('z')
        
    #--------------------------------------------------------------------------
    #Outputs
    return z, zmean, distxy

    #--------------------------------------------------------------------------
