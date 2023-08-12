def interpxyz2xyz(x, y, z, xPoint, yPoint, interpMethod='nearest', dispout='no'):
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

    scientimate.interpxyz2xyz
    =========================

    .. code:: python

        zPoint = scientimate.interpxyz2xyz(x, y, z, xPoint, yPoint, interpMethod='nearest', dispout='no')

    Description
    -----------

    | Interpolate 2d scattered data on given point(s) by down-sampling the input data 
    | For down-sampling, the first point of the k-nearest neighbors calculated from Euclidean distance is used
    | interpxyz2xyz is more suitable for large dataset
    | griddata is more efficient than interpxyz2xyz for regular size dataset

    Inputs
    ------

    x
        Coordinate of data in x direction, as an 1D array
    y
        Coordinate of data in y direction, as an 1D array
    z
        Value of data at (x,y) as z(x,y), as an 1D array
    xPoint
        Coordinate (in x direction) of point that nearest point to that is desired to be found 
    yPoint
        Coordinate (in y direction) of point that nearest point to that is desired to be found 
    interpMethod='nearest'
        | Interpolation method 
        | 'linear': Use default or 'linear' method to interpolate
        | 'nearest': Use nearest neighbor method to interpolate
    dispout='no'
        | Define to display outputs or not ('yes': display, 'no': not display)
        | '2d': 2 dimensional scatter plot 
        | 'surface': 3 dimensional surface plot 
        | 'no': not display 

    Outputs
    -------

    zPoint
        Value of interpolated data at (xPoint,yPoint) as z(xPoint,yPoint) 

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        x=10*np.random.rand(100)
        y=10*np.random.rand(100)
        z=100*np.random.rand(100)
        xPoint=[2.5,5,7.5]
        yPoint=[3,6,9]
        zPoint=sm.interpxyz2xyz(x,y,z,xPoint,yPoint,'linear','2d')

        x=10*np.random.rand(100)
        y=10*np.random.rand(100)
        z=100*np.random.rand(100)
        xgrid,ygrid=np.meshgrid(np.linspace(np.min(x),np.max(x),100),np.linspace(np.min(y),np.max(y),100))
        zgrid=sm.interpxyz2xyz(x,y,z,xgrid,ygrid,'nearest','no')

        x=10*np.random.rand(100)
        y=10*np.random.rand(100)
        z=y*np.sin(x)-x*np.cos(y)
        xPoint=10*np.random.rand(10,1)
        yPoint=10*np.random.rand(10,1)
        zPoint=sm.interpxyz2xyz(x,y,z,xPoint,yPoint,'nearest','no')

        x=10*np.random.rand(100)
        y=10*np.random.rand(100)
        z=y*np.sin(x)-x*np.cos(y)
        xgrid,ygrid=np.meshgrid(np.linspace(np.min(x),np.max(x),100),np.linspace(np.min(y),np.max(y),100))
        zgrid=sm.interpxyz2xyz(x,y,z,xgrid,ygrid,'nearest','surface')

    References
    ----------


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
    from scipy.spatial import distance
    from scipy import interpolate
    if dispout!='no':
        import matplotlib.pyplot as plt 
        from mpl_toolkits.mplot3d import Axes3D

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
    
    x=type2numpy(x)
    y=type2numpy(y)
    z=type2numpy(z)
    xPoint=type2numpy(xPoint)
    yPoint=type2numpy(yPoint)

    #--------------------------------------------------------------------------
    #Distance calculation method
    distCalcMethod='2d'
    #    | Distance calculation method 
    #    | '2d': use 2d array
    #    | '1d': use 1d array
    #    | 'pdist2': use pdist2 function
    #    | 'vector': use vectorized distance 

    #--------------------------------------------------------------------------
    #Downsampling input data by finding the nearest neighbors

    numofneighbors=1 #Number of nearest neighbors to (xPoint,yPoint) that are desired to be found
    
    #M=len(xPoint(:,1));
    #N=len(xPoint(1,:));
    #xPoint and yPoint are 1d array
    if np.ndim(xPoint)==1: 
        M=(np.shape(xPoint))[0]
        N=1 #N=1 when xPoint and yPoint are 1d array
    #xPoint and yPoint are 2d array
    if np.ndim(xPoint)==2: 
        M,N=np.shape(xPoint)

    #First method, using 2d array
    if distCalcMethod=='2d':
    
        if N==1:
            xnearest=np.ones((M)) #Pre-assigning array
            ynearest=np.ones((M)) #Pre-assigning array
            znearest=np.ones((M)) #Pre-assigning array
            for i in range(0,M,1):
                for j in range(0,N,1):
                    Dist=np.sqrt((x-xPoint[i])**2+(y-yPoint[i])**2) #Calculating distance of (x,y) to (xPoint,yPoint)
                    #DistSort=np.sort(Dist) #Sorting distances
                    #Indx=np.argsort(Dist) #Sorting distances
                    #minDist=np.min(Dist) #Finding minimum distances
                    Indx=np.argmin(Dist) #Finding minimum distances
                    xnearest[i]=x[Indx] #x of the nearest point to (xPoint,yPoint)
                    ynearest[i]=y[Indx] #y of the nearest point to (xPoint,yPoint)
                    znearest[i]=z[Indx] #Value of the nearest point to (xPoint,yPoint)

        elif N>1:
            xnearest=np.ones((M,N)) #Pre-assigning array
            ynearest=np.ones((M,N)) #Pre-assigning array
            znearest=np.ones((M,N)) #Pre-assigning array
            for i in range(0,M,1):
                for j in range(0,N,1):
                    Dist=np.sqrt((x-xPoint[i,j])**2+(y-yPoint[i,j])**2) #Calculating distance of (x,y) to (xPoint,yPoint)
                    #DistSort=np.sort(Dist) #Sorting distances
                    #Indx=np.argsort(Dist) #Sorting distances
                    #minDist=np.min(Dist) #Finding minimum distances
                    Indx=np.argmin(Dist) #Finding minimum distances
                    xnearest[i,j]=x[Indx] #x of the nearest point to (xPoint,yPoint)
                    ynearest[i,j]=y[Indx] #y of the nearest point to (xPoint,yPoint)
                    znearest[i,j]=z[Indx] #Value of the nearest point to (xPoint,yPoint)
    
    #Second method, using 1d array
    elif distCalcMethod=='1d':
    
        #Reshaping grid data to 1D vector
        xPoint1d=np.reshape(xPoint,(np.size(xPoint)))
        yPoint1d=np.reshape(yPoint,(np.size(xPoint)))
    
        xnearest1d=np.ones(len(xPoint1d)) #Pre-assigning array
        ynearest1d=np.ones(len(xPoint1d)) #Pre-assigning array
        znearest1d=np.ones(len(xPoint1d)) #Pre-assigning array
        for i in range(0,len(xPoint1d),1):
            Dist=np.sqrt((x-xPoint1d[i])**2+(y-yPoint1d[i])**2) #Calculating distance of (X,Y) to (XPoint,YPoint)
            #DistSort=np.sort(Dist) #Sorting distances
            #Indx=np.argsort(Dist) #Sorting distances
            #minDist=np.min(Dist) #Finding minimum distances
            Indx=np.argmin(Dist) #Finding minimum distances
            xnearest1d[i]=x[Indx] #x of the nearest point to (xPoint,yPoint)
            ynearest1d[i]=y[Indx] #y of the nearest point to (xPoint,yPoint)
            znearest1d[i]=z[Indx] #Value of the nearest point to (xPoint,yPoint)
    
        if N==1:
            xnearest=xnearest1d.copy()
            ynearest=ynearest1d.copy()
            znearest=znearest1d.copy()
        elif N>1:
            xnearest=np.reshape(xnearest1d,(M,N))
            ynearest=np.reshape(ynearest1d,(M,N))
            znearest=np.reshape(znearest1d,(M,N))
    
    #Third method, using pdist2 function
    elif distCalcMethod=='pdist2':
    
        if N==1:
            xnearest=np.ones((M)) #Pre-assigning array
            ynearest=np.ones((M)) #Pre-assigning array
            znearest=np.ones((M)) #Pre-assigning array
            X=np.column_stack((x,y))
            Y=np.column_stack((xPoint,yPoint))
            Dist=sp.spatial.distance.cdist(X,Y) #Calculating distance of (x,y) to (xPoint,yPoint)
            #Dist=sp.spatial.distance.cdist(np.column_stack((x,y)),np.column_stack((xPoint[i,:],yPoint[i,:]))) #Calculating distance of (x,y) to (xPoint,yPoint)
            #DistSort=np.sort(Dist) #Sorting distances
            #Indx=np.argsort(Dist) #Sorting distances
            #minDist=np.min(Dist,axis=0) #Finding minimum distances
            Indx=np.argmin(Dist,axis=0) #Finding minimum distances
            xnearest=x[Indx] #x of the nearest point to (xPoint,yPoint)
            ynearest=y[Indx] #y of the nearest point to (xPoint,yPoint)
            znearest=z[Indx] #Value of the nearest point to (xPoint,yPoint)
        
        elif N>1:
            xnearest=np.ones((M,N)) #Pre-assigning array
            ynearest=np.ones((M,N)) #Pre-assigning array
            znearest=np.ones((M,N)) #Pre-assigning array
            for j in range(0,N,1):
                X=np.column_stack((x,y))
                Y=np.column_stack((xPoint[:,j],yPoint[:,j]))
                Dist=sp.spatial.distance.cdist(X,Y) #Calculating distance of (x,y) to (xPoint,yPoint)
                #Dist=sp.spatial.distance.cdist(np.column_stack((x,y)),np.column_stack((xPoint[i,:],yPoint[i,:]))) #Calculating distance of (x,y) to (xPoint,yPoint)
                #DistSort=np.sort(Dist) #Sorting distances
                #Indx=np.argsort(Dist) #Sorting distances
                #minDist=np.min(Dist,axis=0) #Finding minimum distances
                Indx=np.argmin(Dist,axis=0) #Finding minimum distances
                xnearest[:,j]=x[Indx] #x of the nearest point to (xPoint,yPoint)
                ynearest[:,j]=y[Indx] #y of the nearest point to (xPoint,yPoint)
                znearest[:,j]=z[Indx] #Value of the nearest point to (xPoint,yPoint)

    #Fourth method, using vectorized distance 
    elif distCalcMethod=='vector':
    
        if N==1:
            xnearest=np.ones(M) #Pre-assigning array
            ynearest=np.ones(M) #Pre-assigning array
            znearest=np.ones(M) #Pre-assigning array
            M1,N1=np.shape(np.column_stack((x,y)))
            M2,N2=np.shape(np.column_stack((xPoint,yPoint)))
            X=np.column_stack((x,y))
            Y=np.column_stack((xPoint,yPoint))
            #M1,N1=np.shape(X)
            #M2,N2=np.shape(Y)
            
            X2=np.sum(X**2,axis=1)
            Y2=np.sum(Y**2,axis=1)
            Dist=np.sqrt(np.ones((M1,1))*np.transpose(Y2)+np.dot(X2[:,None],(np.ones((1,M2))))-2*np.dot(X,np.transpose(Y))) #Calculating distance from (x1,y1) to (xPoint,yPoint)

            #X2=np.sum((np.column_stack((x,y)))**2,axis=1)
            #Y2=np.sum((np.column_stack((xPoint[i,:],yPoint[i,:])))**2,axis=1)
            #Dist=np.sqrt(np.ones((M1,1))*np.transpose(Y2)+np.dot(X2[:,None],(np.ones((1,M2))))-2*np.dot((np.column_stack((x,y))),np.transpose(np.column_stack((xPoint[i,:],yPoint[i,:]))))) #Calculating distance from (x1,y1) to (xPoint,yPoint)
    
            #Dist=sp.spatial.distance.cdist(np.column_stack((x,y)),np.column_stack((xPoint[i,:],yPoint[i,:]))) #Calculating distance of (x,y) to (xPoint,yPoint)
            #DistSort=np.sort(Dist) #Sorting distances
            #Indx=np.argsort(Dist) #Sorting distances
            #minDist=np.min(Dist,axis=0) #Finding minimum distances
            Indx=np.argmin(Dist,axis=0) #Finding minimum distances
            xnearest=x[Indx] #x of the nearest point to (xPoint,yPoint)
            ynearest=y[Indx] #y of the nearest point to (xPoint,yPoint)
            znearest=z[Indx] #Value of the nearest point to (xPoint,yPoint)

        elif N>1:
            xnearest=np.ones((M,N)) #Pre-assigning array
            ynearest=np.ones((M,N)) #Pre-assigning array
            znearest=np.ones((M,N)) #Pre-assigning array
            M1,N1=np.shape(np.column_stack((x,y)))
            M2,N2=np.shape(np.column_stack((xPoint[:,1],yPoint[:,1])))
            for j in range(0,N,1):
                X=np.column_stack((x,y))
                Y=np.column_stack((xPoint[:,j],yPoint[:,j]))
                #M1,N1=np.shape(X)
                #M2,N2=np.shape(Y)
                
                X2=np.sum(X**2,axis=1)
                Y2=np.sum(Y**2,axis=1)
                Dist=np.sqrt(np.ones((M1,1))*np.transpose(Y2)+np.dot(X2[:,None],(np.ones((1,M2))))-2*np.dot(X,np.transpose(Y))) #Calculating distance from (x1,y1) to (xPoint,yPoint)
    
                #X2=np.sum((np.column_stack((x,y)))**2,axis=1)
                #Y2=np.sum((np.column_stack((xPoint[i,:],yPoint[i,:])))**2,axis=1)
                #Dist=np.sqrt(np.ones((M1,1))*np.transpose(Y2)+np.dot(X2[:,None],(np.ones((1,M2))))-2*np.dot((np.column_stack((x,y))),np.transpose(np.column_stack((xPoint[i,:],yPoint[i,:]))))) #Calculating distance from (x1,y1) to (xPoint,yPoint)
        
                #Dist=sp.spatial.distance.cdist(np.column_stack((x,y)),np.column_stack((xPoint[i,:],yPoint[i,:]))) #Calculating distance of (x,y) to (xPoint,yPoint)
                #DistSort=np.sort(Dist) #Sorting distances
                #Indx=np.argsort(Dist) #Sorting distances
                #minDist=np.min(Dist,axis=0) #Finding minimum distances
                Indx=np.argmin(Dist,axis=0) #Finding minimum distances
                xnearest[:,j]=x[Indx] #x of the nearest point to (xPoint,yPoint)
                ynearest[:,j]=y[Indx] #y of the nearest point to (xPoint,yPoint)
                znearest[:,j]=z[Indx] #Value of the nearest point to (xPoint,yPoint)


    #--------------------------------------------------------------------------
    #Interpolating data into (xPoint,yPoint)

    #Interpolating data into grid using default or linear method
    if interpMethod=='linear':
        if N==1:
            if M<=2:
                fun=sp.interpolate.NearestNDInterpolator((xnearest,ynearest),znearest)
                zPoint=fun(xPoint,yPoint)
            elif M>2: #need at least 3 data points for griddata
                zPoint=sp.interpolate.griddata((xnearest,ynearest),znearest,(xPoint,yPoint),method=interpMethod)

        elif N>1:
            #Reshaping xnearest, xnearest and xnearest to 1d array for griddata
            zPoint=sp.interpolate.griddata((np.reshape(xnearest,M*N),np.reshape(ynearest,M*N)),np.reshape(znearest,M*N),(xPoint,yPoint),method=interpMethod)
    
        #Replacing NaN data point resulted from interpolation by nearest point
        if np.sum(np.isnan(zPoint))>0:
            zPoint[np.isnan(zPoint)==1]=znearest[np.isnan(zPoint)==1]
    
    #Interpolating data into grid using nearest neighbor method
    elif interpMethod=='nearest':
        zPoint=znearest.copy() 
    

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='2d':
        
        #Plotting data
        plt.scatter(x,y,label='Input Data')
        plt.scatter(xPoint,yPoint,label='Interpolated Data')
        plt.scatter(xnearest,ynearest,label='KNN')
        
        plt.xlabel('x')
        plt.ylabel('y')
        plt.legend()
        
    elif dispout=='surface':
        
        #Plotting data
        fig=plt.figure()
        ax=fig.gca(projection='3d')
        ax.scatter(x,y,z,color='r',label='Input Data')
        surf=ax.plot_surface(xPoint,yPoint,zPoint,cmap=plt.get_cmap())
        
        fig.colorbar(surf)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        fig.legend()


    #--------------------------------------------------------------------------
    #Outputs
    return zPoint

    #--------------------------------------------------------------------------
