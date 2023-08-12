def findknn(x, y, xpoint, ypoint, numofneighbors=1, distCalcMethod='1d', dispout='no'):
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

    scientimate.findknn
    ===================

    .. code:: python

        indxknn, distknn = scientimate.findknn(x, y, xpoint, ypoint, numofneighbors=1, distCalcMethod='1d', dispout='no')

    Description
    -----------

    Find k-nearest neighbors using Euclidean distance

    Inputs
    ------

    x
        Coordinate of data in x direction
    y
        Coordinate of data in y direction
    xpoint
        Coordinate (in x direction) of point that nearest point to that is disired to be found 
    ypoint
        Coordinate (in y direction) of point that nearest point to that is disired to be found 
    numofneighbors=1
        Number of nearest neighbors to (xpoint,ypoint) that are desired to be found
    distCalcMethod='1d'
        | Distance calculation method 
        | '1d': use 1d array
        | 'pdist2': Use 2d distance function
        | 'vector': Use vectorized distance 
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    indxknn
        | Index of nearest neighbors points
        | returns M*N array where M=length(xpoint) and N=numofneighbors
        | nth row associated with nth point in (xpoint,ypoint)
    distknn
        | Distance of nearest neighbors points
        | returns M*N array where M=length(xpoint) and N=numofneighbors
        | mth row associated with mth point in (xpoint,ypoint)
        | nth column associated with nth nearest neighbors

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        x=10*np.random.rand(100)
        y=10*np.random.rand(100)
        xpoint=np.mean(x)
        ypoint=np.mean(y)
        indxknn,distknn=sm.findknn(x,y,xpoint,ypoint,1,'1d','yes')

        x=10*np.random.rand(100)
        y=10*np.random.rand(100)
        xpoint=[2.5,5,7.5]
        ypoint=[3,6,9]
        indxknn,distknn=sm.findknn(x,y,xpoint,ypoint,10,'1d','yes')

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
    if dispout=='yes':
        import matplotlib.pyplot as plt 

    #--------------------------------------------------------------------------
    #Convert inputs to numpy array

    #Changing type to numpy array
    def type2numpy(variable):
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
    xpoint=type2numpy(xpoint)
    ypoint=type2numpy(ypoint)

    #--------------------------------------------------------------------------
    #Calculating distance from (x,y) to (xpoint,ypoint)

    M=len(x)
    N=len(xpoint)
    
    #Calculating distance
    #First method, using 1d array
    if distCalcMethod=='1d':
    
        #indxknn=np.ones((N,numofneighbors)) #Pre-assigning array
        #distknn=np.ones((N,numofneighbors)) #Pre-assigning array
        distxy=np.zeros((M,N)) #Pre-assigning array
        for i in range (0,N,1):
            distxy[:,i]=np.sqrt((x-xpoint[i])**2+(y-ypoint[i])**2) #Calculating distance from (x,y) to (xpoint,ypoint)
            #distxysort=np.sort(distxy) #Sorting distances
            #Indx=np.argsort(distxy) #Sorting distances
            #xsort=x[Indx] #Sorting X based on distance to (xpoint,ypoint)
            #ysort=y[Indx] #Sorting Y based on distance to (xpoint,ypoint)
            #indxknn[i,:]=Indx[0:numofneighbors] #Storing index of first numofneighbors points
            #distknn[i,:]=DistSort[0:numofneighbors] #Storing distance of first numofneighbors points
    
    #Second method, using pdist2 function
    elif distCalcMethod=='pdist2':
    
        X=np.column_stack((x,y))
        Y=np.column_stack((xpoint,ypoint))
        distxy=sp.spatial.distance.cdist(X,Y) #Calculating distance from (x1,y1) to (x2,y2)
    
    #Third method, using vectorized distance 
    elif distCalcMethod=='vector':
    
        X=np.column_stack((x,y))
        Y=np.column_stack((xpoint,ypoint))
        
        X2=np.sum(X**2,axis=1)
        Y2=np.sum(Y**2,axis=1)
        distxy=np.sqrt(np.ones((M,1))*np.transpose(Y2)+np.dot(X2[:,None],(np.ones((1,N))))-2*np.dot(X,np.transpose(Y))) #Calculating distance from (x1,y1) to (x2,y2)

    #--------------------------------------------------------------------------
    #Finding k-nearest neighbors

    if (numofneighbors<1): numofneighbors=1 #Minimum value for numofneighbors is one
    
    if numofneighbors==1:
        distxymin=np.min(distxy,axis=0) #Finding minimum distance
        Indx=np.argmin(distxy,axis=0) #Finding minimum distance
        indxknn=Indx.copy() #Storing index of first numofneighbors points
        distknn=distxymin.copy() #Storing distance of first numofneighbors points
    else:
        distxysort=np.sort(distxy,axis=0) #Finding minimum distance
        Indx=np.argsort(distxy,axis=0) #Finding minimum distance
        indxknn=np.transpose(Indx[0:numofneighbors,:]) #Storing index of first numofneighbors points
        distknn=np.transpose(distxysort[0:numofneighbors,:]) #Storing distance of first numofneighbors points

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
    
        #Plotting data
        plt.scatter(x,y,label='Data')
        plt.scatter(xpoint,ypoint,label='Point')
        if numofneighbors==1:
            plt.scatter(x[np.int64(indxknn)],y[np.int64(indxknn)])
        else:
            for i in range(0,N,1):
                plt.scatter(x[np.int64(indxknn[i,0:numofneighbors])],y[np.int64(indxknn[i,0:numofneighbors])])
        
        plt.xlabel('x')
        plt.ylabel('y')
        plt.legend()
        
    #--------------------------------------------------------------------------
    #Outputs
    return indxknn, distknn

    #--------------------------------------------------------------------------
