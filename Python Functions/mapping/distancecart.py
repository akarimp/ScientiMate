def distancecart(x1, y1, x2, y2, CalcMethod='1d', dispout='no'):
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

    scientimate.distancecart
    ========================

    .. code:: python

        distxy, theta = scientimate.distancecart(x1, y1, x2, y2, CalcMethod='1d', dispout='no')

    Description
    -----------

    Calculate distance from (x1,y1) to (x2,y2) on cartesian coordinate

    Inputs
    ------

    x1
        x of start point (first point)
    y1
        y of start point (first point)
    x2
        x of end point (last point) 
    y2
        y of end point (last point) 
    CalcMethod='1d'
        | Distance calculation method 
        | '1d': use 1d array
        | 'pdist2': Use 2d distance function
        | 'vector': Use vectorized distance 
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    distxy
        | Distance from (x1,y1) to (x2,y2)
        | returns M*N array where M=length(x1) and N=length(x2)
        | mth row associated with mth point in (x,y)
        | nth column is associated with nth point in (x2,y2)
    theta
        | Angle from start point to end point in (Degree)
        | returns M*N array where M=length(x1) and N=length(x2)
        | mth row associated with mth point in (x,y)
        | nth column is associated with nth point in (x2,y2)

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        x1=10*np.random.rand(100)
        y1=10*np.random.rand(100)
        x2=[2.5,5,7.5]
        y2=[3,6,9]
        distxy,theta=sm.distancecart(x1,y1,x2,y2,'1d','yes')

        x1=10*np.random.rand(100)
        y1=10*np.random.rand(100)
        x2=100*np.random.rand(10)
        y2=100*np.random.rand(10)
        distxy,theta=sm.distancecart(x1,y1,x2,y2,'pdist2','yes')

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

    x1=type2numpy(x1)
    y1=type2numpy(y1)
    x2=type2numpy(x2)
    y2=type2numpy(y2)

    #--------------------------------------------------------------------------
    #Calculating distance from (x1,y1) to (x2,y2)

    M=len(x1)
    N=len(x2)
    
    #First method, using 1d array
    if CalcMethod=='1d':
    
        distxy=np.zeros((M,N)) #Pre-assigning array
        for i in range(0,N,1):
            distxy[:,i]=np.sqrt((x1-x2[i])**2+(y1-y2[i])**2) #Calculating distance from (x1,y1) to (x2,y2)
    
    #Second method, using pdist2 function
    elif CalcMethod=='pdist2':
    
        X=np.column_stack((x1,y1))
        Y=np.column_stack((x2,y2))
        distxy=sp.spatial.distance.cdist(X,Y) #Calculating distance from (x1,y1) to (x2,y2)
    
    #Third method, using vectorized distance 
    elif CalcMethod=='vector':
    
        X=np.column_stack((x1,y1))
        Y=np.column_stack((x2,y2))
        
        X2=np.sum(X**2,axis=1)
        Y2=np.sum(Y**2,axis=1)
        distxy=np.sqrt(np.ones((M,1))*np.transpose(Y2)+np.dot(X2[:,None],(np.ones((1,N))))-2*np.dot(X,np.transpose(Y))) #Calculating distance from (x1,y1) to (x2,y2)
        
    #--------------------------------------------------------------------------
    #Calculating angle of the line between start and end points

    thetarad=np.zeros((M,N)) #Pre-assigning an array
    for i in range(0,N,1):
        thetarad[:,i]=np.arctan2(y2[i]-y1,x2[i]-x1) #Angle in radian

    theta=np.rad2deg(thetarad) #Angle in degree

    #Add 360 to all numbers to have them all positive
    #Use mod(360) to take care of the ones larger than 360
    theta=((theta+360)%360) 

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
        
        #Plotting data
        plt.scatter(x1,y1,label='Start Point')
        plt.scatter(x2,y2,label='End Point')
        
        plt.xlabel('x')
        plt.ylabel('y')
        plt.legend()

    #--------------------------------------------------------------------------
    #Outputs
    return distxy, theta

    #--------------------------------------------------------------------------
