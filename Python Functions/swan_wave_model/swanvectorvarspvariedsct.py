def swanvectorvarspvariedsct(xgrid, ygrid, xpoint, ypoint, Vx, Vy, savedata='no', outfilename='swanwind.wnd', outfilelocation=None, CalcMethod='linear'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2017-12-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.swanvectorvarspvariedsct
    ====================================

    .. code:: python

        swanvectorvariable, Vxpoint, Vypoint = scientimate.swanvectorvarspvariedsct(xgrid, ygrid, xpoint, ypoint, Vx, Vy, savedata='no', outfilename='swanwind.wnd', outfilelocation=None, CalcMethod='linear')

    Description
    -----------

    Generate SWAN file for spatially varied vector variable from scattered input data

    Inputs
    ------

    xgrid
        x (longitude) of output grid points as a [M*N] array
    ygrid
        y (latitude) of output grid points as a [M*N] array
    xpoint
        | x (longitude) of the location (such as meteorological station) that vector variable is known in that location
        | as a [L] array
        | 1st element is 1st point, 2nd element is 2nd point, 3rd element is 3rd point, ...
    ypoint
        | y (latitude) of the location (such as meteorological station) that vector variable is known in that location
        | as a [L] array
        | 1st element is 1st point, 2nd element is 2nd point, 3rd element is 3rd point, ...
    Vx
        | Variable in x direction (x component of input variable) as a [K*L] array
        | K is number of time steps for a time series
        |     1st row (i.e. [1,:]) is 1st time step, 2nd row (i.e. [2,:]) is 2nd time step, ...
        |     Example: windvel(1,:) is the wind velocity data for all points at the first time step
        | L is number of points (such as meteorological stations) that wind velocity is known in those locations
        |     L should be L>=3
        |     1st column (i.e. [:,1]) is 1st point, 2nd column (i.e. [:,2]) is 2nd point, 3rd column (i.e. [:,3]) is 3rd point, ...
        |     Example: windvel(:,1) is the wind velocity data for all time steps at the first point
    Vy
        | Variable in y direction (y component of input variable) as a [K*L] array
        | K is number of time steps for a time series
        |     1st row (i.e. [1,:]) is 1st time step, 2nd row (i.e. [2,:]) is 2nd time step, ...
        |     Example: windvel(1,:) is the wind velocity data for all points at the first time step
        | L is number of points (such as meteorological stations) that wind velocity is known in those locations
        |     L should be L>=3
        |     1st column (i.e. [:,1]) is 1st point, 2nd column (i.e. [:,2]) is 2nd point, 3rd column (i.e. [:,3]) is 3rd point, ...
        |     Example: windvel(:,1) is the wind velocity data for all time steps at the first point
    savedata='no'
        | Define if save data in a file or not
        | 'no': does not save 
        | 'yes': save data as ascii file
    outfilename='swanwind.wnd'
        | Name of output file between ' ' mark, example: 'swanwind.wnd'
        | outfilename should have proper name and extension
    outfilelocation=pwd
        Location of output file between ' ' mark, example: 'C:\' in MATLAB, or 'C:/' in Python
    CalcMethod='linear'
        | Interpolation method 
        | 'linear': Use default or 'linear' method to interpolate
        | 'nearest': Use nearest neighbor method to interpolate

    Outputs
    -------

    swanvectorvariable
        | Spatially varied vector variable data formated for SWAN
        | Vector variable data at each time step is assigned into the grid points
    Vxpoint
        | Nearest interpolated vector variable in x direction at (xpoint,ypoint) as a [K*L] array
        | K is number of time steps for a time series
        | L is number of points (such as meteorological stations) that vector variable in x direction is known in those locations
    Vypoint
        | Nearest interpolated vector variable in y direction at (xpoint,ypoint) as a [K*L] array
        | K is number of time steps for a time series
        | L is number of points (such as meteorological stations) that vector variable in y direction is known in those locations

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        xgrid,ygrid=np.meshgrid(np.linspace(-91,-90,100),np.linspace(28,30,100))
        windvelx=[[10.5,10.55,10.6],[10.64,10.69,10.74],[10.7,10.75,10.8],[10.4,10.45,10.5]] #Data for 4 time steps
        windvely=[[1.5,1.55,1.6],[1.64,1.69,1.74],[1.7,1.75,1.8],[1.4,1.45,1.5]] #Data for 4 time steps
        xpoint=[-90.5,-90.3,-90.7] #Data are known at 3 locations
        ypoint=[29.2,29,28.8] #Data are known at 3 locations
        savedata='no'
        outfilename='swanwind.wnd'
        outfilelocation=None
        CalcMethod='linear'
        swanvectorvariable,windvelxpoint,windvelypoint=sm.swanvectorvarspvariedsct(xgrid,ygrid,xpoint,ypoint,windvelx,windvely,savedata,outfilename,outfilelocation,CalcMethod)


        xgrid,ygrid=np.meshgrid(np.linspace(-91,-90,100),np.linspace(28,30,100))
        windvelx=[10.5,10.55,10.6] #Data for 1 time step
        windvely=[1.5,1.55,1.6] #Data for 1 time step
        xpoint=[-90.5,-90.3,-90.7] #Data are known at 3 locations
        ypoint=[29.2,29,28.8] #Data are known at 3 locations
        savedata='no'
        outfilename='swanwind.wnd'
        outfilelocation=None
        CalcMethod='linear'
        swanvectorvariable,windvelxpoint,windvelypoint=sm.swanvectorvarspvariedsct(xgrid,ygrid,xpoint,ypoint,windvelx,windvely,savedata,outfilename,outfilelocation,CalcMethod)

    References
    ----------

    Booij, N. R. R. C., Ris, R. C., & Holthuijsen, L. H. (1999). 
    A third‐generation wave model for coastal regions: 1. Model description and validation. 
    Journal of geophysical research: Oceans, 104(C4), 7649-7666.

    SWAN Team. (2007). S
    WAN user manual. 
    Delft University of Technology. The Netherlands.

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
    import os

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
    xpoint=type2numpy(xpoint)
    ypoint=type2numpy(ypoint)
    Vx=type2numpy(Vx)
    Vy=type2numpy(Vy)

    #--------------------------------------------------------------------------
    #Assign default values

    if outfilelocation is None: outfilelocation=os.getcwd()

    #--------------------------------------------------------------------------
    #Defining required functions

    def cart2pol(x,y):
        rho=np.sqrt(x**2+y**2)
        theta=np.arctan2(y,x)
        return theta,rho
    
    def pol2cart(theta,rho):
        x=rho*np.cos(theta)
        y=rho*np.sin(theta)
        return x,y

    #--------------------------------------------------------------------------
    #Generating SWAN file

    #If 
    #xgrid has a format of:
    # x1 x2 x3 x4
    # x1 x2 x3 x4
    # x1 x2 x3 x4
    # x1 x2 x3 x4
    #ygrid gas a format of:
    # y1 y1 y1 y1
    # y2 y2 y2 y2
    # y3 y3 y3 y3
    # y4 y4 y4 y4
    #Then SWAN wind file (.wnd) for wind with three time steps has a format of (RowColumn,TimeStep):
    # U(11,1) U(12,1) U(13,1) U(14,1)
    # U(21,1) U(22,1) U(23,1) U(24,1)
    # U(31,1) U(32,1) U(33,1) U(34,1)
    # U(41,1) U(42,1) U(43,1) U(44,1)
    # V(11,1) V(12,1) V(13,1) V(14,1)
    # V(21,1) V(22,1) V(23,1) V(24,1)
    # V(31,1) V(32,1) V(33,1) V(34,1)
    # V(41,1) V(42,1) V(43,1) V(44,1)
    # U(11,2) U(12,2) U(13,2) U(14,2)
    # U(21,2) U(22,2) U(23,2) U(24,2)
    # U(31,2) U(32,2) U(33,2) U(34,2)
    # U(41,2) U(42,2) U(43,2) U(44,2)
    # V(11,2) V(12,2) V(13,2) V(14,2)
    # V(21,2) V(22,2) V(23,2) V(24,2)
    # V(31,2) V(32,2) V(33,2) V(34,2)
    # V(41,2) V(42,2) V(43,2) V(44,2)
    # U(11,3) U(12,3) U(13,3) U(14,3)
    # U(21,3) U(22,3) U(23,3) U(24,3)
    # U(31,3) U(32,3) U(33,3) U(34,3)
    # U(41,3) U(42,3) U(43,3) U(44,3)
    # V(11,3) V(12,3) V(13,3) V(14,3)
    # V(21,3) V(22,3) V(23,3) V(24,3)
    # V(31,3) V(32,3) V(33,3) V(34,3)
    # V(41,3) V(42,3) V(43,3) V(44,3)
    
    M,N=np.shape(xgrid)

    #Vx is 2d array
    if np.ndim(Vx)==2: 
        K,L=np.shape(Vx)

    #Vx is 1d array
    elif np.ndim(Vx)==1: 
        #If Vx has 1 time step, then it should be 1d-row array
        K=1 #K=1 when Vx is 1d array
        L=(np.shape(Vx))[0]

    #Checking if Vx is 1d-row (1 time step)
    #If Vx has 1 time step, then it should be 1d-row array
    if np.ndim(Vx)==1: 
        Vx=np.reshape(Vx,(1,np.size(Vx)))
    elif ((np.ndim(Vx)==2) and ((K==1) or (L==1))):
        Vx=np.reshape(Vx,(1,np.size(Vx)))

    #Checking if Vy is 1d-row (1 time step)
    #If Vy has 1 time step, then it should be 1d-row array
    if np.ndim(Vy)==1: 
        Vy=np.reshape(Vy,(1,np.size(Vy)))
    elif ((np.ndim(Vy)==2) and ((K==1) or (L==1))):
        Vy=np.reshape(Vy,(1,np.size(Vy)))

    Vxpoint=np.zeros((K,L)) #Pre-assigning an array
    Vypoint=np.zeros((K,L)) #Pre-assigning an array
    
    #K is number of time steps for a time series
    for i in range(0,K,1):
        
        Vtimestepx=Vx[i,:]
        Vtimestepy=Vy[i,:]
        
        #x direction
        #Interpolating data into grid using default or linear method
        if CalcMethod=='linear':
            Vgridx=sp.interpolate.griddata((xpoint,ypoint),Vtimestepx,(xgrid,ygrid))
        
            #Replacing NaN data point resulted from default method with ones from nearest method
            if np.sum(np.isnan(Vgridx))>0:
                Vgridxnearest=sp.interpolate.griddata((xpoint,ypoint),Vtimestepx,(xgrid,ygrid),method='nearest')
                Vgridx[np.isnan(Vgridx)==True]=Vgridxnearest[np.isnan(Vgridx)==True] 
        
        #Interpolating data into grid using nearest neighbor method
        elif CalcMethod=='nearest':
            Vgridx=sp.interpolate.griddata((xpoint,ypoint),Vtimestepx,(xgrid,ygrid),method='nearest')

        #y direction
        #Interpolating data into grid using default or linear method
        if CalcMethod=='linear':
            Vgridy=sp.interpolate.griddata((xpoint,ypoint),Vtimestepy,(xgrid,ygrid))
        
            #Replacing NaN data point resulted from default method with ones from nearest method
            if np.sum(np.isnan(Vgridy))>0:
                Vgridynearest=sp.interpolate.griddata((xpoint,ypoint),Vtimestepy,(xgrid,ygrid),method='nearest')
                Vgridy[np.isnan(Vgridy)==True]=Vgridynearest[np.isnan(Vgridy)==True] 
        
        #Interpolating data into grid using nearest neighbor method
        elif CalcMethod=='nearest':
            Vgridy=sp.interpolate.griddata((xpoint,ypoint),Vtimestepy,(xgrid,ygrid),method='nearest')

        #Creating the wind
        Vgrid=np.concatenate((Vgridx,Vgridy))
        if i==0:
            swanvectorvariable=Vgrid.copy()
        else:
            swanvectorvariable=np.concatenate((swanvectorvariable,Vgrid))
        
        #Reshaping array into [:,1]
        xgridrshp=np.reshape(xgrid,(M*N,1)) #Reshape xgrid into [:,1]
        ygridrshp=np.reshape(ygrid,(M*N,1)) #Reshape ygrid into [:,1]
        Vgridxrshp=np.reshape(Vgridx,(M*N,1)) #Reshape Vgridx into [:,1]
        Vgridyrshp=np.reshape(Vgridy,(M*N,1)) #Reshape Vgridy into [:,1]
    
        #Export interpolated velocity at each station location
        for j in range(0,L,1):
            distxy=np.sqrt((xgridrshp-xpoint[j])**2+(ygridrshp-ypoint[j])**2) #Calculating distance from (x,y) to (xpoint,ypoint)
            Indx=np.argmin(distxy,axis=0) #Finding minimum distance
            Vxpoint[i,j]=Vgridxrshp[Indx,0] #Nearest interpolated wind velocity in x direction at each point of (xpoint,ypoint)
            Vypoint[i,j]=Vgridyrshp[Indx,0] #Nearest interpolated wind velocity in y direction at each point of (xpoint,ypoint)
        

    #--------------------------------------------------------------------------
    #Saving data

    if savedata=='yes':
    
        #Changing directory to saving directory
        currentFolder=os.getcwd()
        os.chdir(outfilelocation)
    
        #Saving data
        np.savetxt(outfilename,swanvectorvariable,delimiter=' ')
    
        #Changing directory to working directory
        os.chdir(currentFolder)

    #--------------------------------------------------------------------------
    #Outputs
    return swanvectorvariable, Vxpoint, Vypoint

    #--------------------------------------------------------------------------
