def swanwindspvariedgrid(xgrid, ygrid, xpointgrid, ypointgrid, windvelgrid, winddirgrid, winddirtype='mete', windvelmin=0, savedata='no', outfilename='swanwind.wnd', outfilelocation=None, CalcMethod='linear'):
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

    scientimate.swanwindspvariedgrid
    ================================

    .. code:: python

        swanwind = scientimate.swanwindspvariedgrid(xgrid, ygrid, xpointgrid, ypointgrid, windvelgrid, winddirgrid, winddirtype='mete', windvelmin=0, savedata='no', outfilename='swanwind.wnd', outfilelocation=None, CalcMethod='linear')

    Description
    -----------

    Generate SWAN wind file for spatially varied wind from gridded input data

    Inputs
    ------

    xgrid
        x (longitude) of output grid points as a [M*N] array
    ygrid
        y (latitude) of output grid points as a [M*N] array
    xpointgrid
        x (longitude) of the locations that wind is known in those locations as a [K*L] array
    ypointgrid
        y (latitude) of the locations that wind is known in those locations as a [K*L] array
    windvelgrid
        | Wind velocity at (xpointgrid,ypointgrid) as a [K*L*P] array
        | P is number of time steps for a time series
    winddirgrid
        | Wind direction at (xpointgrid,ypointgrid) as a [K*L*P] array in (Degree)
        | P is number of time steps for a time series
    winddirtype='mete'
        | Define wind direction type
        | 'mete': meteorological wind direction 
        |     Meteorological direction represents a direction wind comes from and is measured counter-clockwise from the North
        |     0 (degree): from North, 90 (degree): from East, 180 (degree): from South, 270 (degree): from West
        | 'trig': trigonometric wind direction
    windvelmin=0
        Minimum allowed wind velocity
    savedata='no'
        | Define if save data in a file or not
        | 'no': does not save 
        | 'yes': save data as ascii file
    outfilename='swanwind.wnd'
        | Name of output file between ' ' mark, example: 'swanwind.wnd'
        | outfilename should have '.wnd' extension
    outfilelocation=pwd
        Location of output file between ' ' mark, example: 'C:\' in MATLAB, or 'C:/' in Python
    CalcMethod='linear'
        | Interpolation method 
        | 'linear': Use default or 'linear' method to interpolate
        | 'nearest': Use nearest neighbor method to interpolate

    Outputs
    -------

    swanwind
        | Spatially varied wind velocity data formated for SWAN
        | Wind velocity data at each time step is assigned into the grid points

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        xgrid,ygrid=np.meshgrid(np.linspace(-91,-90,100),np.linspace(28,30,100))
        xpointgrid,ypointgrid=np.meshgrid(np.linspace(-92,-89,100),np.linspace(27,31,100))
        windvelgrid=10+(12-10)*np.random.rand(100,100,4) #Data for 4 time steps
        winddirgrid=60+(65-60)*np.random.rand(100,100,4) #Data for 4 time steps
        winddirtype='mete'
        windvelmin=2.5
        savedata='no'
        outfilename='swanwind.wnd'
        outfilelocation=None
        CalcMethod='linear'
        swanwind=sm.swanwindspvariedgrid(xgrid,ygrid,xpointgrid,ypointgrid,windvelgrid,winddirgrid,winddirtype,windvelmin,savedata,outfilename,outfilelocation,CalcMethod)


        xgrid,ygrid=np.meshgrid(np.linspace(-91,-90,100),np.linspace(28,30,100))
        xpointgrid,ypointgrid=np.meshgrid(np.linspace(-92,-89,100),np.linspace(27,31,100))
        windvelgrid=10+(12-10)*np.random.rand(100,100) #Data for 1 time step
        winddirgrid=60+(65-60)*np.random.rand(100,100) #Data for 1 time step
        winddirtype='mete'
        windvelmin=2.5
        savedata='no'
        outfilename='swanwind.wnd'
        outfilelocation=None
        CalcMethod='linear'
        swanwind=sm.swanwindspvariedgrid(xgrid,ygrid,xpointgrid,ypointgrid,windvelgrid,winddirgrid,winddirtype,windvelmin,savedata,outfilename,outfilelocation,CalcMethod)

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
    xpointgrid=type2numpy(xpointgrid)
    ypointgrid=type2numpy(ypointgrid)
    windvelgrid=type2numpy(windvelgrid)
    winddirgrid=type2numpy(winddirgrid)

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
    #Generating SWAN wind file (.wnd)

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
    
    #Checking minimum value of wind velocity
    if (windvelmin<0): windvelmin=0
    
    M,N=np.shape(xgrid)

    #windvelgrid is 3d array
    if np.ndim(windvelgrid)==3: 
        K,L,P=np.shape(windvelgrid)

    #windvelgrid is 2d array
    elif np.ndim(windvelgrid)==2: 
        K,L=np.shape(windvelgrid)
        P=1
    
    #Checking if windvelgrid is 2d
    #If windvelgrid has 1 time step, then it should be 2d array of [K*L*1]
    if np.ndim(windvelgrid)==2: 
        windvelgrid=np.reshape(windvelgrid,(K,L,1))

    #Checking if winddirgrid is 2d
    #If winddirgrid has 1 time step, then it should be 2d array of [K*L*1]
    if np.ndim(winddirgrid)==2: 
        winddirgrid=np.reshape(winddirgrid,(K,L,1))

    #Reshaping array into [K*L,1]
    xpoint=np.reshape(xpointgrid,(K*L))
    ypoint=np.reshape(ypointgrid,(K*L))
    
    #P is number of time steps for a time series
    for i in range(0,P,1):
        
        #Reshaping array into [K*L,1]
        windveltimestep=np.reshape(windvelgrid[:,:,i],(K*L))
        winddirtimestep=np.reshape(winddirgrid[:,:,i],(K*L))
        
        #Checking for the minimum wind velocity
        windveltimestep[windveltimestep<windvelmin]=windvelmin
        
        #Meteorological direction
        if winddirtype=='mete':
            #Convert wind data from wind source to wind destination
            #Converting meteorological direction to trigonometric direction
            winddirtimesteptrig=-winddirtimestep+270
    
            #Add 360 to all numbers to have them all positive
            #Use mod(360) to take care of the ones larger than 360
            winddirtimesteptrig=((winddirtimesteptrig+360)%360) 
    
        #Trigonometric direction
        elif winddirtype=='trig':
            winddirtimesteptrig=winddirtimestep.copy()
        
    
        #Converting wind velocity into x and y components
        winddirtimesteptrigrad=np.deg2rad(winddirtimesteptrig) #Converting wind direction to Radian
        windveltimestepx,windveltimestepy=pol2cart(winddirtimesteptrigrad,windveltimestep) #Calculating x and y components of wind velocity
        
        #x direction
        #Interpolating data into grid using default or linear method
        if CalcMethod=='linear':
            windvelgridx=sp.interpolate.griddata((xpoint,ypoint),windveltimestepx,(xgrid,ygrid))
        
            #Replacing NaN data point resulted from default method with ones from nearest method
            if np.sum(np.isnan(windvelgridx))>0:
                windvelgridxnearest=sp.interpolate.griddata((xpoint,ypoint),windveltimestepx,(xgrid,ygrid),method='nearest')
                windvelgridx[np.isnan(windvelgridx)==True]=windvelgridxnearest[np.isnan(windvelgridx)==True] 
        
        #Interpolating data into grid using nearest neighbor method
        elif CalcMethod=='nearest':
            windvelgridx=sp.interpolate.griddata((xpoint,ypoint),windveltimestepx,(xgrid,ygrid),method='nearest')
    
        #y direction
        #Interpolating data into grid using default or linear method
        if CalcMethod=='linear':
            windvelgridy=sp.interpolate.griddata((xpoint,ypoint),windveltimestepy,(xgrid,ygrid))
        
            #Replacing NaN data point resulted from default method with ones from nearest method
            if np.sum(np.isnan(windvelgridy))>0:
                windvelgridynearest=sp.interpolate.griddata((xpoint,ypoint),windveltimestepy,(xgrid,ygrid),method='nearest')
                windvelgridy[np.isnan(windvelgridy)==True]=windvelgridynearest[np.isnan(windvelgridy)==True] 
        
        #Interpolating data into grid using nearest neighbor method
        elif CalcMethod=='nearest':
            windvelgridy=sp.interpolate.griddata((xpoint,ypoint),windveltimestepy,(xgrid,ygrid),method='nearest')
    
        #Creating the wind
        windgrid=np.concatenate((windvelgridx,windvelgridy))
        if i==0:
            swanwind=windgrid.copy()
        else:
            swanwind=np.concatenate((swanwind,windgrid))
        

    #--------------------------------------------------------------------------
    #Saving data

    if savedata=='yes':
    
        #Changing directory to saving directory
        currentFolder=os.getcwd()
        os.chdir(outfilelocation)
    
        #Saving data
        np.savetxt(outfilename,swanwind,delimiter=' ')
    
        #Changing directory to working directory
        os.chdir(currentFolder)

    #--------------------------------------------------------------------------
    #Outputs
    return swanwind

    #--------------------------------------------------------------------------
