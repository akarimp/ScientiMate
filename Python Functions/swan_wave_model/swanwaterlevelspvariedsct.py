def swanwaterlevelspvariedsct(xgrid, ygrid, xpoint, ypoint, waterlevel, savedata='no', outfilename='swanwaterlevel.wl', outfilelocation=None, CalcMethod='linear'):
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

    scientimate.swanwaterlevelspvariedsct
    =====================================

    .. code:: python

        swanwaterlevel, waterlevelpoint = scientimate.swanwaterlevelspvariedsct(xgrid, ygrid, xpoint, ypoint, waterlevel, savedata='no', outfilename='swanwaterlevel.wl', outfilelocation=None, CalcMethod='linear')

    Description
    -----------

    | Generate SWAN water level file for spatially varied water level from scattered input data
    | This function can be used for any other scalar variable as well

    Inputs
    ------

    xgrid
        x (longitude) of output grid points as a [M*N] array
    ygrid
        y (latitude) of output grid points as a [M*N] array
    xpoint
        | x (longitude) of the location (such as meteorological station) that water level (or any scalar variable) is known in that location
        | as a [L] array
        | 1st element is 1st point, 2nd element is 2nd point, 3rd element is 3rd point, ...
    ypoint
        | y (latitude) of the location (such as meteorological station) that water level (or any scalar variable) is known in that location
        | as a [L] array
        | 1st element is 1st point, 2nd element is 2nd point, 3rd element is 3rd point, ...
    waterlevel
        | Water level (or any scalar variable) at (xpoint,ypoint) as a [K*L] array
        | K is number of time steps for a time series
        |     1st row (i.e. [1,:]) is 1st time step, 2nd row (i.e. [2,:]) is 2nd time step, ...
        |     Example: waterlevel(1,:) is the water level data for all points at the first time step
        | L is number of points (such as meteorological stations) that water level (or any scalar variable) is known in those locations
        |     L should be L>=3
        |     1st column (i.e. [:,1]) is 1st point, 2nd column (i.e. [:,2]) is 2nd point, 3rd column (i.e. [:,3]) is 3rd point, ...
        |     Example: waterlevel(:,1) is the water level data for all time steps at the first point
    savedata='no'
        | Define if save data in a file or not
        | 'no': does not save 
        | 'yes': save data as ascii file
    outfilename='swanwaterlevel.wl'
        | Name of output file between ' ' mark, example: 'swanwaterlevel.wl'
        | outfilename should have '.wl' extension
        | In case of using other scalar varable than water level, use proper name and extension
    outfilelocation=pwd
        Location of output file between ' ' mark, example: 'C:\' in MATLAB, or 'C:/' in Python
    CalcMethod='linear'
        | Interpolation method 
        | 'linear': Use default or 'linear' method to interpolate
        | 'nearest': Use nearest neighbor method to interpolate

    Outputs
    -------

    swanwaterlevel
        | Spatially varied water level data (or any scalar variable) formated for SWAN
        | Water level data at each time step is assigned into the grid points
    waterlevelpoint
        | Nearest interpolated water level (or any scalar variable) at (xpoint,ypoint) as a [K*L] array
        | K is number of time steps for a time series
        | L is number of points (such as meteorological stations) that water level (or any scalar variable) is known in those locations

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        xgrid,ygrid=np.meshgrid(np.linspace(-91,-90,100),np.linspace(28,30,100))
        xpoint=[-90.5,-90.3,-90.7] #Data are known at 3 locations
        ypoint=[29.2,29,28.8] #Data are known at 3 locations
        waterlevel=[[0.5,0.55,0.6],[0.64,0.69,0.74],[0.7,0.75,0.8],[0.4,0.45,0.5]] #Data for 4 time steps
        savedata='no'
        outfilename='swanwaterlevel.wl'
        outfilelocation=None
        CalcMethod='linear'
        swanwaterlevel,waterlevelpoint=sm.swanwaterlevelspvariedsct(xgrid,ygrid,xpoint,ypoint,waterlevel,savedata,outfilename,outfilelocation,CalcMethod)


        xgrid,ygrid=np.meshgrid(np.linspace(-91,-90,100),np.linspace(28,30,100))
        xpoint=[-90.5,-90.3,-90.7] #Data are known at 3 locations
        ypoint=[29.2,29,28.8] #Data are known at 3 locations
        waterlevel=[0.5,0.55,0.6] #Data for 1 time step
        savedata='no'
        outfilename='swanwaterlevel.wl'
        outfilelocation=None
        CalcMethod='linear'
        swanwaterlevel,waterlevelpoint=sm.swanwaterlevelspvariedsct(xgrid,ygrid,xpoint,ypoint,waterlevel,savedata,outfilename,outfilelocation,CalcMethod)

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
    waterlevel=type2numpy(waterlevel)

    #--------------------------------------------------------------------------
    #Assign default values

    if outfilelocation is None: outfilelocation=os.getcwd()

    #--------------------------------------------------------------------------
    #Generating SWAN water level file (.wl)

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
    #Then SWAN water level file (.wl) for waterlevel with three time steps has a format of (RowColumn,TimeStep):
    # wl(11,1) wl(12,1) wl(13,1) wl(14,1)
    # wl(21,1) wl(22,1) wl(23,1) wl(24,1)
    # wl(31,1) wl(32,1) wl(33,1) wl(34,1)
    # wl(41,1) wl(42,1) wl(43,1) wl(44,1)
    # wl(11,2) wl(12,2) wl(13,2) wl(14,2)
    # wl(21,2) wl(22,2) wl(23,2) wl(24,2)
    # wl(31,2) wl(32,2) wl(33,2) wl(34,2)
    # wl(41,2) wl(42,2) wl(43,2) wl(44,2)
    # wl(11,3) wl(12,3) wl(13,3) wl(14,3)
    # wl(21,3) wl(22,3) wl(23,3) wl(24,3)
    # wl(31,3) wl(32,3) wl(33,3) wl(34,3)
    # wl(41,3) wl(42,3) wl(43,3) wl(44,3)
    
    M,N=np.shape(xgrid)

    #waterlevel is 2d array
    if np.ndim(waterlevel)==2: 
        K,L=np.shape(waterlevel)

    #waterlevel is 1d array
    elif np.ndim(waterlevel)==1: 
        #If waterlevel has 1 time step, then it should be 1d-row array
        K=1 #K=1 when waterlevel is 1d array
        L=(np.shape(waterlevel))[0]

    #Checking if waterlevel is 1d-row (1 time step)
    #If waterlevel has 1 time step, then it should be 1d-row array
    if np.ndim(waterlevel)==1: 
        waterlevel=np.reshape(waterlevel,(1,np.size(waterlevel)))
    elif ((np.ndim(waterlevel)==2) and ((K==1) or (L==1))):
        waterlevel=np.reshape(waterlevel,(1,np.size(waterlevel)))
    
    waterlevelpoint=np.zeros((K,L)) #Pre-assigning an array

    #K is number of time steps for a time series
    for i in range(0,K,1):
        
        waterleveltimestep=waterlevel[i,:]
        
        #Interpolating data into grid using default or linear method
        if CalcMethod=='linear':
            zgrid=sp.interpolate.griddata((xpoint,ypoint),waterleveltimestep,(xgrid,ygrid))
        
            #Replacing NaN data point resulted from default method with ones from nearest method
            if np.sum(np.isnan(zgrid))>0:
                zgridnearest=sp.interpolate.griddata((xpoint,ypoint),waterleveltimestep,(xgrid,ygrid),method='nearest')
                zgrid[np.isnan(zgrid)==True]=zgridnearest[np.isnan(zgrid)==True] 
        
        #Interpolating data into grid using nearest neighbor method
        elif CalcMethod=='nearest':
            zgrid=sp.interpolate.griddata((xpoint,ypoint),waterleveltimestep,(xgrid,ygrid),method='nearest')

        #Creating the water level
        if i==0:
            swanwaterlevel=zgrid.copy()
        else:
            swanwaterlevel=np.concatenate((swanwaterlevel,zgrid))
        
        #Reshaping array into [:,1]
        xgridrshp=np.reshape(xgrid,(M*N,1)) #Reshape xgrid into [:,1]
        ygridrshp=np.reshape(ygrid,(M*N,1)) #Reshape ygrid into [:,1]
        zgridrshp=np.reshape(zgrid,(M*N,1)) #Reshape zgrid (waterlevelgrid) into [:,1]
        
        #Export nearest interpolated water level at each station location
        for j in range(0,L,1):
            distxy=np.sqrt((xgridrshp-xpoint[j])**2+(ygridrshp-ypoint[j])**2) #Calculating distance from (x,y) to (xpoint,ypoint)
            Indx=np.argmin(distxy,axis=0) #Finding minimum distance
            waterlevelpoint[i,j]=zgridrshp[Indx,0] #Nearest interpolated water level at each point of (xpoint,ypoint)
        

    #--------------------------------------------------------------------------
    #Saving data

    if savedata=='yes':
    
        #Changing directory to saving directory
        currentFolder=os.getcwd()
        os.chdir(outfilelocation)
    
        #Saving data
        np.savetxt(outfilename,swanwaterlevel,delimiter=' ')
    
        #Changing directory to working directory
        os.chdir(currentFolder)

    #--------------------------------------------------------------------------
    #Outputs
    return swanwaterlevel, waterlevelpoint

    #--------------------------------------------------------------------------
