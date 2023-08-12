def swandepthgrid(xgrid, ygrid, zgrid, zmin=None, zmax=None, nanreplacement='no', savedata='no', xyoutfilename='swangrid.xy', zoutfilename='swandepth.dep', outfilelocation=None, dispout='no'):
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

    scientimate.swandepthgrid
    =========================

    .. code:: python

        swandepth, swangrid, ncellx, ncelly, ngridx, ngridy = scientimate.swandepthgrid(xgrid, ygrid, zgrid, zmin=None, zmax=None, nanreplacement='no', savedata='no', xyoutfilename='swangrid.xy', zoutfilename='swandepth.dep', outfilelocation=None, dispout='no')

    Description
    -----------

    Generate SWAN depth file and its associated x-y grid file

    Inputs
    ------

    xgrid
        x (longitude) of grid points as a [M*N] array
    ygrid
        y (latitude) of grid points as a [M*N] array
    zgrid
        | z (elevation) of grid points as a [M*N] array
        | z>0 is land, z<0 is water, z=0 is water surface
    zmin=nanmin(zgrid)
        | Minimum z (elevation) to be considered
        | All z<zmin would be set to NaN
    zmax=nanmax(zgrid)
        | Maximum z (elevation) to be considered
        | All z>zmax would be set to NaN
        | Example: to remove all land values from z data, set zmax=0
    nanreplacement='no'
        | Replace NaN values with nanreplacement
        | By seting up an exception value equal to nanreplacement in SWAN input file, SWAN disregards these data points
        | If nanreplacement='no', then it does not replace NaN values
        | Example, nanreplacement=-999 will replace all NaN values with -999
        |     Then if -999 is set as an exception in SWAN input file, SWAN disregards all -999 values
    savedata='no'
        | Define if save data in a file or not
        | 'no': does not save 
        | 'yes': save data as ascii file
    xyoutfilename='swangrid.xy'
        | Name of output grid file between ' ' mark, example: 'swangrid.xy'
        | xyoutfilename should have '.xy' extension
    zoutfilename='swandepth.dep'
        | Name of output depth file between ' ' mark, example: 'swandepth.dep'
        | zoutfilename should have '.dep' extension
    outfilelocation=pwd
        Location of output file between ' ' mark, example: 'C:\' in MATLAB, or 'C:/' in Python
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    swandepth
        | Depth data formated for SWAN
        | Note: SWAN (Delft3D) needs depth data not bathymetry data. 
        |     It means value below water surface should be positive and above surface negative.
        | Note: All NaN values are replaced with nanreplacement
        |     Set up an exception value equal to nanreplacement in SWAN to disregard these data points
    swangrid
        Grid data formated for SWAN
    ncellx
        | Number of cells in x direction (ncellx=ngridx-1)
        | In SWAN, ncellx is equal to a number of meshes in computational grid in x direction 
    ncelly
        | Number of cells in y direction (ncelly=ngridy-1)
        | In SWAN, ncelly is equal to a number of meshes in computational grid in y direction 
    ngridx
        Number of grid points in x direction
    ngridy
        Number of grid points in y direction

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np
        import scipy as sp
        from scipy import interpolate

        x=10.*np.random.rand(1000)
        y=10*np.random.rand(1000)
        z=x**2+y**2-np.mean(x**2+y**2)
        xgrid,ygrid=np.meshgrid(np.linspace(np.min(x),np.max(x),100),np.linspace(np.min(y),np.max(y),100))
        zgrid=sp.interpolate.griddata((x,y),z,(xgrid,ygrid))
        zmin=np.nanmin(zgrid)
        zmax=np.nanmax(zgrid)
        nanreplacement=-999
        savedata='no'
        xyoutfilename='swangrid.xy'
        zoutfilename='swandepth.dep'
        outfilelocation=None
        swandepth,swangrid,ncellx,ncelly,ngridx,ngridy=sm.swandepthgrid(xgrid,ygrid,zgrid,zmin,zmax,nanreplacement,savedata,xyoutfilename,zoutfilename,outfilelocation,'yes')

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
    if dispout=='yes':
        import matplotlib.pyplot as plt 
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
    zgrid=type2numpy(zgrid)

    #--------------------------------------------------------------------------
    #Assign default values

    if zmin is None: zmin=np.nanmin(zgrid)
    if zmax is None: zmax=np.nanmax(zgrid)
    if outfilelocation is None: outfilelocation=os.getcwd()

    #--------------------------------------------------------------------------
    #Checking input data

    #Makeing sure all z data are in range of [zmin,zmax]
    zgrid[zgrid<zmin]=zmin
    zgrid[zgrid>zmax]=zmax

    #--------------------------------------------------------------------------
    #Generating SWAN grid file (.xy)

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
    #Then SWAN x-y grid file (.xy) for this xgrid and ygrid has a format of:
    # x1 x2 x3 x4
    # x1 x2 x3 x4
    # x1 x2 x3 x4
    # x1 x2 x3 x4
    # y1 y1 y1 y1
    # y2 y2 y2 y2
    # y3 y3 y3 y3
    # y4 y4 y4 y4
    
    #Generating grid file for SWAN, X=xgrid, Y=ygrid
    swangrid=np.concatenate((xgrid,ygrid))

    #--------------------------------------------------------------------------
    #Generating SWAN depth file (.dep)

    #SWAN (Delft3D) needs depth data not bathymetry data. 
    #It means value below water surface should be positive and above surface negative.
    swandepth=-zgrid

    #--------------------------------------------------------------------------
    #Creating exception points by replacing NaN with nanreplacement 
    #SWAN disregards data point with value equal to nanreplacement

    if nanreplacement!='no':
        xgrid[np.isnan(xgrid)==True]=nanreplacement
        ygrid[np.isnan(ygrid)==True]=nanreplacement
        swandepth[np.isnan(swandepth)==True]=nanreplacement

    #--------------------------------------------------------------------------
    #Number of grid points and cells 

    #First method
    #ngridx: number of grid points in x direction
    #ngridy: number of grid points in y direction
    ngridy,ngridx=np.shape(xgrid)
    
    ncellx=ngridx-1 #Number of cells in x direction
    ncelly=ngridy-1 #Number of cells in y direction
    
    
    #Second method
    #ngridx=len(xgrid[0,:]) #Number of grid points in x direction (Number of cell in x direction is ngridx-1)
    #ngridy=len(ygrid[:,0]) #Number of grid points in y direction (Number of cell in y direction is ngridy-1)
    
    #ncellx=len(xgrid[0,:])-1 #Number of cells in x direction (Number of cell in x direction is ngridx-1)
    #ncelly=len(ygrid[:,0])-1 #Number of cells in y direction (Number of cell in y direction is ngridy-1)

    #--------------------------------------------------------------------------
    #Saving data

    if savedata=='yes':
    
        #Changing directory to saving directory
        currentFolder=os.getcwd()
        os.chdir(outfilelocation)
    
        #Saving data
        np.savetxt(xyoutfilename,swangrid,delimiter=' ')
        np.savetxt(zoutfilename,swandepth,delimiter=' ')
    
        #Changing directory to working directory
        os.chdir(currentFolder)

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
        
        plt.subplot(1,2,1)
        #plt.contourf(xgrid,ygrid,zgrid)
    
        #plt.pcolormesh(xgrid,ygrid,zgrid)
    
        plt.imshow(zgrid,extent=[np.nanmin(xgrid),np.nanmax(xgrid),np.nanmin(ygrid),np.nanmax(ygrid)],cmap=plt.cm.get_cmap(),aspect='auto',origin='lower')
    
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Bathymetry')
        plt.colorbar()
        
    
        plt.subplot(1,2,2)
        #plt.contourf(xgrid,ygrid,swandepth)
    
        #plt.pcolormesh(xgrid,ygrid,swandepth)
    
        plt.imshow(swandepth,extent=[np.nanmin(xgrid),np.nanmax(xgrid),np.nanmin(ygrid),np.nanmax(ygrid)],cmap=plt.cm.get_cmap(),aspect='auto',origin='lower')
    
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Depth')
        plt.colorbar()
    
    #--------------------------------------------------------------------------
    #Outputs
    return swandepth, swangrid, ncellx, ncelly, ngridx, ngridy

    #--------------------------------------------------------------------------
