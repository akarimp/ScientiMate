def interpxyz2grid(x, y, z, gridsize=100, gridsizetype='points', xmin=None, xmax=None, ymin=None, ymax=None, zmin=None, zmax=None, RetainRatio='all', interpMethod='nearest', dispout='no'):
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

    scientimate.interpxyz2grid
    ==========================

    .. code:: python

        xgrid, ygrid, zgrid = scientimate.interpxyz2grid(x, y, z, gridsize=100, gridsizetype='points', xmin=None, xmax=None, ymin=None, ymax=None, zmin=None, zmax=None, RetainRatio='all', interpMethod='nearest', dispout='no')

    Description
    -----------

    Interpolate x (longitude), y (latitude) and z (elevation) data into a defined mesh

    Inputs
    ------

    x
        x (longitude) data extracted from xyz file
    y
        y (latitude) data extracted from xyz file
    z
        z (elevation) data extracted from xyz file
    gridsize=100
        | Grid size in x (longitude) and y (latitude) directions to interpolate elevation data on them
        | if gridsizetype='length' then gridsize is a distance between grid points
        | if gridsizetype='points' then gridsize is number of grid points in each direction
    gridsizetype='points'
        | Grid size type 
        | 'points': gridsize is considered as number of grid points in each direction
        | 'length': gridsize is considered as length between grid points
    xmin=nanmin(x)
        Minimum x (longitude) of domain to be interpolated
    xmax=nanmax(x)
        Maximum x (longitude) of domain to be interpolated
    ymin=nanmin(y)
        Minimum y (latitude) of domain to be interpolated
    ymax=nanmax(y)
        Maximum y (latitude) of domain to be interpolated
    zmin=nanmin(z)
        | Minimum z (elevation) of domain to be interpolated
        | All z<zmin would be set to zmin
    zmax=nanmax(z)
        | Maximum z (elevation) of domain to be interpolated
        | All z>zmax would be set to zmax
    RetainRatio='all'
        | Define to down sample input data or not 
        | 'all': data are not down sampled
        | value between 0 and 1: percentage of retaining data
        | RetainRatio=0.8 : 80% of data are retained
    interpMethod='nearest'
        | Interpolation method
        | 'linear': Use default or 'linear' method to interpolate
        | 'nearest': Use nearest neighbor method to interpolate
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    xgrid
        Interpolated x (longitude) data on defined mesh
    ygrid
        Interpolated y (latitude) data on defined mesh
    zgrid
        Interpolated z (elevation) data on defined mesh

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        x=10.*np.random.rand(1000)
        y=10.*np.random.rand(1000)
        z=x**2+y**2
        xgrid,ygrid,zgrid=sm.interpxyz2grid(x,y,z,100,'points',np.nanmin(x),np.nanmax(x),np.nanmin(y),np.nanmax(y),np.nanmin(z),np.nanmax(z),'all','nearest','yes')

        x=(-90-(-91))*np.random.rand(1000)+(-91)
        y=(31-(30))*np.random.rand(1000)+(30)
        z=x**2+y**2
        xgrid,ygrid,zgrid=sm.interpxyz2grid(x,y,z,0.005,'length',np.nanmin(x),np.nanmax(x),np.nanmin(y),np.nanmax(y),np.nanmin(z),np.nanmax(z),'all','linear','yes')

    References
    ----------

    Geospatial data

    * https://www.mathworks.com/help/map/finding-geospatial-data.html
    * https://maps.ngdc.noaa.gov/viewers/wcs-client/
    * https://www.ngdc.noaa.gov/mgg/global/global.html
    * https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/
    * https://www.ngdc.noaa.gov/mgg/image/2minrelief.html
    * https://www.ngdc.noaa.gov/mgg/coastal/crm.html
    * https://viewer.nationalmap.gov/launch/
    * https://earthexplorer.usgs.gov
    * http://www.shadedrelief.com/cleantopo2/index.html

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
        import matplotlib as mpl
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
    
    x=type2numpy(x)
    y=type2numpy(y)
    z=type2numpy(z)

    #--------------------------------------------------------------------------
    #Assigning default values

    if xmin is None: xmin=np.nanmin(x)
    if xmax is None: xmax=np.nanmax(x)
    if ymin is None: ymin=np.nanmin(y)
    if ymax is None: ymax=np.nanmax(y)
    if zmin is None: zmin=np.nanmin(z)
    if zmax is None: zmax=np.nanmax(z)

    #--------------------------------------------------------------------------
    #Defining domain area and associated data

    #Assigning new x, y and z
    xdata=x.copy() #xdata
    ydata=y.copy() #ydata
    zdata=z.copy() #zdata
    
    #Defining domain
    #domain='all'
    #    | Define a domain to be interpolated
    #    | 'all': all data are interpolated 
    #    | 'domain': only data within a defined domain are interpolated
    if ((xmin==np.nanmin(x)) & (xmax==np.nanmax(x)) & (ymin==np.nanmin(y)) & (ymax==np.nanmax(y)) & (zmin==np.nanmin(z)) & (zmax==np.nanmax(z))):
        domain='all'
    else:
        domain='domain'

    #Defining domain boundaries
    if domain=='all':
    
        xmin=np.nanmin(x) 
        xmax=np.nanmax(x) 
        ymin=np.nanmin(y) 
        ymax=np.nanmax(y) 
        zmin=np.nanmin(z)
        zmax=np.nanmax(z)
    
    #Retaining data within a defined domain
    elif domain=='domain':
        
        #Finding a location data that are within [xmin:xmax] range
        xIndx=np.int_((np.nonzero((x>=xmin) & (x<=xmax)))[0])
        
        #Retaining data that are within [xmin:xmax] range
        x1=xdata[xIndx] #x data
        y1=ydata[xIndx] #y data
        z1=zdata[xIndx] #z data
        
        #Finding a location data that are within [ymin:ymax] range
        yIndx=np.int_((np.nonzero((y1>=ymin) & (y1<=ymax)))[0])
        
        #Retaining data that are within [ymin:ymax] range
        x2=x1[yIndx] #xdata
        y2=y1[yIndx] #ydata
        z2=z1[yIndx] #zdata
        
        #Finding a location data that are within [zmin:zmax] range
        zIndx=np.int_((np.nonzero((z2>=zmin) & (z2<=zmax)))[0])
        
        #Retaining data that are within [zmin:zmax] range
        x3=x2[zIndx] #xdata
        y3=y2[zIndx] #ydata
        z3=z2[zIndx] #zdata
    
        #Assigning new x, y and z
        xdata=x3.copy() #xdata
        ydata=y3.copy() #ydata
        zdata=z3.copy() #zdata
    
    #--------------------------------------------------------------------------
    #Down sampling the data

    if RetainRatio=='all':
        #Do nothing
        dummy=1
    
    elif str.isnumeric(str(RetainRatio))==1:
    
        #Defining index number to down sample data
        if ((RetainRatio<=0) or (RetainRatio>1)): RetainRatio=1
        M,N=np.shape(x)
        dsIndx=np.int64(np.linspace(0,M-1,int(M*RetainRatio)))
    
        x4=xdata.copy() #xdata
        y4=ydata.copy() #ydata
        z4=zdata.copy() #zdata
    
        #Retaining down sampled data
        xdata=x4[dsIndx] #xdata
        ydata=y4[dsIndx] #ydata
        zdata=z4[dsIndx] #zdata
    
    #--------------------------------------------------------------------------
    #Interpolating elevation data

    #Calculating grid size in each direction using grid length
    if gridsizetype=='length':
    
        #Checking x and y range to generate correct number of grid points
        if ((xmax-xmin)%gridsize)!=0:
            xmax=xmax-((xmax-xmin)%gridsize)
    
        if ((ymax-ymin)%gridsize)!=0:
            ymax=ymax-((ymax-ymin)%gridsize)
    
        gridsizex=gridsize #Grid size in lat (x) direction
        gridsizey=gridsize #Grid size in lon (y) direction
    
        ngridx=((xmax-xmin)/gridsizex)+1 #Number of grid points in each direction
        ngridy=((ymax-ymin)/gridsizey)+1 #Number of grid points in each direction

        #Generating grid
        xgrid,ygrid=np.meshgrid(np.arange(xmin,xmax+gridsizex,gridsizex),np.arange(ymin,ymax+gridsizey,gridsizey))

    #Calculating grid size in each direction using number grid points 
    elif gridsizetype=='points':
    
        ngridx=gridsize #Number of grid points in each direction
        ngridy=gridsize #Number of grid points in each direction

        gridsizex=(xmax-xmin)/(ngridx-1) #Grid size in x direction
        gridsizey=(ymax-ymin)/(ngridy-1) #Grid size in y direction
        
        #Generating grid
        xgrid,ygrid=np.meshgrid(np.linspace(xmin,xmax,ngridx),np.linspace(ymin,ymax,ngridy))
   
    #Interpolating data into grid
    zgrid=sp.interpolate.griddata((xdata,ydata),zdata,(xgrid,ygrid),method=interpMethod)

    #Replacing NaN data point resulted from interpolation by nearest point
    if np.sum(np.isnan(zgrid))>0:
        zgridnearest=sp.interpolate.griddata((xdata,ydata),zdata,(xgrid,ygrid),method='nearest')
        zgrid[np.isnan(zgrid)==True]=zgridnearest[np.isnan(zgrid)==True]
    
    
    #Makeing sure all z data are in range of [zmin,zmax]
    zgrid[zgrid<zmin]=zmin
    zgrid[zgrid>zmax]=zmax

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
        
        #plt.contourf(xgrid,ygrid,zgrid)
        
        #plt.pcolormesh(xgrid,ygrid,zgrid,shading='nearest')

        plt.imshow(zgrid,extent=[np.nanmin(xgrid),np.nanmax(xgrid),np.nanmin(ygrid),np.nanmax(ygrid)],interpolation=None,cmap=mpl.cm.get_cmap(),aspect='auto',origin='lower')

        plt.xlabel('x')
        plt.ylabel('y')
        plt.colorbar()
        
    #--------------------------------------------------------------------------
    #Outputs
    return xgrid, ygrid, zgrid

    #--------------------------------------------------------------------------
