def plot3dtopo(x, y, z, \
    zcolormap='topocmap', ncolor=256, waterlandcmapratio='none', \
    plotfontsize=12, axisfontsize=12, xaxislabel='x', yaxislabel='y', cbarlabel='z', dispcolorbar='no', dispgrid='no', \
    gridsize=100, gridsizetype='points', domain='all', xmin=None, xmax=None, ymin=None, ymax=None, zmin=None, zmax=None, RetainRatio='all', interpMethod='nearest', dispout='imagesc'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2022-05-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.plot3dtopo
    ======================

    .. code:: python

        xgrid, ygrid, zgrid, cmap_ncolor, cmap_ncolor_mpl = scientimate.plot3dtopo(x, y, z, \
            zcolormap='topocmap', ncolor=256, waterlandcmapratio='none', \
            plotfontsize=12, axisfontsize=12, xaxislabel='x', yaxislabel='y', cbarlabel='z', dispcolorbar='no', dispgrid='no', \
            gridsize=100, gridsizetype='points', domain='all', xmin=None, xmax=None, ymin=None, ymax=None, zmin=None, zmax=None, RetainRatio='all', interpMethod='nearest', dispout='imagesc')

    Description
    -----------

    Plot x (longitude), y (latitude) and z (elevation) data into a defined mesh

    Inputs
    ------

    x
        | x (longitude) data extracted from xyz file
        | Set x=[] if it is not available
        | It may be 1d or 2d array
    y
        | y (latitude) data extracted from xyz file
        | Set y=[] if it is not available
        | It may be 1d or 2d array
    z
        z (elevation) data extracted from xyz file
    zcolormap='topocmap'
        | Colormap for z data
        | Topographic (Water/Land) colormaps:
        | 'topocmap': colormap developed by Arash Karimpour
        | 'topocmaprelief': colormap developed by Arash Karimpour
        | 'topocmapocean': colormap developed by Arash Karimpour
        | 'topocmapearth': colormap developed by Arash Karimpour
        | 'blueocean': colormap developed by Arash Karimpour
        | 'blueoceansea': colormap developed by Arash Karimpour
        | 'greenearth': colormap developed by Arash Karimpour
        | 'greenearthland': colormap developed by Arash Karimpour
        | 'blgrtopocmap': colormap developed by Arash Karimpour
        | 'blrdtopocmap': colormap developed by Arash Karimpour
        | 'grayearth': colormap developed by Arash Karimpour
        | 'etopo1': ETOPO1 colormap, https://www.ngdc.noaa.gov/mgg/global/global.html
        | 'gmtglobe': GMT_globe colormap, https://www.giss.nasa.gov/tools/panoply/colorbars/
        | 'gmtrelief': GMT_relief colormap, https://www.giss.nasa.gov/tools/panoply/colorbars/
        | 'aendekerk': Colormap from  Florian Aendekerk, http://www.mathworks.com/matlabcentral/fileexchange/63590-landseacolormap-m-
        | Any other available color map such as 'cool', 'winter', etc can be used
        | Colormap can be defined by user as [n*3] array in RGB color format between 0 and 255
    ncolor=256
        Number of colors to be used in colormap
    waterlandcmapratio='none'
        | Scale of sea_colormap over land_colormap (between 0 to 1)
        | 'none': no scaling, original color map is used
        | 'auto': put z=0 at center of colormap
        | if equal to a number between 0 to 1: 
        | assigning ratio equal to waterlandcmapratio of the colormap to negative (water) values 
        | and the rest to positive (land) values 
        | waterlandcmapratio=0.55 : assign 55# of colormap to negative (water) values
    plotfontsize=12
        Size of plot fonts
    axisfontsize=12
        Size of axis fonts
    xaxislabel='x'
        x axis label
    yaxislabel='y'
        y axis label
    cbarlabel='z'
        Colorbar label
    dispcolorbar='no'
        Define to display colorbar or not ('yes': display, 'no': not display)
    dispgrid='no'
        | Define to display grid lines or not ('yes': display, 'no': not display)
        | 'yes': divide x axis and y axis to 13 intervals
        | if equal to a number:
        | divide x axis and y axis to dispgrid intervals  
        | dispgrid=10 : divide x axis and y axis to 10 intervals
    gridsize=100
        | Grid size in x (longitude) and y (latitude) directions to interpolate elevation data on them
        | if gridsizetype='length' then gridsize is a distance between grid points
        | if gridsizetype='points' then gridsize is number of grid points in each direction
    gridsizetype='points'
        | Grid size type 
        | 'number': gridsize is considered as number of grid points in each direction
        | 'length': gridsize is considered as length between grid points
    xmin=nanmin(x)
        Minimum x (longitude) of domain to be plotted
    xmax=nanmax(x)
        Maximum x (longitude) of domain to be plotted
    ymin=nanmin(y)
        Minimum y (latitude) of domain to be plotted
    ymax=nanmax(y)
        Maximum y (latitude) of domain to be plotted
    zmin=nanmin(z)
        | Minimum z (elevation) of domain to be plotted
        | All z<zmin would be set to zmin
    zmax=nanmax(z)
        | Maximum z (elevation) of domain to be plotted
        | All z>zmax would be set to zmax
    RetainRatio='all'
        | Define to down sample input data or not 
        | 'all': data are not down sampled
        | value between 0 and 1: percentage of retaining data
        | RetainRatio=0.8 : 80# of data are retained
    interpMethod='nearest'
        | Interpolation method 
        | 'linear': Use default or 'linear' method to interpolate
        | 'nearest': Use nearest neighbor method to interpolate
    dispout='imagesc'
        | Define to display outputs or not
        | 'imagesc': 2 dimensional plot using imagesc or imshow
        | 'pcolor': 2 dimensional plot using pcolor
        | 'contour': 2 dimensional contour plot, number of contour=ncolor
        | 'surface': 3 dimensional surface plot 
        | 'no': not display 

    Outputs
    -------

    xgrid
        Interpolated x (longitude) data on defined mesh
    ygrid
        Interpolated y (latitude) data on defined mesh
    zgrid
        Interpolated z (elevation) data on defined mesh
    cmap_ncolor
        | Colormap for z levels with ncolor number of colors in RGB color format between 0 and 1
        | To convert 0-1 scale tp 0-255 scale, multiply cmap_ncolor values by 255
    cmap_ncolor_mpl
        Colormap for z levels with ncolor number of colors in Matplotlib format

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        x=10*np.random.rand(1000)
        y=10*np.random.rand(1000)
        z=x**2+y**2
        xgrid,ygrid,zgrid,cmap_ncolor,cmap_ncolor_mpl=sm.plot3dtopo(x,y,z,\
            'topocmap',256,'none',\
            12,12,'x','y','z','no','no',\
            100,'points',np.nanmin(x),np.nanmax(x),np.nanmin(y),np.nanmax(y),np.nanmin(z),np.nanmax(z),'all','nearest','imagesc')

        x=10*np.random.rand(1000)
        y=10*np.random.rand(1000)
        z=x**2+y**2
        xgrid,ygrid=np.meshgrid(np.linspace(np.nanmin(x),np.nanmax(x),100),np.linspace(np.nanmin(y),np.nanmax(y),100))
        zgrid=sp.interpolate.griddata((x,y),z,(xgrid,ygrid))
        xgrid,ygrid,zgrid,cmap_ncolor,cmap_ncolor_mpl=sm.plot3dtopo(xgrid,ygrid,zgrid,\
            'topocmap',256,'none',\
            12,12,'x','y','z','no','no',\
            100,'points',np.nanmin(x),np.nanmax(x),np.nanmin(y),np.nanmax(y),np.nanmin(z),np.nanmax(z),'all','nearest','imagesc')

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

    Colormap

    * http://colorbrewer2.org
    * http://matplotlib.org/cmocean/
    * https://matplotlib.org/users/colormaps.html
    * http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml
    * https://www.giss.nasa.gov/tools/panoply/colorbars/
    * http://jdherman.github.io/colormap/

    Color

    * http://htmlcolorcodes.com

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
    import matplotlib as mpl 
    import matplotlib.pyplot as plt 
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import colors 

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
    zcolormap=type2numpy(zcolormap)

    #--------------------------------------------------------------------------
    #Assigning default values

    if xmin is None: xmin=np.nanmin(x)
    if xmax is None: xmax=np.nanmax(x)
    if ymin is None: ymin=np.nanmin(y)
    if ymax is None: ymax=np.nanmax(y)
    if zmin is None: zmin=np.nanmin(z)
    if zmax is None: zmax=np.nanmax(z)

    #--------------------------------------------------------------------------
    #Check input data

    #Create x if x=[]
    if np.size(x)==0:
        if np.ndim(z)==1:
            x=np.linspace(0,1,np.size(z))
        else:
            x,_=np.meshgrid(np.linspace(0,1,len(z[1,:])),np.linspace(0,1,len(z[:,1])))

    #Create y if y=[]
    if np.size(y)==0:
        if np.ndim(z)==1:
            y=np.linspace(0,1,np.size(z))
        else:
            _,y=np.meshgrid(np.linspace(0,1,len(z[1,:])),np.linspace(0,1,len(z[:,1])))


    #Input data shape
    if np.ndim(x)==1:
        Mx=np.shape(x)
        Nx=1   
    else:
        Mx,Nx=np.shape(x)

    if np.ndim(y)==1:
        My=np.shape(y)
        Ny=1   
    else:
        My,Ny=np.shape(y)

    if np.ndim(z)==1:
        Mz=np.shape(z)
        Nz=1   
    else:
        Mz,Nz=np.shape(z)


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
    
    elif ((isinstance(RetainRatio,int)==True) or (isinstance(RetainRatio,float)==True)):
    
        #Defining index number to down sample data
        if ((RetainRatio<=0) or (RetainRatio>1)): RetainRatio=1
        M=len(x)
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

    #Check if inputs are 1-d or 2-d arrays
    if (Mx==1 or Nx==1):

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

        #Calculating grid size in each direction using number grid points 
        elif gridsizetype=='points':
        
            ngridx=gridsize #Number of grid points in each direction
            ngridy=gridsize #Number of grid points in each direction

            gridsizex=(xmax-xmin)/(ngridx-1) #Grid size in x direction
            gridsizey=(ymax-ymin)/(ngridy-1) #Grid size in y direction
        
        #Generating grid
        #xgrid,ygrid=np.meshgrid(np.arange(xmin,xmax+gridsizex,gridsizex),np.arange(ymin,ymax+gridsizey,gridsizey))
        xgrid,ygrid=np.meshgrid(np.linspace(xmin,xmax,ngridx),np.linspace(ymin,ymax,ngridy))
        
        #Interpolating data into grid
        zgrid=sp.interpolate.griddata((xdata,ydata),zdata,(xgrid,ygrid),method=interpMethod)
    
        #Replacing NaN data point resulted from interpolation by nearest point
        if np.sum(np.isnan(zgrid))>0:
            zgridnearest=sp.interpolate.griddata((xdata,ydata),zdata,(xgrid,ygrid),method='nearest')
            zgrid[np.isnan(zgrid)==1]=zgridnearest[np.isnan(zgrid)==1]
          

    else:
        xgrid=x.copy()
        ygrid=y.copy()
        zgrid=z.copy()


    #Makeing sure all z data are in range of [zmin,zmax]
    zgrid[zgrid<zmin]=zmin
    zgrid[zgrid>zmax]=zmax

    #--------------------------------------------------------------------------
    #Colormaps

    #'topocmap' colormap
    #Colors from: http://htmlcolorcodes.com
    if zcolormap=='topocmap':
    
        #Water colormap
        topocmap_water=np.array([[21, 67, 96],[27, 79, 114],[40, 116, 166],[52, 152, 219],[133, 193, 233],[214, 234, 248]])/255
    
        #Land colormap
        topocmapgreen=np.array([[56, 142, 60]])/255 #Green: 1/30 land colormap
        topocmapyellow=np.array([[252, 243, 207],[249, 231, 159],[247, 220, 111],[241, 196, 15]])/255 #Yellow: 3/30 land colormap
        topocmapbrown=np.array([[243, 156, 18],[230, 126, 34],[211, 84, 0],[186, 74, 0],[160, 64, 0],[135, 54, 0],[110, 44, 0]])/255
        topocmapwhite=np.array([[236, 240, 241]])/255 #White: 1/30 land colormap
    
        #Interpolating brwon color from 7 colors to 24 colors
        topocmapbrown1=np.zeros((24,3))
        for i in range(0,3,1):
            topocmapbrown1[:,i]=np.interp(np.linspace(1,24,24),np.linspace(1,24,len(topocmapbrown[:,0])),topocmapbrown[:,i]) #Brown: 24/30 land colormap
    
        #Assigning water and land colormap
        cmap_water1=topocmap_water.copy()
        cmap_land1=np.vstack((topocmapgreen,topocmapyellow,topocmapbrown1,topocmapwhite))
    
        #Convert to 0.55# sea and 0.45# land 
        #Mariana Trench depth= 10994 m
        #Mount Everest height= 8848 m
        #(max depth)/(max depth+max height)~=0.55
        cmap_water=np.zeros((55,3))
        cmap_land=np.zeros((45,3))
        for i in range(0,3,1):
            cmap_water[:,i]=np.interp(np.linspace(1,55,55),np.linspace(1,55,len(cmap_water1[:,0])),cmap_water1[:,i])
            cmap_land[:,i]=np.interp(np.linspace(1,45,45),np.linspace(1,45,len(cmap_land1[:,0])),cmap_land1[:,i])
    
        #Assigning negative (water) and positive (land) colormap
        cmap_z=np.vstack((cmap_water,cmap_land))
    
    #'topocmaprelief' colormap
    #Colors from globalrelief.py
    elif zcolormap=='topocmaprelief':

        #Water colormap
        topocmap_water=mpl.cm.get_cmap('bone', 128)
        topocmap_water=topocmap_water(range(64,128,1))

        #Land colormap
        topocmapsummer=mpl.cm.get_cmap('summer', 10)
        topocmapsummer=topocmapsummer(range(0,10,1))

        topocmapcopper=mpl.cm.get_cmap('copper_r', 54)
        topocmapcopper=topocmapcopper(range(0,54,1))

        #Assigning water and land colormap
        cmap_water1=topocmap_water.copy()
        cmap_land1=np.vstack((topocmapsummer,topocmapcopper))

        #Convert to 0.55# sea and 0.45# land
        #Mariana Trench depth= 10994 m
        #Mount Everest height= 8848 m
        #(max depth)/(max depth+max height)~=0.55
        cmap_water=np.zeros((55,3))
        cmap_land=np.zeros((45,3))
        for i in range(0,3,1):
            cmap_water[:,i]=np.interp(np.linspace(1,55,55),np.linspace(1,55,len(cmap_water1[:,0])),cmap_water1[:,i])
            cmap_land[:,i]=np.interp(np.linspace(1,45,45),np.linspace(1,45,len(cmap_land1[:,0])),cmap_land1[:,i])

        #Assigning negative (water) and positive (land) colormap
        cmap_z=np.vstack((cmap_water,cmap_land))

    #'topocmapocean' colormap
    #Colors from: http://htmlcolorcodes.com
    elif zcolormap=='topocmapocean':
    
        #Water colormap
        topocmap_water=np.array([[21, 67, 96],[27, 79, 114],[40, 116, 166],[52, 152, 219],[133, 193, 233],[214, 234, 248]])/255
    
        cmap_z=topocmap_water.copy()
    
        #Assigning half of the colormap to negative (water) values and other half to positive (land) values 
        M1,N1=np.shape(cmap_z)
        cmap_water=cmap_z[0:int(M1/2),:]
        cmap_land=cmap_z[int(M1/2):,:]
    
    #'topocmapearth' colormap
    #Colors from: http://htmlcolorcodes.com
    elif zcolormap=='topocmapearth':
    
        #Land colormap
        topocmapgreen=np.array([[56, 142, 60]])/255 #Green: 1/30 land colormap
        topocmapyellow=np.array([[252, 243, 207],[249, 231, 159],[247, 220, 111],[241, 196, 15]])/255 #Yellow: 3/30 land colormap
        topocmapbrown=np.array([[243, 156, 18],[230, 126, 34],[211, 84, 0],[186, 74, 0],[160, 64, 0],[135, 54, 0],[110, 44, 0]])/255
        topocmapwhite=np.array([[236, 240, 241]])/255 #White: 1/30 land colormap
    
        #Interpolating brwon color from 7 colors to 24 colors
        topocmapbrown1=np.zeros((24,3))
        for i in range(0,3,1):
            topocmapbrown1[:,i]=np.interp(np.linspace(1,24,24),np.linspace(1,24,len(topocmapbrown[:,0])),topocmapbrown[:,i]) #Brown: 24/30 land colormap
    
        #Assigning water and land colormap
        cmap_land1=np.vstack((topocmapgreen,topocmapyellow,topocmapbrown1,topocmapwhite))
    
        cmap_z=cmap_land1.copy()
    
        #Assigning half of the colormap to negative (water) values and other half to positive (land) values 
        M1,N1=np.shape(cmap_z)
        cmap_water=cmap_z[0:int(M1/2),:]
        cmap_land=cmap_z[int(M1/2):,:]
    
    #'blueocean' colormap
    #Colors from: http://htmlcolorcodes.com
    elif zcolormap=='blueocean':
    
        #Water colormap
        blueoceanwater=np.array([[52, 152, 219],[133, 193, 233]])/255
    
        #Land colormap
        topocmapgreen=np.array([[56, 142, 60]])/255 #Green: 1/30 land colormap
        topocmapyellow=np.array([[252, 243, 207],[249, 231, 159],[247, 220, 111],[241, 196, 15]])/255 #Yellow: 3/30 land colormap
        topocmapbrown=np.array([[243, 156, 18],[230, 126, 34],[211, 84, 0],[186, 74, 0],[160, 64, 0],[135, 54, 0],[110, 44, 0]])/255
        topocmapwhite=np.array([[236, 240, 241]])/255 #White: 1/30 land colormap
    
        #Interpolating brwon color from 7 colors to 24 colors
        topocmapbrown1=np.zeros((24,3))
        for i in range(0,3,1):
            topocmapbrown1[:,i]=np.interp(np.linspace(1,24,24),np.linspace(1,24,len(topocmapbrown[:,0])),topocmapbrown[:,i]) #Brown: 24/30 land colormap
    
        #Assigning water and land colormap
        cmap_water1=blueoceanwater.copy()
        cmap_land1=np.vstack((topocmapgreen,topocmapyellow,topocmapbrown1,topocmapwhite))
    
        #Convert to 0.55# sea and 0.45# land 
        #Mariana Trench depth= 10994 m
        #Mount Everest height= 8848 m
        #(max depth)/(max depth+max height)~=0.55
        cmap_water=np.zeros((55,3))
        cmap_land=np.zeros((45,3))
        for i in range(0,3,1):
            cmap_water[:,i]=np.interp(np.linspace(1,55,55),np.linspace(1,55,len(cmap_water1[:,0])),cmap_water1[:,i])
            cmap_land[:,i]=np.interp(np.linspace(1,45,45),np.linspace(1,45,len(cmap_land1[:,0])),cmap_land1[:,i])
    
        #Assigning negative (water) and positive (land) colormap
        cmap_z=np.vstack((cmap_water,cmap_land))
    
    #'blueoceansea' colormap
    #Colors from: http://htmlcolorcodes.com
    elif zcolormap=='blueoceansea':
    
        #Water colormap
        blueoceanwater=np.array([[52, 152, 219],[133, 193, 233]])/255

        #At least 4 colors requred in colormap to guarantee at least two colors for each of cmap_water and cmap_land colormaps
        blueoceanwater4colors=np.zeros((4,3))
        for i in range(0,3,1):
            blueoceanwater4colors[:,i]=np.interp(np.linspace(1,4,4),np.linspace(1,4,len(blueoceanwater[:,0])),blueoceanwater[:,i])
    
        cmap_z=blueoceanwater4colors.copy()
    
        #Assigning half of the colormap to negative (water) values and other half to positive (land) values 
        M1,N1=np.shape(cmap_z)
        cmap_water=cmap_z[0:int(M1/2),:]
        cmap_land=cmap_z[int(M1/2):,:]
    
    #'greenearth' colormap
    #Colors from: http://htmlcolorcodes.com
    elif zcolormap=='greenearth':
    
        #Mariana Trench depth= 10994 m
        #Mount Everest height= 8848 m
        #(max depth)/(max depth+max height)~=0.55
    
        #Water colormap
        topocmap_water=np.array([[21, 67, 96],[27, 79, 114],[40, 116, 166],[52, 152, 219],[133, 193, 233],[214, 234, 248]])/255
    
        #Land colormap
        greenearthgreen=np.array([[67, 160, 71],[255, 241, 118]])/255 #Green
    
        #Assigning water and land colormap
        cmap_water1=topocmap_water.copy()
        cmap_land1=greenearthgreen.copy()
    
        #Convert to 0.55# sea and 0.45# land 
        #Mariana Trench depth= 10994 m
        #Mount Everest height= 8848 m
        #(max depth)/(max depth+max height)~=0.55
        cmap_water=np.zeros((55,3))
        cmap_land=np.zeros((45,3))
        for i in range(0,3,1):
            cmap_water[:,i]=np.interp(np.linspace(1,55,55),np.linspace(1,55,len(cmap_water1[:,0])),cmap_water1[:,i])
            cmap_land[:,i]=np.interp(np.linspace(1,45,45),np.linspace(1,45,len(cmap_land1[:,0])),cmap_land1[:,i])
    
        #Assigning negative (water) and positive (land) colormap
        cmap_z=np.vstack((cmap_water,cmap_land))
    
    #'greenearthland' colormap
    #Colors from: http://htmlcolorcodes.com
    elif zcolormap=='greenearthland':
    
        #Mariana Trench depth= 10994 m
        #Mount Everest height= 8848 m
        #(max depth)/(max depth+max height)~=0.55
    
        #Land colormap
        greenearthgreen=np.array([[67, 160, 71],[255, 241, 118]])/255 #Green

        #At least 4 colors requred in colormap to guarantee at least two colors for each of cmap_water and cmap_land colormaps
        greenearthgreen4colors=np.zeros((4,3))
        for i in range(0,3,1):
            greenearthgreen4colors[:,i]=np.interp(np.linspace(1,4,4),np.linspace(1,4,len(greenearthgreen[:,0])),greenearthgreen[:,i])
    
        cmap_z=greenearthgreen4colors.copy()
    
        #Assigning half of the colormap to negative (water) values and other half to positive (land) values 
        M1,N1=np.shape(cmap_z)
        cmap_water=cmap_z[0:int(M1/2),:]
        cmap_land=cmap_z[int(M1/2):,:]
    
    #'blgrtopocmap' colormap
    #Colors from: http://htmlcolorcodes.com
    elif zcolormap=='blgrtopocmap':
    
        #Mariana Trench depth= 10994 m
        #Mount Everest height= 8848 m
        #(max depth)/(max depth+max height)~=0.55
    
        #Water colormap
        blgrtopocmap_water=np.array([[52, 152, 219],[133, 193, 233]])/255
    
        #Land colormap
        blgrtopocmapgreen=np.array([[67, 160, 71],[255, 241, 118]])/255 #Green
    
        #Assigning water and land colormap
        cmap_water1=blgrtopocmap_water.copy()
        cmap_land1=blgrtopocmapgreen.copy()
    
        #Convert to 0.55# sea and 0.45# land 
        #Mariana Trench depth= 10994 m
        #Mount Everest height= 8848 m
        #(max depth)/(max depth+max height)~=0.55
        cmap_water=np.zeros((55,3))
        cmap_land=np.zeros((45,3))
        for i in range(0,3,1):
            cmap_water[:,i]=np.interp(np.linspace(1,55,55),np.linspace(1,55,len(cmap_water1[:,0])),cmap_water1[:,i])
            cmap_land[:,i]=np.interp(np.linspace(1,45,45),np.linspace(1,45,len(cmap_land1[:,0])),cmap_land1[:,i])
    
        #Assigning negative (water) and positive (land) colormap
        cmap_z=np.vstack((cmap_water,cmap_land))
    
    #'blrdtopocmap' colormap
    #Colors from: http://htmlcolorcodes.com
    elif zcolormap=='blrdtopocmap':
    
        #Mariana Trench depth= 10994 m
        #Mount Everest height= 8848 m
        #(max depth)/(max depth+max height)~=0.55
    
        #Water colormap
        blrdtopocmap_water=np.array([[41, 128, 185],[212, 230, 241]])/255 #Blue
    
        #Land colormap
        blrdtopocmap_land=np.array([[242, 215, 213],[192, 57, 43]])/255 #Red
    
        #Assigning water and land colormap
        cmap_water1=blrdtopocmap_water.copy()
        cmap_land1=blrdtopocmap_land.copy()
    
        #Convert to 0.55# sea and 0.45# land 
        #Mariana Trench depth= 10994 m
        #Mount Everest height= 8848 m
        #(max depth)/(max depth+max height)~=0.55
        cmap_water=np.zeros((55,3))
        cmap_land=np.zeros((45,3))
        for i in range(0,3,1):
            cmap_water[:,i]=np.interp(np.linspace(1,55,55),np.linspace(1,55,len(cmap_water1[:,0])),cmap_water1[:,i])
            cmap_land[:,i]=np.interp(np.linspace(1,45,45),np.linspace(1,45,len(cmap_land1[:,0])),cmap_land1[:,i])
    
        #Assigning negative (water) and positive (land) colormap
        cmap_z=np.vstack((cmap_water,cmap_land))
    
    #'grayearth' colormap
    #Colors from: http://htmlcolorcodes.com
    elif zcolormap=='grayearth':
    
        #Mariana Trench depth= 10994 m
        #Mount Everest height= 8848 m
        #(max depth)/(max depth+max height)~=0.55
    
        #Water colormap
        topocmap_water=np.array([[21, 67, 96],[27, 79, 114],[40, 116, 166],[52, 152, 219],[133, 193, 233],[214, 234, 248]])/255
    
        #Land colormap
        grayearthgray=np.array([[179, 182, 183],[66, 73, 73]])/255 #Gray
    
        #Assigning water and land colormap
        cmap_water1=topocmap_water.copy()
        cmap_land1=grayearthgray.copy()
    
        #Convert to 0.55# sea and 0.45# land 
        #Mariana Trench depth= 10994 m
        #Mount Everest height= 8848 m
        #(max depth)/(max depth+max height)~=0.55
        cmap_water=np.zeros((55,3))
        cmap_land=np.zeros((45,3))
        for i in range(0,3,1):
            cmap_water[:,i]=np.interp(np.linspace(1,55,55),np.linspace(1,55,len(cmap_water1[:,0])),cmap_water1[:,i])
            cmap_land[:,i]=np.interp(np.linspace(1,45,45),np.linspace(1,45,len(cmap_land1[:,0])),cmap_land1[:,i])
    
        #Assigning negative (water) and positive (land) colormap
        cmap_z=np.vstack((cmap_water,cmap_land))
    
    #ETOPO1 colormap
    #https://www.ngdc.noaa.gov/mgg/global/global.html
    elif zcolormap=='etopo1':
    
        #http://svn.idldev.com/vis/trunk/src/color/cpt-city/ngdc/ETOPO1.cpt
        # Colormap used in the ETOPO1 global relief map:
        # http://ngdc.noaa.gov/mgg/global/global.html
        #
        # Above sea level is a modified version of GMT_globe.cpt, 
        # Designed by Lester M. Anderson (CASP, UK) lester.anderson@casp.cam.ac.uk,
        # Modified by Jesse Varner and Elliot Lim (NOAA/NGDC) with a smaller band of white for the highest elevations.
        # The ocean colors are adapted from GMT_haxby.cpt, popularized by Bill Haxby, LDEO
        # COLOR_MODEL = RGB
    
        #ETOPO1 elavation range
        #min elevation= -10898.0
        #max elevation= 8271.0
    
        etopo1=np.array([[10,0,121],[26,0,137],[38,0,152],[27,3,166],[16,6,180],[5,9,193],[0,14,203],[0,22,210],[0,30,216],[0,39,223],\
            [12,68,231],[26,102,240],[19,117,244],[14,133,249],[21,158,252],[30,178,255],[43,186,255],[55,193,255],[65,200,255],[79,210,255],\
            [94,223,255],[138,227,255],[188,230,255],[51,102,0],[51,204,102],[187,228,146],[255,220,185],[243,202,137],[230,184,88],[217,166,39],[168,154,31],\
            [164,144,25],[162,134,19],[159,123,13],[156,113,7],[153,102,0],[162,89,89],[178,118,118],[183,147,147],[194,176,176],[204,204,204],\
            [229,229,229],[255,255,255]])/255
    
        etopo1water=np.array([[10,0,121],[26,0,137],[38,0,152],[27,3,166],[16,6,180],[5,9,193],[0,14,203],[0,22,210],[0,30,216],[0,39,223],\
            [12,68,231],[26,102,240],[19,117,244],[14,133,249],[21,158,252],[30,178,255],[43,186,255],[55,193,255],[65,200,255],[79,210,255],\
            [94,223,255],[138,227,255],[188,230,255]])/255
    
        etopo1land=np.array([[51,102,0],[51,204,102],[187,228,146],[255,220,185],[243,202,137],[230,184,88],[217,166,39],[168,154,31],\
            [164,144,25],[162,134,19],[159,123,13],[156,113,7],[153,102,0],[162,89,89],[178,118,118],[183,147,147],[194,176,176],[204,204,204],\
            [229,229,229],[255,255,255]])/255
    
        #Assigning negative (water) and positive (land) colormap
        cmap_water=etopo1water.copy()
        cmap_land=etopo1land.copy()
        cmap_z=etopo1.copy()
    
    #GMT_globe colormap
    #https://www.giss.nasa.gov/tools/panoply/colorbars/
    elif zcolormap=='gmtglobe':
    
        #https://www.giss.nasa.gov/tools/panoply/colorbars/
        #https://www.giss.nasa.gov/tools/panoply/colorbars/gmt/GMT_globe.cpt
        # $Id: GMT_globe.cpt,v 1.1.1.1 2000/12/28 01:23:45 gmt Exp $
        #
        # Colormap using in global relief maps
        # Bathymetry colours manually redefined for blue-shade effect and
        # new topography colour scheme for use with Generic Mapping Tools.
        # Designed by Lester M. Anderson (CASP, UK) lester.anderson@casp.cam.ac.uk
        # COLOR_MODEL = RGB
    
        gmtglobe=np.array([[153,0,255],[153,0,255],[153,0,255],[136,17,255],[119,34,255],[102,51,255],[85,68,255],[68,85,255],[51,102,255],[34,119,255],\
            [17,136,255],[0,153,255],[27,164,255],[54,175,255],[81,186,255],[108,197,255],[134,208,255],[161,219,255],[188,230,255],[215,241,255],\
            [241,252,255],[51,102,0],[51,204,102],[187,228,146],[255,220,185],[243,202,137],[230,184,88],[217,166,39],[168,154,31],[164,144,25],\
            [162,134,19],[159,123,13],[156,113,7],[153,102,0],[162,89,89],[178,118,118],[183,147,147],[194,176,176],[204,204,204],[229,229,229],\
            [242,242,242],[255,255,255]])/255
    
        gmtglobewater=np.array([[153,0,255],[153,0,255],[153,0,255],[136,17,255],[119,34,255],[102,51,255],[85,68,255],[68,85,255],[51,102,255],[34,119,255],\
            [17,136,255],[0,153,255],[27,164,255],[54,175,255],[81,186,255],[108,197,255],[134,208,255],[161,219,255],[188,230,255],[215,241,255],\
            [241,252,255]])/255
    
        gmtglobeland=np.array([[51,102,0],[51,204,102],[187,228,146],[255,220,185],[243,202,137],[230,184,88],[217,166,39],[168,154,31],[164,144,25],\
            [162,134,19],[159,123,13],[156,113,7],[153,102,0],[162,89,89],[178,118,118],[183,147,147],[194,176,176],[204,204,204],[229,229,229],\
            [242,242,242],[255,255,255]])/255
    
        #Assigning negative (water) and positive (land) colormap
        cmap_water=gmtglobewater.copy()
        cmap_land=gmtglobeland.copy()
        cmap_z=gmtglobe.copy()
    
    #GMT_relief colormap
    #https://www.giss.nasa.gov/tools/panoply/colorbars/
    elif zcolormap=='gmtrelief':
    
        #https://www.giss.nasa.gov/tools/panoply/colorbars/
        #https://www.giss.nasa.gov/tools/panoply/colorbars/gmt/GMT_relief.cpt
        #	$Id: GMT_relief.cpt,v 1.1.1.1 2000/12/28 01:23:45 gmt Exp $
        #
        # Colortable for whole earth relief used in Wessel topomaps
        # Designed by P. Wessel and F. Martinez, SOEST
        # COLOR_MODEL = RGB
    
        gmtrelief=np.array([[0,0,0],[0,5,25],[0,10,50],[0,80,125],[0,150,200],[86,197,184],[172,245,168],[211,250,211],[250,255,255],[70,120,50],[120,100,50],[146,126,60],\
            [198,178,80],[250,230,100],[250,234,126],[252,238,152],[252,243,177],[253,249,216],[255,255,255]])/255
    
        gmtreliefwater=np.array([[0,0,0],[0,5,25],[0,10,50],[0,80,125],[0,150,200],[86,197,184],[172,245,168],[211,250,211],[250,255,255]])/255
    
        gmtreliefland=np.array([[70,120,50],[120,100,50],[146,126,60],[198,178,80],[250,230,100],[250,234,126],[252,238,152],[252,243,177],[253,249,216],[255,255,255]])/255
    
        #Assigning negative (water) and positive (land) colormap
        cmap_water=gmtreliefwater.copy()
        cmap_land=gmtreliefland.copy()
        cmap_z=gmtrelief.copy()
    
    #Colormap from  Florian Aendekerk
    #http://www.mathworks.com/matlabcentral/fileexchange/63590-landseacolormap-m-
    elif zcolormap=='aendekerk':
    
        aendekerkcmap=np.array([[0.00000, 0.00000, 0.40000],[0.02524, 0.03810, 0.42857],[0.05048, 0.07619, 0.45714],[0.07571, 0.11429, 0.48571],[0.10095, 0.15238, 0.51429],\
            [0.12619, 0.19048, 0.54286],[0.15143, 0.22857, 0.57143],[0.17667, 0.26667, 0.60000],[0.20190, 0.30476, 0.62857],[0.22714, 0.34286, 0.65714],\
            [0.25238, 0.38095, 0.68571],[0.27762, 0.41905, 0.71429],[0.30286, 0.45714, 0.74286],[0.32810, 0.49524, 0.77143],[0.35333, 0.53333, 0.80000],\
            [0.37857, 0.57143, 0.82857],[0.40381, 0.60952, 0.85714],[0.42905, 0.64762, 0.88571],[0.45429, 0.68571, 0.91429],[0.47952, 0.72381, 0.94286],\
            [0.50476, 0.76190, 0.97143],[0.53000, 0.80000, 1.00000],[0.00000, 0.50000, 0.00000],[1.00000, 1.00000, 0.80000],[1.00000, 0.86667, 0.53333],\
            [1.00000, 0.73333, 0.26667],[1.00000, 0.60000, 0.00000],[1.00000, 0.60000, 0.00000],[0.95000, 0.56250, 0.00000],[0.90000, 0.52500, 0.00000],\
            [0.85000, 0.48750, 0.00000],[0.80000, 0.45000, 0.00000],[0.75000, 0.41250, 0.00000],[0.70000, 0.37500, 0.00000],[0.65000, 0.33750, 0.00000],\
            [0.60000, 0.30000, 0.00000],[0.55000, 0.26250, 0.00000],[0.50000, 0.22500, 0.00000],[0.45000, 0.18750, 0.00000],[0.40000, 0.15000, 0.00000]])
        
        aendekerkcmap_water=np.array([[0.00000, 0.00000, 0.40000],[0.02524, 0.03810, 0.42857],[0.05048, 0.07619, 0.45714],[0.07571, 0.11429, 0.48571],[0.10095, 0.15238, 0.51429],\
            [0.12619, 0.19048, 0.54286],[0.15143, 0.22857, 0.57143],[0.17667, 0.26667, 0.60000],[0.20190, 0.30476, 0.62857],[0.22714, 0.34286, 0.65714],\
            [0.25238, 0.38095, 0.68571],[0.27762, 0.41905, 0.71429],[0.30286, 0.45714, 0.74286],[0.32810, 0.49524, 0.77143],[0.35333, 0.53333, 0.80000],\
            [0.37857, 0.57143, 0.82857],[0.40381, 0.60952, 0.85714],[0.42905, 0.64762, 0.88571],[0.45429, 0.68571, 0.91429],[0.47952, 0.72381, 0.94286],\
            [0.50476, 0.76190, 0.97143],[0.53000, 0.80000, 1.00000]])
        
        aendekerkcmap_land=np.array([[0.00000, 0.50000, 0.00000],[1.00000, 1.00000, 0.80000],[1.00000, 0.86667, 0.53333],\
            [1.00000, 0.73333, 0.26667],[1.00000, 0.60000, 0.00000],[1.00000, 0.60000, 0.00000],[0.95000, 0.56250, 0.00000],[0.90000, 0.52500, 0.00000],\
            [0.85000, 0.48750, 0.00000],[0.80000, 0.45000, 0.00000],[0.75000, 0.41250, 0.00000],[0.70000, 0.37500, 0.00000],[0.65000, 0.33750, 0.00000],\
            [0.60000, 0.30000, 0.00000],[0.55000, 0.26250, 0.00000],[0.50000, 0.22500, 0.00000],[0.45000, 0.18750, 0.00000],[0.40000, 0.15000, 0.00000]])
    
        #Assigning negative (water) and positive (land) colormap
        cmap_water=aendekerkcmap_water.copy()
        cmap_land=aendekerkcmap_land.copy()
        cmap_z=aendekerkcmap.copy()
    
    #User input colormap
    elif ((type(zcolormap) is list) or (type(zcolormap) is np.ndarray)):
    
        cmap_z=zcolormap/255 #Converting RGB to values between 0 and 1
        cmap_z[cmap_z<0]=0
        cmap_z[cmap_z>1]=1
    
        #Assigning half of the colormap to negative (water) values and other half to positive (land) values 
        M1,N1=np.shape(cmap_z)
        cmap_water=cmap_z[0:int(M1/2),:]
        cmap_land=cmap_z[int(M1/2):,:]
    
    #Pre-defined colormap
    else:
    
        pltcmap=mpl.cm.get_cmap(zcolormap)
        cmap_z=[]
        for i in range(pltcmap.N):
            cmap_z.append(pltcmap(i)[0:3])
        cmap_z=np.array(cmap_z).astype('float') #Converting a list to numpy array
    
        #Assigning half of the colormap to negative (water) values and other half to positive (land) values 
        M1,N1=np.shape(cmap_z)
        cmap_water=cmap_z[0:int(M1/2),:]
        cmap_land=cmap_z[int(M1/2):,:]
    
    #--------------------------------------------------------------------------
    #Interpolating a colormap to include ncolor colors

    #Define contour levels
    if waterlandcmapratio=='auto':

        zmingrid=np.nanmin(zgrid)
        zmaxgrid=np.nanmax(zgrid)

        if ((zmingrid<0) and (zmaxgrid>0)):
            num_z_levels=ncolor+1
            z_range=np.max(zgrid)-np.min(zgrid)
            delta_z=(z_range/(num_z_levels))
            z_levels_neg=np.arange(-1,int(np.min(zgrid)/delta_z)-1,-1)*delta_z
            z_levels_pos=np.arange(0,int(np.max(zgrid)/delta_z)+1,1)*delta_z
            z_levels=np.concatenate((z_levels_neg[::-1],z_levels_pos))
        else:
            z_levels=ncolor+1

    #Interpolate colormap
    if waterlandcmapratio=='auto':
    
        #z data have both negative and positive values
        if ((zmingrid<0) and (zmaxgrid>0)):
    
            #Defining negative (water) colormap
            ncolor_neg=len(z_levels_neg)
            cmap_ncolor_neg=np.zeros((ncolor_neg,3))
            for i in range(0,3,1):
                cmap_ncolor_neg[:,i]=np.interp(np.linspace(1,ncolor_neg,ncolor_neg),np.linspace(1,ncolor_neg,len(cmap_water[:,0])),cmap_water[:,i])
    
            #Defining positive (land) colormap
            ncolor_pos=len(z_levels_pos)-1
            cmap_ncolor_pos=np.zeros((ncolor_pos,3))
            for i in range(0,3,1):
                cmap_ncolor_pos[:,i]=np.interp(np.linspace(1,ncolor_pos,ncolor_pos),np.linspace(1,ncolor_pos,len(cmap_land[:,0])),cmap_land[:,i])
    
            #Defining colormap that include ncolor colors
            cmap_ncolor=np.vstack((cmap_ncolor_neg,cmap_ncolor_pos))
    
            #Scale z data
            zgrid[zgrid>np.max(z_levels)]=np.max(z_levels)
            zgrid[zgrid<np.min(z_levels)]=np.min(z_levels)

        #z data onle have positive values
        elif ((zmingrid>=0) and (zmaxgrid>=0)):
        
            ncolor_neg=0
            ncolor_pos=ncolor
    
            #Defining land colormap
            cmap_ncolor_pos=np.zeros((ncolor_pos,3))
            for i in range(0,3,1):
                cmap_ncolor_pos[:,i]=np.interp(np.linspace(1,ncolor_pos,ncolor_pos),np.linspace(1,ncolor_pos,len(cmap_land[:,0])),cmap_land[:,i])
    
            #Defining colormap that include ncolor colors
            cmap_ncolor=cmap_ncolor_pos.copy()
    
        #z data only have negative values
        elif ((zmingrid<=0) and (zmaxgrid<=0)):
    
            ncolor_neg=ncolor
            ncolor_pos=0
    
            #Defining water colormap
            cmap_ncolor_neg=np.zeros((ncolor_neg,3))
            for i in range(0,3,1):
                cmap_ncolor_neg[:,i]=np.interp(np.linspace(1,ncolor_neg,ncolor_neg),np.linspace(1,ncolor_neg,len(cmap_water[:,0])),cmap_water[:,i])
    
            #Defining colormap that include ncolor colors
            cmap_ncolor=cmap_ncolor_neg.copy()
    
    
    #Colormap interpolated to include ncolor colors without considering water_land ratio
    elif waterlandcmapratio=='none':
    
        cmap_ncolor=np.zeros((ncolor,3))
        for i in range(0,3,1):
            cmap_ncolor[:,i]=np.interp(np.linspace(1,ncolor,ncolor),np.linspace(1,ncolor,len(cmap_z[:,0])),cmap_z[:,i])
    
    #Assigning waterlandcmapratio of the colormap to negative (water) values and the rest to positive (land) values 
    elif ((isinstance(waterlandcmapratio,int)==True) or (isinstance(waterlandcmapratio,float)==True)):
    
        #Checking waterlandcmapratio value
        if ((waterlandcmapratio<=0) or (waterlandcmapratio>1)): waterlandcmapratio=0.5
    
        #Assigning water_land ratio to to negative (water) and positive (land) colormaps
        #M1,N1=np.shape(cmap_z)
        #cmap_waterscaled=cmap_z[0:int(M1/2),:]
        #cmap_landscaled=cmap_z[int(M1/2):,:]
    
        #Defining water colormap
        ncolor_neg=(int(ncolor*waterlandcmapratio))
        cmap_ncolor_neg=np.zeros((ncolor_neg,3))
        for i in range(0,3,1):
            cmap_ncolor_neg[:,i]=np.interp(np.linspace(1,ncolor_neg,ncolor_neg),np.linspace(1,ncolor_neg,len(cmap_water[:,0])),cmap_water[:,i])
    
        #Defining land colormap
        ncolor_pos=ncolor-ncolor_neg
        #ncolor_pos=ncolor-(fix(M1/2)*waterlandcmapratio)
        cmap_ncolor_pos=np.zeros((ncolor_pos,3))
        for i in range(0,3,1):
            cmap_ncolor_pos[:,i]=np.interp(np.linspace(1,ncolor_pos,ncolor_pos),np.linspace(1,ncolor_pos,len(cmap_land[:,0])),cmap_land[:,i])
    
        #Defining colormap that include ncolor colors
        cmap_ncolor=np.vstack((cmap_ncolor_neg,cmap_ncolor_pos))
    
    #Generating matplotlib colormap 
    #Use matplotlib.colors.ListedColormap or matplotlib.colors.LinearSegmentedColormap
    #cmap_ncolor_mpl=colors.LinearSegmentedColormap.from_list(name='my_colormap',colors=cmap_ncolor,N=ncolor)
    cmap_ncolor_mpl=colors.ListedColormap(cmap_ncolor,name='my_colormap',N=None)
    
    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='imagesc': #2 dimensional plot using imagesc or imshow
        
        plt.imshow(zgrid,extent=[np.nanmin(xgrid),np.nanmax(xgrid),np.nanmin(ygrid),np.nanmax(ygrid)],interpolation=None,cmap=cmap_ncolor_mpl,aspect='auto',origin='lower')
    
    elif dispout=='pcolor': #2 dimensional plot using pcolor
    
        plt.pcolormesh(xgrid,ygrid,zgrid,cmap=cmap_ncolor_mpl,shading='nearest')

    elif dispout=='contour': #2 dimensional contour plot
    
        zmingrid=np.nanmin(zgrid)
        zmaxgrid=np.nanmax(zgrid)
        plt.contourf(xgrid,ygrid,zgrid,z_levels,cmap=cmap_ncolor_mpl)
    
    elif dispout=='surface': #Surface plot
        
        fig=plt.figure()
        fig3d=fig.gca(projection='3d')
        #fig3d.view_init(elev=50,azim=-20)
        #surf1=fig3d.plot_surface(xgrid,ygrid,zgrid,cmap=cmap_ncolor_mpl,edgecolor='none')
        surf1=fig3d.plot_surface(xgrid,ygrid,zgrid,cmap=cmap_ncolor_mpl,rstride=1,cstride=1)
    
    
    #Setting plot properties
    if ((dispout=='imagesc') or (dispout=='pcolor') or (dispout=='contour')):
    
        #Setting axis limits
        plt.xlim(np.nanmin(xgrid),np.nanmax(xgrid))
        plt.ylim(np.nanmin(ygrid),np.nanmax(ygrid))
    
        #Plotting colorbar
        if dispcolorbar=='yes':
            cbar=plt.colorbar()
            cbar.set_label(cbarlabel,fontsize=axisfontsize)
            #cbar.ax.set_title(cbarlabel,fontsize=axisfontsize)
            cbar.ax.tick_params(labelsize=axisfontsize)
            #c=np.max(np.abs([zmin,zmax]))
    
        #Plotting grid lines
        if dispgrid=='yes':
            plt.grid()
            plt.xticks(np.linspace(np.nanmin(xgrid),np.nanmax(xgrid),13))
            plt.yticks(np.linspace(np.nanmin(ygrid),np.nanmax(ygrid),13))
        elif ((isinstance(dispgrid,int)==True) or (isinstance(dispgrid,float)==True)):
            plt.grid()
            plt.xticks(np.linspace(np.nanmin(xgrid),np.nanmax(xgrid),dispgrid+1))
            plt.yticks(np.linspace(np.nanmin(ygrid),np.nanmax(ygrid),dispgrid+1))
    
        #Setting figure font size
        plt.rcParams.update({'font.size':plotfontsize})
    
        #Setting axis label
        plt.xlabel(xaxislabel,fontsize=axisfontsize)
        plt.ylabel(yaxislabel,fontsize=axisfontsize)
        

    elif dispout=='surface':

        #Setting axis limits
        fig3d.set_xlim(np.nanmin(xgrid),np.nanmax(xgrid))
        fig3d.set_ylim(np.nanmin(ygrid),np.nanmax(ygrid))
    
        #Plotting colorbar
        if dispcolorbar=='yes':
            cbar=fig.colorbar(surf1)
            cbar.set_label(cbarlabel,fontsize=axisfontsize)
            #cbar.ax.set_title(cbarlabel,fontsize=axisfontsize)
            cbar.ax.tick_params(labelsize=axisfontsize)
            #c=np.max(np.abs([zmin,zmax]))
    
        #Plotting grid lines
        if dispgrid=='yes':
            fig.grid()
            fig3d.set_xticks(np.linspace(np.nanmin(xgrid),np.nanmax(xgrid),13))
            fig3d.set_yticks(np.linspace(np.nanmin(ygrid),np.nanmax(ygrid),13))
        elif ((isinstance(dispgrid,int)==True) or (isinstance(dispgrid,float)==True)):
            fig.grid()
            fig3d.set_xticks(np.linspace(np.nanmin(xgrid),np.nanmax(xgrid),dispgrid+1))
            fig3d.set_yticks(np.linspace(np.nanmin(ygrid),np.nanmax(ygrid),dispgrid+1))
    
        #Setting figure font size
        mpl.rcParams.update({'font.size':plotfontsize})
    
        #Setting axis label
        fig3d.set_xlabel(xaxislabel,fontsize=axisfontsize)
        fig3d.set_ylabel(yaxislabel,fontsize=axisfontsize)

    #--------------------------------------------------------------------------
    #Outputs
    return xgrid, ygrid, zgrid, cmap_ncolor, cmap_ncolor_mpl

    #--------------------------------------------------------------------------
