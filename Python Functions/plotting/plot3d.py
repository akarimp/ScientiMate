def plot3d(x, y, z, plottype='imagesc', cmapcolors='blue'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2019-02-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.plot3d
    ==================

    .. code:: python

        scientimate.plot3d(x, y, z, plottype='imagesc', cmapcolors='blue')

    Description
    -----------

    Plot x , y, and z data in 2-d/3-d contour/surface plot

    Inputs
    ------

    x
        | x data
        | Set x=[] if it is not available
        | It may be 1d or 2d array
    y
        | y data
        | Set y=[] if it is not available
        | It may be 1d or 2d array
    z
        | z data
        | It may be 1d or 2d array
    plottype='imagesc'
        | Plot type
        | 'imagesc': 2 dimensional plot using imagesc or imshow
        | 'pcolor': 2 dimensional plot using pcolor
        | 'contour': 2 dimensional contour plot, number of contour=32
        | 'contourf': 2 dimensional filled contour plot, number of contour=32
        | 'surface': 3 dimensional surface plot 
    cmapcolors='blue'
        | Colormap style
        | 'blue': blue colormap
        | 'red': red colormap
        | 'green': green colormap
        | 'yellow': yellow colormap
        | 'purple': purple colormap
        | 'brown': brown colormap
        | 'gray': gray colormap
        | 'blue_red': blue-red colormap
        | 'red_blue': red-blue colormap
        | 'blue_green': blue-green colormap
        | 'green_blue': green-blue colormap
        | 'green_yellow': green-yellow colormap
        | 'yellow_green': yellow-green colormap
        | 'red_yellow': red-yellow colormap
        | 'yellow_red': yellow-red colormap
        | 'cyclic': cyclic/oscillation colormap 
        | 'seq': sequential colormap
        | User-defined colors may be used to generate colormap
        | User-defined colors should be defined as (M,3) array in RGB color format
        | At least two colors should be defined, i.e. M should be equal or larger than 2
        | User-defined colors values should be between 0 and 255
        | Any available colormap name such as 'cool', 'winter', etc also can be used

    Outputs
    -------


    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        x,y=np.meshgrid(np.linspace(-10,10,50),np.linspace(-10,10,50))
        r=np.sqrt(x**2+y**2)+1e-10 #Add 1e-10 to prevent divide by 0
        z=np.sin(r)/r
        sm.plot3d(x,y,z,'surface','purple')

        x,y=np.meshgrid(np.linspace(-10,10,21),np.linspace(-10,10,21))
        z=(np.sin(x)+np.sin(y))/(x+y+1e-10) #Add 1e-10 to prevent divide by 0
        sm.plot3d(x,y,z,'pcolor','purple')

        x=10*np.random.rand(1000)
        y=10*np.random.rand(1000)
        z=x**2+y**2
        sm.plot3d(x,y,z,'pcolor','blue')

    References
    ----------

    Colormap

    * https://matplotlib.org/tutorials/colors/colormaps.html
    * https://www.mathworks.com/help/PYTHON/ref/colormap.html
    * http://colorbrewer2.org
    * http://matplotlib.org/cmocean/
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

    #--------------------------------------------------------------------------
    #Assigning default values

    if x is None: x=np.array([]) #Check if x is empty
    if y is None: y=np.array([]) #Check if y is empty

    #--------------------------------------------------------------------------
    #Set plot drawing size

    plotfontsize=18 #Size of plot fonts
    axisfontsize=18 #Size of axis fonts


    #--------------------------------------------------------------------------
    #Input variable names

    #At the moment there in no way to find an exaxt name of input arguments that passed to function
    #Python insert module provide functions to chack some input properties
    xaxislabel='x' #x axis label
    yaxislabel='y' #y axis label
    zaxislabel='z' #z axis label
    cbarlabel='z' #Colorbar label


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

    #Defining domain boundaries
    xmin=np.nanmin(x)
    xmax=np.nanmax(x)
    ymin=np.nanmin(y)
    ymax=np.nanmax(y)
    zmin=np.nanmin(z)
    zmax=np.nanmax(z)


    #--------------------------------------------------------------------------

    ncolor=32 #Number of colors to be used in colormap
    gridsize=100 #Grid size in x and y directions to interpolate elevation data on them
    interpMethod='nearest' #Interpolation method ('linear','nearest')

    #--------------------------------------------------------------------------
    #Interpolating z data

    #Check if inputs are 1-d or 2-d arrays
    if (Mx==1 or Nx==1):

        #Calculating grid size in each direction using number grid points 
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
    #Colors

    light_blue=np.array([212, 230, 241])
    mid_blue=np.array([41, 128, 185])
    dark_blue=np.array([21, 67, 96])

    light_red=np.array([242, 215, 213])
    mid_red=np.array([192, 57, 43])
    dark_red=np.array([100, 30, 22])

    light_green=np.array([200, 230, 201])
    mid_green=np.array([76, 175, 80])
    dark_green=np.array([27, 94, 32])

    light_yellow=np.array([255, 249, 196])
    mid_yellow=np.array([255, 235, 59])
    dark_yellow=np.array([245, 127, 23])

    light_purple=np.array([235, 222, 240])
    mid_purple=np.array([155, 89, 182])
    dark_purple=np.array([81, 46, 95])

    light_brown=np.array([246, 221, 204])
    mid_brown=np.array([211, 84, 0])
    dark_brown=np.array([110, 44, 0])

    light_gray=np.array([250, 250, 250])
    dark_gray=np.array([1, 1, 1])


    #--------------------------------------------------------------------------
    #Colormaps

    #zclmap
    #'blue': blue colormap
    if cmapcolors=='blue':
        zclmap=np.vstack((light_blue,mid_blue,dark_blue))/255

    #'red': red colormap
    elif cmapcolors=='red':
        zclmap=np.vstack((light_red,mid_red,dark_red))/255

    #'green': green colormap
    elif cmapcolors=='green':
        zclmap=np.vstack((light_green,mid_green,dark_green))/255

    #'yellow': yellow colormap
    elif cmapcolors=='yellow':
        zclmap=np.vstack((light_yellow,mid_yellow,dark_yellow))/255

    #'purple': purple colormap
    elif cmapcolors=='purple':
        zclmap=np.vstack((light_purple,mid_purple,dark_purple))/255

    #'brown': brown colormap, 'cooper'
    elif cmapcolors=='brown':
        zclmap=np.vstack((light_brown,mid_brown,dark_brown))/255

    #'gray': gray colormap, 'gray'
    elif cmapcolors=='gray':
        zclmap=np.vstack((light_gray,dark_gray))/255

    #'blue_red': blue-red colormap, 'cool'
    elif cmapcolors=='blue_red':
        zclmap=np.vstack((mid_blue,mid_red))/255

    #'red_blue': red-blue colormap, 'cool-1'
    elif cmapcolors=='red_blue':
        zclmap=np.vstack((mid_red,mid_blue))/255

    #'blue_green': blue-green colormap, 'winter'
    elif cmapcolors=='blue_green':
        zclmap=np.vstack((mid_blue,mid_green))/255

    #'green_blue': green-blue colormap, 'winter-1'
    elif cmapcolors=='green_blue':
        zclmap=np.vstack((mid_green,mid_blue))/255

    #'green_yellow': green-yellow colormap, 'summer'
    elif cmapcolors=='green_yellow':
        zclmap=np.vstack((mid_green,mid_yellow))/255

    #'yellow_green': yellow-green colormap, 'summer-1'
    elif cmapcolors=='yellow_green':
        zclmap=np.vstack((mid_yellow,mid_green))/255

    #'red_yellow': red-yellow colormap, 'autumn'
    elif cmapcolors=='red_yellow':
        zclmap=np.vstack((mid_red,mid_yellow))/255

    #'yellow_red': yellow-red colormap, 'autumn-1'
    elif cmapcolors=='yellow_red':
        zclmap=np.vstack((mid_yellow,mid_red))/255

    #'cyclic': cyclic/oscillation colormap, 'hsv'
    elif cmapcolors=='cyclic':
        #colormap('hsv')
        #zclmap=np.vstack((mid_red,light_red,light_blue,mid_blue,mid_blue,light_blue,light_red,mid_red))/255
        zclmap=np.vstack((light_red,light_blue,mid_blue,mid_red))/255

    #'seq': sequential colormap, 'line'
    elif cmapcolors=='seq':
        zclmap=mpl.cm.get_cmap()

    #User input colormap
    elif ((type(cmapcolors) is list) or (type(cmapcolors) is np.ndarray)):
        zclmap=cmapcolors/255 #Converting RGB to values between 0 and 1
        zclmap[zclmap<0]=0
        zclmap[zclmap>1]=1

    #Pre-defined colormap such as 'cool', 'winter', etc
    else:
        pltcmap=mpl.cm.get_cmap(cmapcolors)
        zclmap=[]
        for i in range(pltcmap.N):
            zclmap.append(pltcmap(i)[0:3])
        zclmap=np.array(zclmap).astype('float') #Converting a list to numpy array


    #--------------------------------------------------------------------------
    #Interpolating a colormap to include ncolor zclmap

    if cmapcolors=='seq':
        #Line colors cycle in Hex
        line_color_cycle_hex=mpl.rcParams['axes.prop_cycle'].by_key()['color']
        
        #Convert Hex colors to RGB
        line_color_cycle_rgb=[]
        for c in line_color_cycle_hex:
            rgb_tuple=mpl.colors.to_rgb(c)
            rgb_list=list(rgb_tuple)
            line_color_cycle_rgb.append(rgb_list)
        
        line_color_cycle_rgb=np.array(line_color_cycle_rgb) #Convert to numpy array
        cmap_ncolor=np.tile(line_color_cycle_rgb,(int(ncolor/10)+1,1)) #Repeat enough colors for Nxdata data series 

    else:
        cmap_ncolor=np.zeros((ncolor,3))
        for i in range(0,3,1):
            cmap_ncolor[:,i]=np.interp(np.linspace(1,ncolor,ncolor),np.linspace(1,ncolor,len(zclmap[:,0])),zclmap[:,i])

        #Generating matplotlib colormap 
        #Use matplotlib.colors.ListedColormap or matplotlib.colors.LinearSegmentedColormap
        #cmap_ncolor_plt=colors.LinearSegmentedColormap.from_list(name='my_colormap',colors=cmap_ncolor,N=ncolor)
        cmap_ncolor_plt=colors.ListedColormap(cmap_ncolor,name='my_colormap',N=None)


    #--------------------------------------------------------------------------
    #2 dimensional plot using imagesc or imshow

    if ((plottype=='imagesc') or (plottype=='imagesc_grid')):
        
        plt.imshow(zgrid,extent=[np.nanmin(xgrid),np.nanmax(xgrid),np.nanmin(ygrid),np.nanmax(ygrid)],interpolation=None,cmap=cmap_ncolor_plt,aspect='auto',origin='lower')


    #--------------------------------------------------------------------------
    #2 dimensional plot using pcolor

    if ((plottype=='pcolor') or (plottype=='pcolor_grid')):

        plt.pcolormesh(xgrid,ygrid,zgrid,cmap=cmap_ncolor_plt,shading='nearest')


    #--------------------------------------------------------------------------
    #2 dimensional contour plot

    if ((plottype=='contour') or (plottype=='contour_grid')):

        zmingrid=np.nanmin(zgrid)
        zmaxgrid=np.nanmax(zgrid)
        zlevels=np.linspace(zmingrid,zmaxgrid,ncolor) #Deviding z levels between zmin and zmax to ncolor
        plt.contour(xgrid,ygrid,zgrid,zlevels,cmap=cmap_ncolor_plt)


    #--------------------------------------------------------------------------
    #2 dimensional contour plot

    if ((plottype=='contourf') or (plottype=='contourf_grid')):

        zmingrid=np.nanmin(zgrid)
        zmaxgrid=np.nanmax(zgrid)
        zlevels=np.linspace(zmingrid,zmaxgrid,ncolor) #Deviding z levels between zmin and zmax to ncolor
        plt.contourf(xgrid,ygrid,zgrid,zlevels,cmap=cmap_ncolor_plt)


    #--------------------------------------------------------------------------
    #Surface plot

    if plottype=='surface':
        
        fig=plt.figure()
        fig3d=fig.gca(projection='3d')
        #fig3d.view_init(elev=50,azim=-20)
        #surf1=fig3d.plot_surface(xgrid,ygrid,zgrid,cmap=cmap_ncolor_plt,edgecolor='none')
        surf1=fig3d.plot_surface(xgrid,ygrid,zgrid,cmap=cmap_ncolor_plt,rstride=1,cstride=1)


    #--------------------------------------------------------------------------
    #Set plot properties

    if ((plottype=='imagesc') or (plottype=='pcolor') or (plottype=='contour') or (plottype=='contourf')):

        #Set axis limits
        plt.xlim((np.nanmin(xgrid),np.nanmax(xgrid)))
        plt.ylim((np.nanmin(ygrid),np.nanmax(ygrid)))

        #Plot colorbar
        cbar=plt.colorbar()
        cbar.set_label(cbarlabel,fontsize=axisfontsize)
        #cbar.ax.set_title(cbarlabel,fontsize=axisfontsize)
        cbar.ax.tick_params(labelsize=axisfontsize)
        #c=np.max(np.abs([zmin,zmax]))

        #Plot grid lines
        if ((plottype=='imagesc_grid') or (plottype=='pcolor_grid') \
            or (plottype=='contour_grid') or (plottype=='contourf_grid')):
            plt.grid()

        #Set figure font size
        plt.rcParams.update({'font.size':plotfontsize})

        #Set axis label
        plt.xlabel(xaxislabel,fontsize=axisfontsize)
        plt.ylabel(yaxislabel,fontsize=axisfontsize)
    

    elif plottype=='surface':

        #Set axis limits
        fig3d.set_xlim(np.nanmin(xgrid),np.nanmax(xgrid))
        fig3d.set_ylim(np.nanmin(ygrid),np.nanmax(ygrid))
        fig3d.set_zlim(np.nanmin(zgrid),np.nanmax(zgrid))
    
        #Plot colorbar
        cbar=fig.colorbar(surf1)
        cbar.set_label(cbarlabel,fontsize=axisfontsize)
        #cbar.ax.set_title(cbarlabel,fontsize=axisfontsize)
        cbar.ax.tick_params(labelsize=axisfontsize)
        #c=np.max(np.abs([zmin,zmax]))
    
        #Set figure font size
        mpl.rcParams.update({'font.size':plotfontsize})
    
        #Set axis label
        fig3d.set_xlabel(xaxislabel,fontsize=axisfontsize)
        fig3d.set_ylabel(yaxislabel,fontsize=axisfontsize)

    #--------------------------------------------------------------------------
