def plot3dhillshades(x, y, z, plottype):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2020-02-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.plot3dhillshades
    ============================

    .. code:: python

        scientimate.plot3dhillshades(x, y, z, plottype)

    Description
    -----------

    Plot hillshades (shaded relief) of x (longitude), y (latitude) and z (elevation) data

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
    plottype='imagesc'
        | Plot type
        | 'imagesc': 2 dimensional plot using imagesc or imshow
        | 'pcolor': 2 dimensional plot using pcolor
        | 'contour': 2 dimensional contour plot, number of contour=ncolor
        | 'surface': 3 dimensional surface plot 

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
        sm.plot3dhillshades(x,y,z,'imagesc')

        x,y=np.meshgrid(np.linspace(-10,10,21),np.linspace(-10,10,21))
        z=(np.sin(x)+np.sin(y))/(x+y+1e-10) #Add 1e-10 to prevent divide by 0
        sm.plot3dhillshades(x,y,z,'imagesc')

        x=10*np.random.rand(1000)
        y=10*np.random.rand(1000)
        z=x**2+y**2
        sm.plot3dhillshades(x,y,z,'imagesc')


    References
    ----------

    * https://pro.arcgis.com/en/pro-app/tool-reference/3d-analyst/how-hillshade-works.htm
    * https://pro.arcgis.com/en/pro-app/tool-reference/3d-analyst/applying-a-z-factor.htm
    * http://mike.teczno.com/img/hillshade.py
    * http://geospatialpython.com/2013/12/python-and-elevation-data-creating.html
    * https://github.com/ThomasLecocq/geophysique.be/blob/master/2014-02-25#20Shaded#20Relief#20Map#20in#20Python.ipynb
    * https://www.geophysique.be/2014/02/25/shaded-relief-map-in-python/
    * https://www.mathworks.com/matlabcentral/fileexchange/14863-hillshade
    * http://www.reliefshading.com/

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

    ncolor=32 #Number of colors to be used in colormap
    plotfontsize=18 #Size of plot fonts
    axisfontsize=18 #Size of axis fonts
    xaxislabel='x' #x axis label
    yaxislabel='y' #y axis label
    cbarlabel='z' #Colorbar label
    dispcolorbar='no' #Define to display colorbar or not ('yes': display, 'no': not display)
    dispgrid='no' #Define to display grid lines or not ('yes': display, 'no': not display)
    gridsize=100 #Grid size in x (longitude) and y (latitude) directions to interpolate elevation data on them
    gridsizetype='points' #Grid size type 
    xmin=np.nanmin(x) #Minimum x (longitude) of domain to be plotted
    xmax=np.nanmax(x) #Maximum x (longitude) of domain to be plotted
    ymin=np.nanmin(y) #Minimum y (latitude) of domain to be plotted
    ymax=np.nanmax(y) #Maximum y (latitude) of domain to be plotted
    zmin=np.nanmin(z) #Minimum z (elevation) of domain to be plotted
    zmax=np.nanmax(z) #Maximum z (elevation) of domain to be plotted
    interpMethod='nearest' #Interpolation method 

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
    #Calculate Hillshade using ArcGIS method

    #Cell size
    dx=xgrid[0,1]-xgrid[0,0]
    dy=ygrid[1,0]-ygrid[0,0]

    #Altitude and Azimuth of illumination source
    altitude_deg=45 #Angle of illumination source above the horizon (default 45)
    azimuth_deg=315 #Angular direction of the sun, measured from north in clockwise from 0 to 360 (default 315), Exp: azimuth of 90 degrees is east

    #Calculate z-factor
    #z-factor is a conversion factor that adjusts units of measure for vertical (elevation) units 
    #when they are different from horizontal coordinate (x,y) units of the input surface
    #If (x,y) units are decimal degrees and z-units are meters, then z-factor is:
    #known_latitude=[0,10,20,30,40,50,60,70,80]
    #known_z_factor=[0.00000898,0.00000912,0.00000956,0.00001036,0.00001171,0.00001395,0.00001792,0.00002619,0.00005156] #z-factor if z is in meter
    #z_factor=interp1(known_latitude,known_z_factor,ygrid,'linear',1)
    z_factor=1

    #Compute illumination angle
    zenith_deg=90-altitude_deg #Angle from zenith point (directly overhead) to the direction of illumination source
    zenith_rad=np.deg2rad(zenith_deg) #Convert to radian

    #Compute illumination direction
    azimuth_math=360-azimuth_deg+90 #Convert from geographic unit (compass direction) to a mathematic unit (right angle)
    if azimuth_math>360 : azimuth_math=azimuth_math-360 #Make sure azimuth_math is equal or less than 360
    azimuth_rad=np.deg2rad(azimuth_math) #Convert to radian

    #Compute slope and aspect
    dz_dy,dz_dx=np.gradient(zgrid,dx,dy) #Rate of change in the x and y directions
    slope_rad=np.arctan(z_factor*np.sqrt(dz_dx**2+dz_dy**2)) #Steepest downhill descent from each cell in the surface
    aspect_rad=np.arctan2(dz_dy,-dz_dx) #Direction that steepest downslope direction is facing
    aspect_rad[aspect_rad<0]=2*np.pi+aspect_rad[aspect_rad<0] #Make sure aspect_rad is equal or larger than 0
    aspect_rad[(dz_dx==0) & (dz_dy>0)]=np.pi/2 #Check for dz_dx==0 and dz_dy>0
    aspect_rad[(dz_dx==0) & (dz_dy<0)]=2*np.pi-np.pi/2 #Check for dz_dx==0 and dz_dy<0

    #Calculate hillshade
    hillshade=255*((np.cos(zenith_rad)*np.cos(slope_rad))+(np.sin(zenith_rad)*np.sin(slope_rad)*np.cos(azimuth_rad-aspect_rad)))
    hillshade[hillshade<0]=0 #Check for hillshade<0

    #-----------------------------------------------------------------------------
    #Example

    #import numpy as np
    #zgrid=[[2450,2461,2483],[2452,2461,2483],[2447,2455,2477]]
    #dx=5
    #dy=5
    #altitude_deg=45
    #azimuth_deg=315
    #z_factor=1
    #zenith_deg=90-altitude_deg #Answer: 45
    #zenith_rad=np.deg2rad(zenith_deg) #Answer: 0.7857142857
    #azimuth_math=360-azimuth_deg+90 #Answer: 135
    #if azimuth_math>360 : azimuth_math=azimuth_math-360 #Answer: 135
    #azimuth_rad=np.deg2rad(azimuth_math) #Answer: 2.3571428571
    #dz_dy,dz_dx=np.gradient(zgrid,dx,dy) #Answer for center point: dz_dx=3.125, dz_dy=-0.525
    #slope_rad=np.arctan(z_factor*np.sqrt(dz_dx**2+dz_dy**2)) #Answer for center point: 1.26511
    #aspect_rad=np.arctan2(dz_dy,-dz_dx) #Answer for center point: -2.9751469600412
    #aspect_rad[aspect_rad<0]=2*np.pi+aspect_rad[aspect_rad<0] #Answer for center point: 3.310567
    #aspect_rad[(dz_dx==0) & (dz_dy>0)]=np.pi/2 #Answer for center point: 3.310567
    #aspect_rad[(dz_dx==0) & (dz_dy<0)]=2*np.pi-np.pi/2 #Answer for center point: 3.310567
    #hillshade=255*((np.cos(zenith_rad)*np.cos(slope_rad))+(np.sin(zenith_rad)*np.sin(slope_rad)*np.cos(azimuth_rad-aspect_rad))) #Answer for center point: 153.82
    #hillshade[hillshade<0]=0 #Answer for center point: 153.82

    #--------------------------------------------------------------------------
    #Displaying results

    if plottype=='imagesc': #2 dimensional plot using imagesc or imshow
        
        plt.imshow(hillshade,extent=[np.nanmin(xgrid),np.nanmax(xgrid),np.nanmin(ygrid),np.nanmax(ygrid)],interpolation=None,cmap=mpl.cm.get_cmap('gray'),aspect='auto',origin='lower')
    
    elif plottype=='pcolor': #2 dimensional plot using pcolor
    
        plt.pcolormesh(xgrid,ygrid,hillshade,cmap=mpl.cm.get_cmap('gray'),shading='nearest')

    elif plottype=='contour': #2 dimensional contour plot
    
        zmingrid=np.nanmin(hillshade)
        zmaxgrid=np.nanmax(hillshade)
        zlevels=np.linspace(zmingrid,zmaxgrid,ncolor) #Deviding z levels between zmin and zmax to ncolor
        plt.contourf(xgrid,ygrid,hillshade,zlevels,cmap=mpl.cm.get_cmap('gray'))
    
    elif plottype=='surface': #Surface plot
        
        fig=plt.figure()
        fig3d=fig.gca(projection='3d')
        #fig3d.view_init(elev=50,azim=-20)
        #surf1=fig3d.plot_surface(xgrid,ygrid,hillshade,cmap=mpl.cm.get_cmap('gray'),edgecolor='none')
        surf1=fig3d.plot_surface(xgrid,ygrid,hillshade,cmap=mpl.cm.get_cmap('gray'),rstride=1,cstride=1)
    
    
    #Setting plot properties
    if ((plottype=='imagesc') or (plottype=='pcolor') or (plottype=='contour')):
    
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
        

    elif plottype=='surface':

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
    #return xgrid, ygrid, hillshade

    #--------------------------------------------------------------------------
