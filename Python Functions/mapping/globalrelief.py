def globalrelief(filepath, xmin=-180, xmax=180, ymin=-90, ymax=90, dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2020-09-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.globalrelief
    ========================

    .. code:: python

        x, y, z, xyz, xgrid, ygrid, zgrid = scientimate.globalrelief(filepath, xmin=-180, xmax=180, ymin=-90, ymax=90, dispout='no')

    Description
    -----------

    | Return x (longitude), y (latitude) and z (elevation) data from ETOPO1 Global Relief Model (Amante & Eakins, 2009) interpolated on 0.08 degree grid
    | ETOPO1 is 1 arc-minute global relief, however, this data are interpolated on 4.8 arc-minute (0.08 degree)
    | This data are obtained from ETOPO1 Global Relief Bedrock (grid-registered)
    | ETOPO1 horizontal datum: WGS 84 geographic
    | ETOPO1 vertical datum: sea level
    | https://www.ngdc.noaa.gov/mgg/global/global.html

    Inputs
    ------

    filepath
        | Path of the folder that contains 'ETOPO1_XYZ_0_08_Deg_Python.npz' file
        | Download ETOPO1_XYZ_0_08_Deg_Python.npz file from
        | https://github.com/akarimp/ScientiMate/releases/download/2.0/ETOPO1_XYZ_0_08_Deg_Python.npz
        | Example: filepath = r'C:'
    xmin=-180
        | Minimum x (longitude) of domain to be returned in degree 
        | It should be between -180 and 180 degree
    xmax=180
        | Maximum x (longitude) of domain to be returned in degree
        | It should be between -180 and 180 degree
    ymin=-90
        | Minimum y (latitude) of domain to be returned in degree
        | It should be between -90 and 90 degree
    ymax=90
        | Maximum y (latitude) of domain to be returned in degree
        | It should be between -90 and 90 degree
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    x
        Interpolated x (longitude) data in degree
    y
        Interpolated y (latitude) data in degree
    z
        Interpolated z (elevation) data in degree
    xyz
        | Interpolated xyz data
        | xyz is a 3-column array
        | 1st column contains longitude (x) data in (Degree)
        | 2nd column contains latitude (y) data in (Degree)
        | 3rd column contains elevation (z) data in (m)
    xgrid
        Interpolated x (longitude) data on 2d mesh in degree
    ygrid
        Interpolated y (latitude) data on 2d mesh in degree
    zgrid
        Interpolated z (elevation) data on 2d mesh in degree

    Examples
    --------

    .. code:: python

        import scientimate as sm

        #Path of ETOPO1_XYZ_0_08_Deg_Python.npz
        filepath = r'C:'

        #Globe
        x,y,z,xyz,xgrid,ygrid,zgrid=sm.globalrelief(filepath,-180,180,-90,90,'yes')

        #Middle East
        x,y,z,xyz,xgrid,ygrid,zgrid=sm.globalrelief(filepath,24,64,9,43,'yes')

        #North America
        x,y,z,xyz,xgrid,ygrid,zgrid=sm.globalrelief(filepath,-169,-8,5,90,'yes')

    References
    ----------

    ETOPO1 Global Relief Model

        Amante, C. and B.W. Eakins, 2009.
        ETOPO1 1 Arc-Minute Global Relief Model: Procedures, Data Sources and Analysis.
        NOAA Technical Memorandum NESDIS NGDC-24. National Geophysical Data Center, NOAA.
        doi:10.7289/V5C8276M

    * https://www.ngdc.noaa.gov/mgg/global/global.html
    * https://maps.ngdc.noaa.gov/viewers/grid-extract/index.html
    * https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/
    * https://data.nodc.noaa.gov/cgi-bin/iso?id=gov.noaa.ngdc.mgg.dem:316

    GEBCO Global ocean & land terrain models

    * https://www.gebco.net/data_and_products/gridded_bathymetry_data/

    Natural Earth 1:10m Raster Data

    * https://www.naturalearthdata.com/downloads/10m-raster-data/

    Geospatial data

    * https://www.mathworks.com/help/map/finding-geospatial-data.html
    * https://www.ngdc.noaa.gov/mgg/global/etopo2.html
    * https://www.ngdc.noaa.gov/mgg/global/etopo5.HTML
    * https://www.ngdc.noaa.gov/mgg/image/2minrelief.html
    * https://www.ngdc.noaa.gov/mgg/coastal/crm.html
    * https://www.ngdc.noaa.gov/mgg
    * https://viewer.nationalmap.gov/launch
    * https://earthexplorer.usgs.gov
    * http://www.shadedrelief.com/cleantopo2/index.html
    * https://www.usna.edu/Users/oceano/pguth/md_help/html/bathymetry.htm
    * https://www.opentopodata.org
    * https://en.wikipedia.org/wiki/Global_relief_model

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

    import os
    import numpy as np
    import scipy as sp
    from scipy import interpolate
    import matplotlib as mpl
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
    
    #x=type2numpy(x)

    #--------------------------------------------------------------------------
    #Check input values

    if xmin<-180 : xmin=-180
    if xmax>180 : xmax=180
    if ymin<-90 : ymin=-90
    if ymax>90 : ymax=90

    #--------------------------------------------------------------------------
    #Read data

    #Change folder to data folder
    Current_Folder=os.getcwd() #Current path
    os.chdir(filepath)

    #The xyz data
    #xyz = np.loadtxt('ETOPO1_XYZ_0_125_Degree.txt.gz')
    xyz = np.load('ETOPO1_XYZ_0_08_Deg_Python.npz')
    xyz = xyz['ETOPO1_XYZ_0_08_Deg_Python']
    xyz = xyz.astype(np.float64)
    x=xyz[:,0].copy() #x data
    y=xyz[:,1].copy() #y data
    z=xyz[:,2].copy() #z data

    #The z data
    #z = np.loadtxt('ETOPO1_Z_0_125_Deg_Python.txt.gz')
    #z = np.load('ETOPO1_Z_0_125_Deg_Python.npz')
    #z = z['ETOPO1_Z_0_125_Deg_Python']
    #z = z.astype(np.float64)

    #Return to current folder
    os.chdir(Current_Folder)

    #--------------------------------------------------------------------------
    #xyz data

    #The xy 2d and xy data
    dx=0.08
    dy=0.08
    xm=np.arange(-180,180+dx,dx)
    yn=np.arange(-90,90+dy,dy)
    Xnm_temp,Ynm_temp = np.meshgrid(xm,yn)
    #x=np.reshape(Xnm, -1)
    #y=np.reshape(Ynm, -1)

    #z 2d data
    M,N = Xnm_temp.shape
    Xnm = np.reshape(x, (M,N))
    Ynm = np.reshape(y, (M,N))
    Znm = np.reshape(z, (M,N))

    #The xyz data
    #xyz = np.column_stack((x,y,z))

    #--------------------------------------------------------------------------
    #Defining domain area and associated data for 1d data

    #Retain data within a defined domain

    #Finding a location data that are within [xmin:xmax] range
    xIndx=np.int_((np.nonzero((x>=xmin) & (x<=xmax)))[0])
    
    #Retaining data that are within [xmin:xmax] range
    x1=x[xIndx] #x data
    y1=y[xIndx] #y data
    z1=z[xIndx] #z data
    
    #Finding a location data that are within [ymin:ymax] range
    yIndx=np.int_((np.nonzero((y1>=ymin) & (y1<=ymax)))[0])
    
    #Retaining data that are within [ymin:ymax] range
    x2=x1[yIndx] #x data
    y2=y1[yIndx] #y data
    z2=z1[yIndx] #z data
    
    #--------------------------------------------------------------------------
    #Defining domain area and associated data for 2d data

    #Retain data within a 2d defined domain

    #Finding a location data that are within [xmin:xmax] range (bounding box)
    xIndx_2d=np.nonzero((Xnm>=xmin) & (Xnm<=xmax))
    i_start = xIndx_2d[0].min()
    j_start = xIndx_2d[1].min()
    i_end = xIndx_2d[0].max() + 1
    j_end = xIndx_2d[1].max() + 1
    
    #Retaining data that are within [xmin:xmax] range
    Xnm_1 = Xnm[i_start:i_end, j_start:j_end] #x grid data
    Ynm_1 = Ynm[i_start:i_end, j_start:j_end] #x grid data
    Znm_1 = Znm[i_start:i_end, j_start:j_end] #x grid data

    #Finding a location data that are within [ymin:ymax] range (bounding box)
    yIndx_2d=np.nonzero((Ynm_1>=ymin) & (Ynm_1<=ymax))
    i_start = yIndx_2d[0].min()
    j_start = yIndx_2d[1].min()
    i_end = yIndx_2d[0].max() + 1
    j_end = yIndx_2d[1].max() + 1
    
    #Retaining data that are within [ymin:ymax] range
    Xnm_2 = Xnm_1[i_start:i_end, j_start:j_end] #x grid data
    Ynm_2 = Ynm_1[i_start:i_end, j_start:j_end] #x grid data
    Znm_2 = Znm_1[i_start:i_end, j_start:j_end] #x grid data
    
    #--------------------------------------------------------------------------
    #Assigning new x, y and z

    x=x2.copy() #x data
    y=y2.copy() #y data
    z=z2.copy() #z data
    xyz=np.column_stack((x,y,z)) #xyz data

    xgrid=Xnm_2.copy() #x data
    ygrid=Ynm_2.copy() #y data
    zgrid=Znm_2.copy() #z data
    
    #--------------------------------------------------------------------------
    #Colormap, developed by Arash Karimpour
    #https://matplotlib.org/tutorials/colors/colormapnorms.html
    #https://matplotlib.org/tutorials/colors/colormap-manipulation.html

    cmap_water=mpl.cm.get_cmap('bone', 128)
    cmap_water=cmap_water(range(64,128,1))

    cmap_land_1=mpl.cm.get_cmap('summer', 10)
    cmap_land_1=cmap_land_1(range(0,10,1))

    cmap_land_2=mpl.cm.get_cmap('copper_r', 54)
    cmap_land_2=cmap_land_2(range(0,54,1))

    #cmap_z=np.vstack((cmap_water,cmap_land_1,cmap_land_2))

    #Assigning negative (water) and positive (land) colormap
    if (np.min(z)<0 and np.max(z)>0):
        cmap_z=np.vstack((cmap_water,cmap_land_1,cmap_land_2))
        z_lim=np.max([np.abs(np.max(z)),np.abs(np.min(z))])
        z_lim_min=-z_lim
        z_lim_max=z_lim
    elif (np.min(z)>=0 and np.max(z)>0):
        #cmap_z=np.vstack((cmap_land_1,cmap_land_2))
        cmap_z=cmap_land_2.copy()
        z_lim_min=np.min(z)
        z_lim_max=np.max(z)
    elif (np.min(z)<0 and np.max(z)<=0):
        cmap_z=cmap_water.copy()
        z_lim_min=np.min(z)
        z_lim_max=np.max(z)

    cmap_z_mpl=mpl.colors.ListedColormap(cmap_z,name='my_colormap',N=None)

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
        
        #plt.contourf(xgrid,ygrid,zgrid)
        
        plt.pcolormesh(xgrid,ygrid,zgrid,cmap=cmap_z_mpl,vmin=z_lim_min,vmax=z_lim_max,shading='nearest')

        #plt.imshow(zgrid,extent=[np.nanmin(xgrid),np.nanmax(xgrid),np.nanmin(ygrid),np.nanmax(ygrid)],interpolation=None,cmap=mpl.cm.get_cmap(),aspect='auto',origin='lower')

        plt.xlabel('Longitude (degree)')
        plt.ylabel('Latitude (degree)')
        plt.colorbar()
        
    #--------------------------------------------------------------------------
    #Outputs
    return x, y, z, xyz, xgrid, ygrid, zgrid

    #--------------------------------------------------------------------------
