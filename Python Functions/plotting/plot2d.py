def plot2d(x, y=None, plottype='line', cmapcolors='blue', sizestyle='medium'):
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

    scientimate.plot2d
    ==================

    .. code:: python

        scientimate.plot2d(x, y=None, plottype='line', cmapcolors='blue', sizestyle='medium')

    Description
    -----------

    Plot x and y data in 2-d plot

    Inputs
    ------

    x
        x data as a 2-d array of (M,N)
    y=[]
        y data as a 2-d array of (M,N)
    plottype='line'
        | Plot type
        | 'line': line plot
        | 'line_grid': line plot with grid lines
        | 'line_ascend': line plot with ascending line width
        | 'line_ascend_grid': line plot with ascending line width and grid lines
        | 'line_confid': line plot with 95% confidence intervals band (approximated)
        | 'line_confid_grid': line plot with 95% confidence intervals band (approximated) and grid lines
        | 'scatter': scatter plot
        | 'scatter_grid': scatter plot with grid lines
        | 'scatter_ascend': scatter plot with ascending point size
        | 'scatter_ascend_grid': scatter plot with grid lines and ascending point size
        | 'bar': bar plot
        | 'bar_grid': bar plot with grid lines
        | 'bar_stacked': stacked bar plot
        | 'barh': horizontal bar plot
        | 'barh_grid': horizontal bar plot with grid lines
        | 'barh_stacked': horizontal stacked bar plot
        | 'histogram': histogram plot
        | 'histogram_grid': histogram plot with grid lines
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
    sizestyle='medium'
        | Plot drawing size style
        | 'small': small plot size
        | 'medium': medium plot size
        | 'large': large plot size

    Outputs
    -------


    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        x=[0,1]
        y=np.zeros((2,50))
        y[0,:]=np.linspace(2,51,50)
        y[1,:]=np.linspace(2,51,50)
        sm.plot2d(x,y,'line','blue_red','large')

        x=np.linspace(1,10,10)
        y=np.zeros((10,2))
        y[:,0]=1+np.random.rand(10)
        y[:,1]=2+np.random.rand(10)
        sm.plot2d(x,y,'line_confid','blue_red','large')

        x=np.linspace(0,2*np.pi,1000)
        x=np.tile(x[:,np.newaxis],(1,10))
        s=np.arange(1,11)
        s=np.tile(s[np.newaxis,:],(1000,1))
        y=s+np.sin(1.0*np.pi*x)
        sm.plot2d(x,y,'line','cool','large')

        x=np.random.rand(100,3)
        y=np.random.rand(100,3)
        y=np.zeros((100,3))
        y[:,0]=1+2.0*x[:,0]+np.random.rand(100)
        y[:,1]=3+2.0*x[:,1]+np.random.rand(100)
        y[:,2]=5+2.0*x[:,2]+np.random.rand(100)
        sm.plot2d(x,y,'scatter','seq','large')

        x=np.random.rand(100)
        y=np.random.rand(100)
        sm.plot2d(x,y,'scatter_ascend','purple','large')

        x=[[1,1,1],[2,2,2],[3,3,3],[4,4,4]]
        y=[[2,3,8],[2,5,6],[5,7,9],[1,2,3]]
        sm.plot2d(x,y,'bar','purple','medium')

        x=[1,3,5,7,9,11,13,15]
        y=[2,3,9,8,2,5,6,9]
        sm.plot2d(x,y,'bar','purple','medium')

        x=np.random.randn(1000)
        sm.plot2d(x,[],'histogram','purple','medium')

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
    if y is not None: y=type2numpy(y) #Check if y is empty

    #--------------------------------------------------------------------------
    #Assigning default values

    if y is None: y=np.array([]) #Check if y is empty

    #--------------------------------------------------------------------------
    #Set plot drawing size

    if ((plottype=='line') or (plottype=='line_grid') \
        or (plottype=='line_ascend') or (plottype=='line_ascend_grid') \
        or (plottype=='line_confid') or (plottype=='line_confid_grid')):

        if sizestyle=='small':
            drawsize=1 #Size of plot
        elif sizestyle=='medium':
            drawsize=2 #Size of plot
        elif sizestyle=='large':
            drawsize=4 #Size of plot


    elif ((plottype=='scatter') or (plottype=='scatter_grid') \
        or (plottype=='scatter_ascend') or (plottype=='scatter_ascend_grid')):

        if sizestyle=='small':
            drawsize=18 #Size of plot
        elif sizestyle=='medium':
            drawsize=36 #Size of plot
        elif sizestyle=='large':
            drawsize=100 #Size of plot


    elif ((plottype=='bar') or (plottype=='bar_grid') or (plottype=='bar_stacked') \
        or (plottype=='barh') or (plottype=='barh_grid') or (plottype=='barh_stacked') \
        or (plottype=='histogram') or (plottype=='histogram_grid')):

        if sizestyle=='small':
            drawsize=0.5 #Size of plot
        elif sizestyle=='medium':
            drawsize=0.8 #Size of plot
        elif sizestyle=='large':
            drawsize=1 #Size of plot


    plotfontsize=18 #Size of plot fonts
    axisfontsize=18 #Size of axis fonts


    #--------------------------------------------------------------------------
    #Input variable names

    #At the moment there in no way to find an exaxt name of input arguments that passed to function
    #Python insert module provide functions to chack some input properties
    xaxislabel='x' #x axis label
    yaxislabel='y' #y axis label
    #cbarlabel='z' #Colorbar label


    #--------------------------------------------------------------------------
    #Check input data

    #Create x and y if y=[]
    if np.size(y)==0: #for Python list: if not y==True
        Mx=np.shape(x)[0]
        y=x.copy()
        x=np.arange(1,Mx+1,1)


    #Change inputs to 2d array
    if np.ndim(x)==1:
        x=x[:,np.newaxis]

    if np.ndim(y)==1:
        y=y[:,np.newaxis]


    #Input data shape
    Mx,Nx=np.shape(x)
    My,Ny=np.shape(y)


    #--------------------------------------------------------------------------
    #Reshape input data to (M,N)

    #x and y have the same shape of (M,N)
    if ((Mx==My) and (Nx==Ny)):
        xdata=x.copy()
        ydata=y.copy()

    #x has a shape of (M,1) and y has a shape of (M,N)
    elif ((Mx==My) and (Nx==1)):
        xdata=np.tile(x,(1,Ny))
        ydata=y.copy()

    #x has a shape of (M,1) and y has a shape of (N,M)
    elif ((Mx==Ny) and (Nx==1)):
        xdata=np.tile(x,(1,My))
        ydata=np.transpose(y)

    #x has a shape of (1,N) and y has a shape of (M,N)
    elif ((Mx==1) and (Nx==Ny)):
        xdata=np.transpose(x)
        xdata=np.tile(xdata,(1,My))
        ydata=np.transpose(y)

    #x has a shape of (1,N) and y has a shape of (N,M)
    elif ((Mx==1) and (Nx==My)):
        xdata=np.transpose(x)
        xdata=np.tile(xdata,(1,Ny))
        ydata=y.copy()


    #Shape of reshaped data
    Mxdata,Nxdata=np.shape(xdata)
    Mydata,Nydata=np.shape(ydata)


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

    if Nxdata==1:
        ncolor=3 #Number of zclmap to be used in colormap
    else:
        ncolor=Nxdata #Number of zclmap to be used in colormap


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


    if Nxdata==1:
        cmap_ncolor=cmap_ncolor[1:,:] #Remove first (light) color


    #--------------------------------------------------------------------------
    #The line plot

    if ((plottype=='line') or (plottype=='line_grid')):

        for i in range(0,Nxdata,1):
            plt.plot(xdata[:,i],ydata[:,i],color=cmap_ncolor[i,:],linewidth=drawsize)


    #--------------------------------------------------------------------------
    #The line plot with ascending line width

    if ((plottype=='line_ascend') or (plottype=='line_ascend_grid')):

        drawsize=np.linspace(1,4,Nxdata) #Default 'LineWidth' is 0.5
        for i in range(0,Nxdata,1):
            plt.plot(xdata[:,i],ydata[:,i],color=cmap_ncolor[i,:],linewidth=drawsize[i])


    #--------------------------------------------------------------------------
    #The line plot with confidence interval

    if ((plottype=='line_confid') or (plottype=='line_confid_grid')):

        #Plot lines
        for i in range(0,Nxdata,1):
            plt.plot(xdata[:,i],ydata[:,i],color=cmap_ncolor[i,:],linewidth=drawsize)

        #Plot confidence intervals
        for i in range(0,Nxdata,1):
            sem=np.std(ydata[:,i])/np.sqrt(len(ydata[:,i])) #Standard error
            ci_up=ydata[:,i]+1.96*sem #Confidence interval upper level
            ci_low=ydata[:,i]-1.96*sem #Confidence interval lower level
            x_fill=np.concatenate((xdata[:,i],np.flipud(xdata[:,i])))
            y_fill=np.concatenate((ci_up,np.flipud(ci_low)))
            plt.fill_between(x_fill,y_fill,facecolor=cmap_ncolor[i,:],edgecolor=None,alpha=0.4)


    #--------------------------------------------------------------------------
    #The scatter plot

    if ((plottype=='scatter') or (plottype=='scatter_grid')):

        for i in range(0,Nxdata,1):
            plt.scatter(xdata[:,i],ydata[:,i],s=drawsize,color=cmap_ncolor[i,:])


    #--------------------------------------------------------------------------
    #The scatter plot with ascending point size

    if ((plottype=='scatter_ascend') or (plottype=='scatter_ascend_grid')):

        if Nxdata==1:
            drawsize=np.linspace(36,100,Mxdata) #Default size is 36
            plt.scatter(xdata[:,0],ydata[:,0],s=drawsize,color=cmap_ncolor[0,:])
        
        else:
            drawsize=np.linspace(36,100,Nxdata) #Default size is 36
            for i in range(0,Nxdata,1):
                plt.scatter(xdata[:,i],ydata[:,i],s=drawsize,color=cmap_ncolor[i,:])


    #--------------------------------------------------------------------------
    #The bar plot

    if ((plottype=='bar') or (plottype=='bar_grid')):

        if Nxdata==1:
            plt.bar(xdata[:,0],ydata[:,0],width=drawsize,color=cmap_ncolor[0,:],edgecolor=None)
        else:
            for i in range(0,Nxdata,1):
                xdata_shift=-drawsize/2+i/float(Nxdata)*drawsize
                plt.bar(xdata[:,i]+xdata_shift,ydata[:,i],width=drawsize/float(Nxdata),color=cmap_ncolor[i,:],edgecolor=None)


    #--------------------------------------------------------------------------
    #The stacked bar plot

    if plottype=='bar_stacked':

        plt.bar(xdata[:,0],ydata[:,0],width=drawsize,color=cmap_ncolor[0,:],edgecolor=None)
        for i in range(1,Nxdata,1):
            plt.bar(xdata[:,i],ydata[:,i],width=drawsize,bottom=ydata[:,i-1],color=cmap_ncolor[i,:],edgecolor=None)


    #--------------------------------------------------------------------------
    #The horizontal bar plot

    if ((plottype=='barh') or (plottype=='barh_grid')):

        if Nxdata==1:
            plt.barh(xdata[:,0],ydata[:,0],height=drawsize,color=cmap_ncolor[0,:],edgecolor=None)
        else:
            for i in range(0,Nxdata,1):
                xdata_shift=-drawsize/2+i/float(Nxdata)*drawsize
                plt.barh(xdata[:,i]+xdata_shift,ydata[:,i],height=drawsize/float(Nxdata),color=cmap_ncolor[i,:],edgecolor=None)


    #--------------------------------------------------------------------------
    #The horizontal stacked bar plot

    if plottype=='barh_stacked':
    
        plt.barh(xdata[:,0],ydata[:,0],height=drawsize,color=cmap_ncolor[0,:],edgecolor=None)
        for i in range(1,Nxdata,1):
            plt.barh(xdata[:,i],ydata[:,i],height=drawsize,left=ydata[:,i-1],color=cmap_ncolor[i,:],edgecolor=None)


    #--------------------------------------------------------------------------
    #The histogram plot

    if ((plottype=='histogram') or (plottype=='histogram_grid')):
    
        plt.hist(ydata,rwidth=drawsize,color=cmap_ncolor[0,:],edgecolor=None)


    #--------------------------------------------------------------------------
    #Set plot properties

    #Set axis limits
    if ((plottype=='bar') or (plottype=='bar_grid') or (plottype=='bar_stacked') \
        or (plottype=='barh') or (plottype=='barh_grid') or (plottype=='barh_stacked') \
        or (plottype=='histogram') or (plottype=='histogram_grid')):
        dummy=1 #Do nothing
    else:
        plt.xlim((np.nanmin(x),np.nanmax(x)))
        plt.ylim((np.nanmin(y),np.nanmax(y)))

    #Plot grid lines
    if ((plottype=='line_grid') or (plottype=='line_ascend_grid') or (plottype=='line_confid_grid') \
        or (plottype=='scatter_grid') or (plottype=='scatter_ascend_grid') \
        or (plottype=='bar_grid') or (plottype=='barh_grid') or (plottype=='histogram_grid')):
        plt.grid()

    #Set figure font size
    plt.rcParams.update({'font.size':plotfontsize})

    #Set axis label
    if plottype=='histogram':
        plt.xlabel(xaxislabel,fontsize=axisfontsize)
        plt.ylabel('Counts',fontsize=axisfontsize)
    else:
        plt.xlabel(xaxislabel,fontsize=axisfontsize)
        plt.ylabel(yaxislabel,fontsize=axisfontsize)
    
    #Plot legends
    legend_all_names=['$y_1$','$y_2$','$y_3$','$y_4$','$y_5$','$y_6$','$y_7$','$y_8$','$y_9$']
    if ((Nxdata>=2) and (Nxdata<=3)):
        legend_names=legend_all_names[0:Nxdata]
        plt.legend(legend_names)


    #--------------------------------------------------------------------------