def plot2dtimeseries(x, starttime='2000-10-20 00:00:00', endtime='2100-10-20 23:00:00', plottype='line', cmapcolors='blue', sizestyle='medium'):
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

    scientimate.plot2dtimeseries
    ============================

    .. code:: python

        scientimate.plot2dtimeseries(x, starttime='2000-10-20 00:00:00', endtime='2100-10-20 23:00:00', plottype='line', cmapcolors='blue', sizestyle='medium')

    Description
    -----------

    Plot x data in 2-d timeseries

    Inputs
    ------

    x
        x data as a 2-d array of (M,N)
    starttime='2000-10-20 00:00:00'
        | Start date and time for time series generation
        | Format: 'yyyy-mm-dd HH:MM:SS'
    endtime='2100-10-20 23:00:00'
        | End date and time for time series generation
        | Format: 'yyyy-mm-dd HH:MM:SS'
    plottype='bar'
        | Plot type
        | 'line': line plot
        | 'line_grid': line plot with grid lines
        | 'line_subplot': line plot with subplots
        | 'line_subplot_grid': line plot with subplots and grid lines
        | 'bar': bar plot
        | 'bar_grid': bar plot with grid lines
        | 'barh': horizontal bar plot
        | 'barh_grid': horizontal bar plot with grid lines
    cmapcolors='seq'
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

        x=np.random.rand(366)
        starttime='2020-01-01 00:00:00'
        endtime='2020-12-31 00:00:00' 
        sm.plot2dtimeseries(x,starttime,endtime,'line','purple','medium')

        x=np.random.rand(31,3)
        starttime='2020-01-01 00:00:00'
        endtime='2020-01-31 00:00:00'
        sm.plot2dtimeseries(x,starttime,endtime,'line_subplot','seq','medium')

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
    import datetime
    import matplotlib as mpl 
    import matplotlib.pyplot as plt 
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import colors 
    #import matplotlib.dates as mdates

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

    #--------------------------------------------------------------------------
    #Set plot drawing size

    if ((plottype=='line') or (plottype=='line_grid') \
        or (plottype=='line_subplot') or (plottype=='line_subplot_grid')):

        if sizestyle=='small':
            drawsize=1 #Size of plot
        elif sizestyle=='medium':
            drawsize=2 #Size of plot
        elif sizestyle=='large':
            drawsize=4 #Size of plot

    elif ((plottype=='bar') or (plottype=='bar_grid') \
        or (plottype=='barh') or (plottype=='barh_grid')):

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
    xaxislabel='time' #x axis label
    yaxislabel='y' #y axis label


    #--------------------------------------------------------------------------
    #Check input data

    y=np.array([])
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
    #Generate time vector

    timeformat='%Y-%m-%d %H:%M:%S' #Time format
    
    starttime=datetime.datetime.strptime(starttime,timeformat) #The datetime corresponding to date_string
    endtime=datetime.datetime.strptime(endtime,timeformat) #The datetime corresponding to date_string
    
    startdate=starttime.timestamp() #Convert start date and time to timestamp format
    enddate=endtime.timestamp() #Convert end date and time to timestamp format
    
    timedata=np.linspace(startdate,enddate,Mxdata) #Time series in timestamp format
    timedata=[datetime.datetime.fromtimestamp(elm) for elm in timedata] #Time series in datetime format
    dt=timedata[1]-timedata[0] #Delta t
    
    datestring=[elm.strftime(timeformat) for elm in timedata] #Time series in string format


    #--------------------------------------------------------------------------
    #The line plot

    if ((plottype=='line') or (plottype=='line_grid')):

        for i in range(0,Nxdata,1):
            plt.plot(timedata,ydata[:,i],color=cmap_ncolor[i,:],linewidth=drawsize)


    #--------------------------------------------------------------------------
    #The line plot with subplots

    if ((plottype=='line_subplot') or (plottype=='line_subplot_grid')):

        plotfontsize=18 #Size of plot fonts
        axisfontsize=18 #Size of axis fonts
        for i in range(0,Nxdata,1):
            plt.subplot(Nxdata,1,i+1)
            plt.plot(timedata,ydata[:,i],color=cmap_ncolor[i,:],linewidth=drawsize)
            plt.title(str(i+1))
            plt.xlim((timedata[0],timedata[-1]))
            plt.ylim((np.nanmin(y),np.nanmax(y)))
            plt.rcParams.update({'font.size':plotfontsize})
            plt.xlabel(xaxislabel,fontsize=axisfontsize)
            plt.ylabel(yaxislabel,fontsize=axisfontsize)

            ax=plt.gca()
            tick_locator=mpl.dates.AutoDateLocator()
            tick_format=mpl.dates.AutoDateFormatter(tick_locator)
            ax.xaxis.set_major_locator(tick_locator)
            ax.xaxis.set_major_formatter(tick_format)

            #years=mpl.dates.YearLocator() #every year
            #months=mpl.dates.MonthLocator() #every month
            #tick_format=mpl.dates.DateFormatter('%Y-%m-%d')
            #ax.xaxis.set_major_locator(years)
            #ax.xaxis.set_major_formatter(tick_format)
            #ax.xaxis.set_minor_locator(months)


    #--------------------------------------------------------------------------
    #The line plot with subplots

    if plottype=='line_subplot_grid':

        plotfontsize=18 #Size of plot fonts
        axisfontsize=18 #Size of axis fonts
        for i in range(0,Nxdata,1):
            plt.subplot(Nxdata,1,i+1)
            plt.plot(timedata,ydata[:,i],color=cmap_ncolor[i,:],linewidth=drawsize)
            plt.grid()
            plt.title(str(i+1))
            plt.xlim((timedata[0],timedata[-1]))
            plt.ylim((np.nanmin(y),np.nanmax(y)))
            plt.rcParams.update({'font.size':plotfontsize})
            plt.xlabel(xaxislabel,fontsize=axisfontsize)
            plt.ylabel(yaxislabel,fontsize=axisfontsize)

            ax=plt.gca()
            fig=plt.gcf()
            tick_locator=mpl.dates.AutoDateLocator()
            tick_format=mpl.dates.AutoDateFormatter(tick_locator)
            ax.xaxis.set_major_locator(tick_locator)
            ax.xaxis.set_major_formatter(tick_format)
            fig.autofmt_xdate()

            #years=mpl.dates.YearLocator() #every year
            #months=mpl.dates.MonthLocator() #every month
            #tick_format=mpl.dates.DateFormatter('%Y-%m-%d')
            #ax.xaxis.set_major_locator(years)
            #ax.xaxis.set_major_formatter(tick_format)
            #ax.xaxis.set_minor_locator(months)


    #--------------------------------------------------------------------------
    #Plot bar plot

    if ((plottype=='bar') or (plottype=='bar_grid')):

        if Nxdata==1:
            plt.bar(timedata,ydata[:,0],width=drawsize,color=cmap_ncolor[0,:],edgecolor=None)
        else:
            plt.bar(timedata,ydata[:,0],width=drawsize,color=cmap_ncolor[0,:],edgecolor=None)


    #--------------------------------------------------------------------------
    #Plot bar plot

    #The bar plot
    if ((plottype=='barh') or (plottype=='barh_grid')):

        if Nxdata==1:
            plt.barh(timedata,ydata[:,0],height=drawsize,color=cmap_ncolor[0,:],edgecolor=None)
        else:
            plt.barh(timedata,ydata[:,0],height=drawsize,color=cmap_ncolor[0,:],edgecolor=None)


    #--------------------------------------------------------------------------
    #Set plot properties

    if ((plottype!='line_subplot') and (plottype!='line_subplot_grid')):

        #Set time axis ticks
        #set(gca,'XTick',ceil(startdate):dt:floor(enddate))
        if ((plottype=='barh') or (plottype=='barh_grid')):
            ax=plt.gca()
            tick_locator=mpl.dates.AutoDateLocator()
            tick_format=mpl.dates.AutoDateFormatter(tick_locator)
            ax.yaxis.set_major_locator(tick_locator)
            ax.yaxis.set_major_formatter(tick_format)

            #years=mpl.dates.YearLocator() #every year
            #months=mpl.dates.MonthLocator() #every month
            #tick_format=mdates.DateFormatter('%Y-%m-%d')
            #ax.yaxis.set_major_locator(years)
            #ax.yaxis.set_major_formatter(tick_format)
            #ax.yaxis.set_minor_locator(months)
        else:
            ax=plt.gca()
            tick_locator=mpl.dates.AutoDateLocator()
            tick_format=mpl.dates.AutoDateFormatter(tick_locator)
            ax.xaxis.set_major_locator(tick_locator)
            ax.xaxis.set_major_formatter(tick_format)

            #years=mpl.dates.YearLocator() #every year
            #months=mpl.dates.MonthLocator() #every month
            #tick_format=mpl.dates.DateFormatter('%Y-%m-%d')
            #ax.xaxis.set_major_locator(years)
            #ax.xaxis.set_major_formatter(tick_format)
            #ax.xaxis.set_minor_locator(months)
            

        #Set axis limits
        if ((plottype=='barh') or (plottype=='barh_grid')):
            plt.xlim((np.nanmin(y),np.nanmax(y)))
            plt.ylim((timedata[0],timedata[-1]))
        else:
            plt.xlim((timedata[0],timedata[-1]))
            plt.ylim((np.nanmin(y),np.nanmax(y)))

        #Plot grid lines
        if ((plottype=='line_grid') or (plottype=='line_subplot_grid') \
            or (plottype=='bar_grid') or (plottype=='barh_grid')):
            plt.grid()

        #Set figure font size
        plt.rcParams.update({'font.size':plotfontsize})

        #Plot axis label
        if ((plottype=='barh') or (plottype=='barh_grid')):
            plt.xlabel(yaxislabel,fontsize=axisfontsize)
            plt.ylabel(xaxislabel,fontsize=axisfontsize)
        else:
            plt.xlabel(xaxislabel,fontsize=axisfontsize)
            plt.ylabel(yaxislabel,fontsize=axisfontsize)

        #Plot legends
    legend_all_names=['$y_1$','$y_2$','$y_3$','$y_4$','$y_5$','$y_6$','$y_7$','$y_8$','$y_9$']
    if ((Nxdata>=2) and (Nxdata<=3)):
        legend_names=legend_all_names[0:Nxdata]
        plt.legend(legend_names)


    #--------------------------------------------------------------------------
