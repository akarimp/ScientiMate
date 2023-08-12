def probability2d(x, y, binedgex=None, binedgey=None, dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2017-06-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.probability2d
    =========================

    .. code:: python

        fdensityxy, fdensityx, fdensityy, fdensitycumulativex, fdensitycumulativey, bincenterx, bincentery, xmean, ymean, xstd, ystd = scientimate.probability2d(x, y, binedgex=None, binedgey=None, dispout='no')

    Description
    -----------

    Calculate 2D (joint) probability density distribution for two given datasets

    Inputs
    ------

    x
        First input dataset 
    y
        Second input dataset 
    binedgex
        | Bin edges for x data  
        | length(binedgex)=number of bin +1   
        | If there are N bins in histogram/distribution, then values in binedgex are as:   
        | 1st value: left edge of 1st bin, 2nd value: left edge of 2nd bin, ...   
        | (N)th value: left edge of last bin, (N+1)th value: right edge of last bin   
    binedgey
        | Bin edge for y data  
        | length(binedgey)=number of bin +1   
        | If there are N bins in histogram/distribution, then values in binedgey are as:   
        | 1st value: left edge of 1st bin, 2nd value: left edge of 2nd bin, ...   
        | (N)th value: left edge of last bin, (N+1)th value: right edge of last bin   
    dispout='no'
        | Define to display outputs or not
        | 'no': not display 
        | 'bar': bar plot
        | 'imagesc': 2 dimensional plot using imagesc or imshow
        | 'pcolor': 2 dimensional plot using pcolor
        | 'contour': 2 dimensional contour plot, number of contour=32
        | 'contourf': 2 dimensional filled contour plot, number of contour=32
        | 'surface': 3 dimensional surface plot 
        | 'binscatter': 2 dimensional histogram plot using imagesc or imshow

    Outputs
    -------

    fdensityxy
        2D (joint) probability density distribution for x-y data
    fdensityx
        Probability density distribution for x data
    fdensityy
        Probability density distribution for y data
    fdensitycumulativex
        Cumulative probability density distribution for x data
    fdensitycumulativey
        Cumulative probability density distribution for y data
    bincenterx
        Bin center for x data
    bincentery
        Bin center for y data
    xmean
        Mean value of x data
    ymean
        Mean value of y data
    xstd
        Standard deviation of x data
    ystd
        Standard deviation of y data

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        x=(-0.1+(0.1-(-0.1)))*np.random.randn(1024*2)
        y=(-0.2+(0.2-(-0.2)))*np.random.randn(1024*2)
        binedgex=np.linspace(min(x),max(x),11)
        binedgey=np.linspace(min(y),max(y),11)
        fdensityxy,fdensityx,fdensityy,fdensitycumulativex,fdensitycumulativey,bincenterx,bincentery,xmean,ymean,xstd,ystd=sm.probability2d(x,y,binedgex,binedgey,'surface')

        x=np.random.randn(100000)
        y=1.5*x+np.random.randn(100000)
        binedgex=np.linspace(min(x),max(x),101)
        binedgey=np.linspace(min(y),max(y),101)
        sm.probability2d(x,y,binedgex,binedgey,'binscatter')

    References
    ----------


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
    if dispout!='no':
        import matplotlib.pyplot as plt 
        from mpl_toolkits.mplot3d import Axes3D

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

    #--------------------------------------------------------------------------
    #Assign default values

    if binedgex is not None: binedgex=type2numpy(binedgex)
    if binedgex is None: binedgex=np.linspace(np.nanmin(x),np.nanmax(x),11)

    if binedgey is not None: binedgey=type2numpy(binedgey)
    if binedgey is None: binedgey=np.linspace(np.nanmin(x),np.nanmax(x),11)

    #--------------------------------------------------------------------------
    #Removing NaN and Inf values

    #Total number of elements in dataset
    NTotal=len(x)

    #Finding NaN and Inf in x dataset
    Indx1=np.int_((np.nonzero((np.isnan(x)!=1) & (np.isinf(x)!=1)))[0])
    x1=x[Indx1]
    y1=y[Indx1]
    
    #Finding NaN and Inf in y dataset
    Indx2=np.int_((np.nonzero((np.isnan(y1)!=1) & (np.isinf(y1)!=1)))[0])
    x2=x1[Indx2]
    y2=y1[Indx2]
    
    #Datasets without NaN and Inf
    x=x2.copy()
    y=y2.copy()

    #Total number of elements in dataset after removal of NaN and Inf
    N=len(x)
    
    #Total number of NaN and Inf elements in dataset
    NNaNInf=NTotal-N

    #--------------------------------------------------------------------------
    #Calculating probability density distribution for x data

    binwidthx=np.diff(binedgex) #Bin width for x data, length(binwidth) is one element less than length(binedge) 
    bincenterx=binedgex[0:-1]+binwidthx/2 #Bin center for x data
    
    NPointsx,binedgeoutx=np.histogram(x,binedgex) #Calculating number of data points in each bin for x data
    histareax=NPointsx*binwidthx #Area under each element of histogram for x data
    fdensityx=NPointsx/(np.sum(histareax)) #Calculating probability density distribution for x data
    #Note: np.sum(fdensityx*binwidthx)=1

    #--------------------------------------------------------------------------
    #Calculating cumulative probability density distribution for x data

    fdensitycumulativex=np.zeros(len(bincenterx)) #Pre-assigning vector
    for i in range(0,len(bincenterx),1):
        fdensitycumulativex[i]=np.sum(fdensityx[0:i+1]*binwidthx[0:i+1])

    #--------------------------------------------------------------------------
    #Calculating statistical properties for x data

    xmean=np.mean(x) #Mean value of x data
    xstd=np.std(x,ddof=1) #Standard deviation of x data

    #--------------------------------------------------------------------------
    #Calculating probability density distribution for y data

    binwidthy=np.diff(binedgey) #Bin width for y data, length(binwidth) is one element less than length(binedge) 
    bincentery=binedgey[0:-1]+binwidthy/2 #Bin center for y data
    
    NPointsy,binedgeouty=np.histogram(y,binedgey) #Calculating number of data points in each bin for y data
    histareay=NPointsy*binwidthy #Area under each element of histogram for y data
    fdensityy=NPointsy/(np.sum(histareay)) #Calculating probability density distribution for y data
    #Note: np.sum(fdensityy*binwidthy)=1

    #--------------------------------------------------------------------------
    #Calculating cumulative probability density distribution for y data

    fdensitycumulativey=np.zeros(len(bincentery)) #Pre-assigning vector
    for i in range(0,len(bincentery),1):
        fdensitycumulativey[i]=np.sum(fdensityy[0:i+1]*binwidthy[0:i+1])

    #--------------------------------------------------------------------------
    #Calculating statistical properties for y data

    ymean=np.mean(y) #Mean value of y data
    ystd=np.std(y,ddof=1) #Standard deviation of y data

    #--------------------------------------------------------------------------
    #Calculating 2D (joint) probability density distribution for x and y

    NPointsxy,binedgeoutx,binedgeouty=np.histogram2d(x,y,[binedgex,binedgey])
    
    binwidthxGrid,binwidthyGrid=np.meshgrid(binwidthx,binwidthy)
    histxyvolume=NPointsxy*binwidthxGrid.T*binwidthyGrid.T #Volume under each element of histogram
    fdensityxy=NPointsxy/(np.sum(histxyvolume)) #Calculating 2D (joint) probability density distribution for x-y data
    #Note: np.sum(fdensityxy*binwidthxGrid.T*binwidthyGrid.T)=1
    
    #fdensityxy=np.zeros((len(bincenterx),len(bincentery))) #Pre-assigning vector
    #for i in range(0,len(bincenterx),1):
    #    for j in range(0,len(bincentery),1):
    #        fdensityxy[i,j]=NPointsxy[i,j]/(np.sum(histxyvolume)) #Calculating 2D (joint) probability density distribution for x-y data
    #        #Note: np.sum(fdensityxy*binwidthxGrid.T*binwidthyGrid.T)=1

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='bar': #Bar plot
    
        val=[xmean,xstd,ymean,ystd,NTotal,N,NNaNInf]
        name=['x mean','x st. deviation','y mean','y st. dev.','N Input','N Used','N NaN Inf']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            #print('{0:10}= {1:0.10f}'.format(name[i],val[i])) 
            print('     {0:15}     {1:+0.10g}'.format(name[i],val[i])) 

        #Plotting x Data
        plt.subplot(2,2,1)
        plt.bar(bincenterx,fdensityx,width=0.8*binwidthx)
        plt.title('Probability Density Distribution')
        plt.ylabel('Probability Density, fx')
        plt.xlabel('x Data')
        
        plt.subplot(2,2,3)
        plt.plot(bincenterx,fdensitycumulativex)
        plt.title('Cumulative Probability Density Distribution')
        plt.ylabel('Cumulative Probability Density, fx')
        plt.xlabel('x Data')
        
        #Plotting y Data
        plt.subplot(2,2,2)
        plt.bar(bincentery,fdensityy,width=0.8*binwidthy)
        plt.title('Probability Density Distribution') #Title of plot
        plt.ylabel('Probability Density, fy')
        plt.xlabel('y Data')
    
        plt.subplot(2,2,4)
        plt.plot(bincentery,fdensitycumulativey)
        plt.title('Cumulative Probability Density Distribution')
        plt.ylabel('Cumulative Probability Density, fy')
        plt.xlabel('y Data')
    
    elif dispout=='imagesc': #Imagesc plot
        
        val=[xmean,xstd,ymean,ystd,NTotal,N,NNaNInf]
        name=['x mean','x st. deviation','y mean','y st. dev.','N Input','N Used','N NaN Inf']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            print('     {0:15}     {1:+0.10g}'.format(name[i],val[i])) 
    
        bincenterxGrid,bincenteryGrid=np.meshgrid(bincenterx,bincentery)
        plt.imshow(fdensityxy.T,extent=[np.nanmin(bincenterxGrid),np.nanmax(bincenterxGrid),np.nanmin(bincenteryGrid),np.nanmax(bincenteryGrid)],interpolation=None,cmap=plt.get_cmap(),aspect='auto',origin='lower')
        plt.colorbar()
        plt.title('2D Probability Density Distribution')
        plt.xlabel('x Data')
        plt.ylabel('y Data')

    elif dispout=='pcolor': #Pcolor plot

        val=[xmean,xstd,ymean,ystd,NTotal,N,NNaNInf]
        name=['x mean','x st. deviation','y mean','y st. dev.','N Input','N Used','N NaN Inf']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            print('     {0:15}     {1:+0.10g}'.format(name[i],val[i])) 
    
        bincenterxGrid,bincenteryGrid=np.meshgrid(bincenterx,bincentery)
        plt.pcolormesh(bincenterxGrid,bincenteryGrid,fdensityxy.T,cmap=plt.get_cmap())    
        plt.colorbar()
        plt.title('2D Probability Density Distribution')
        plt.xlabel('x Data')
        plt.ylabel('y Data')

    elif dispout=='contour': #Contour plot
    
        val=[xmean,xstd,ymean,ystd,NTotal,N,NNaNInf]
        name=['x mean','x st. deviation','y mean','y st. dev.','N Input','N Used','N NaN Inf']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            print('     {0:15}     {1:+0.10g}'.format(name[i],val[i])) 
    
        bincenterxGrid,bincenteryGrid=np.meshgrid(bincenterx,bincentery)
        count1=plt.contour(bincenterxGrid,bincenteryGrid,fdensityxy.T,cmap=plt.get_cmap())
        plt.colorbar(count1)
        plt.title('2D Probability Density Distribution')
        plt.xlabel('x Data')
        plt.ylabel('y Data')

    elif dispout=='contourf': #Contourf plot
    
        val=[xmean,xstd,ymean,ystd,NTotal,N,NNaNInf]
        name=['x mean','x st. deviation','y mean','y st. dev.','N Input','N Used','N NaN Inf']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            print('     {0:15}     {1:+0.10g}'.format(name[i],val[i])) 
    
        bincenterxGrid,bincenteryGrid=np.meshgrid(bincenterx,bincentery)
        count1=plt.contourf(bincenterxGrid,bincenteryGrid,fdensityxy.T,cmap=plt.get_cmap())
        plt.colorbar(count1)
        plt.title('2D Probability Density Distribution')
        plt.xlabel('x Data')
        plt.ylabel('y Data')

    elif dispout=='surface': #Surface plot
    
        val=[xmean,xstd,ymean,ystd,NTotal,N,NNaNInf]
        name=['x mean','x st. deviation','y mean','y st. dev.','N Input','N Used','N NaN Inf']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            print('     {0:15}     {1:+0.10g}'.format(name[i],val[i])) 
    
        bincenterxGrid,bincenteryGrid=np.meshgrid(bincenterx,bincentery)

        fig=plt.figure()
        fig3d=fig.gca(projection='3d')
        #fig3d.view_init(elev=50,azim=-20)
        surf1=fig3d.plot_surface(bincenterxGrid,bincenteryGrid,fdensityxy.T,cmap=plt.get_cmap(),edgecolor='none')
        fig.colorbar(surf1)
        fig3d.set_title('2D Probability Density Distribution')
        fig3d.set_xlabel('x Data')
        fig3d.set_ylabel('y Data')
        fig3d.set_zlabel('2D Probability Density, f_{x,y}')

    elif dispout=='binscatter': #Binscatter plot
        
        val=[xmean,xstd,ymean,ystd,NTotal,N,NNaNInf]
        name=['x mean','x st. deviation','y mean','y st. dev.','N Input','N Used','N NaN Inf']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            print('     {0:15}     {1:+0.10g}'.format(name[i],val[i])) 
    
        plt.imshow(NPointsxy.T,extent=[np.nanmin(binedgeoutx),np.nanmax(binedgeoutx),np.nanmin(binedgeouty),np.nanmax(binedgeouty)],interpolation=None,cmap=plt.get_cmap(),aspect='auto',origin='lower')
        plt.colorbar()
        plt.title('2D Histogram')
        plt.xlabel('x Data')
        plt.ylabel('y Data')

    #--------------------------------------------------------------------------
    #Outputs
    return fdensityxy, fdensityx, fdensityy, fdensitycumulativex, fdensitycumulativey, bincenterx, bincentery, xmean, ymean, xstd, ystd

    #--------------------------------------------------------------------------
