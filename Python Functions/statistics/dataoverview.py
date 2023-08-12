def dataoverview(x, y=[], z=[]):
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

    scientimate.dataoverview
    ========================

    .. code:: python

        scientimate.dataoverview(x, y=[], z=[])

    Description
    -----------

    Display an overview of the input data

    Inputs
    ------

    x
        x data
    y=[]
        y data (Optional)
    z=[]
        z data (Optional)

    Outputs
    -------


    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        x=(-0.1+(0.1-(-0.1)))*np.random.randn(1024*2)
        sm.dataoverview(x)

        x=(-0.1+(0.1-(-0.1)))*np.random.randn(1024*2)
        y=x+(-0.01+(0.01-(-0.01)))*np.random.randn(1024*2)
        sm.dataoverview(x,y)

        x=(-0.1+(0.1-(-0.1)))*np.random.randn(1024*2)
        y=(-0.2+(0.2-(-0.2)))*np.random.randn(1024*2)
        z=x**2+y**2+0.1+(-0.1+(0.1-(-0.1)))*np.random.rand(1024*2)
        sm.dataoverview(x,y,z)

    References
    ----------

    * https://docs.python.org/3/library/string.html

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
    z=type2numpy(z)

    #--------------------------------------------------------------------------
    #Removing NaN and Inf values

    #Total number of elements in dataset
    NTotal=len(x)

    if ((np.size(x)>0) and (np.size(y)==0) and (np.size(z)==0)): #Only x is available
    
        #Finding NaN and Inf in x dataset
        Indx1=np.int_((np.nonzero((np.isnan(x)!=1) & (np.isinf(x)!=1)))[0])
        x1=x[Indx1]
    
        #Datasets without NaN and Inf
        x=x1.copy()
    
    elif ((np.size(x)>0) and (np.size(y)>0) and (np.size(z)==0)): #x and y are available
    
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
    
    elif ((np.size(x)>0) and (np.size(y)>0) and (np.size(z)>0)): #x, y and z are available
    
        #Finding NaN and Inf in x dataset
        Indx1=np.int_((np.nonzero((np.isnan(x)!=1) & (np.isinf(x)!=1)))[0])
        #Indx1=~(np.isnan(x) | np.isinf(x))
        x1=x[Indx1]
        y1=y[Indx1]
        z1=z[Indx1]
    
        #Finding NaN and Inf in y dataset
        Indx2=np.int_((np.nonzero((np.isnan(y1)!=1) & (np.isinf(y1)!=1)))[0])
        #Indx2=~(np.isnan(y1) | np.isinf(1))
        x2=x1[Indx2]
        y2=y1[Indx2]
        z2=z1[Indx2]
    
        #Finding NaN and Inf in z dataset
        Indx3=np.int_((np.nonzero((np.isnan(z2)!=1) & (np.isinf(z2)!=1)))[0])
        #Indx3=~(np.isnan(z2) | np.isinf(z2))
        x3=x2[Indx3]
        y3=y2[Indx3]
        z3=z2[Indx3]
    
        #Datasets without NaN and Inf
        x=x3.copy()
        y=y3.copy()
        z=z3.copy()

    #Total number of elements in dataset after removal of NaN and Inf
    N=len(x)
    
    #Total number of NaN and Inf elements in dataset
    NNaNInf=NTotal-N

    #--------------------------------------------------------------------------
    #Calculating statistical properties for x, y and z data

    if ((np.size(x)>0) and (np.size(y)==0) and (np.size(z)==0)): #Only x is available
    
        xmean=np.mean(x) #Mean value of x data
        xmin=np.min(x) #Minimum value of x data
        xmax=np.max(x) #Maximum value of x data
        xstd=np.std(x,ddof=1) #Standard deviation of x data
        xSE=xstd/np.sqrt(len(x)) #Standard error of x data
        xCV=xstd/xmean #Coefficient of variation of x data
    
    elif ((np.size(x)>0) and (np.size(y)>0) and (np.size(z)==0)): #x and y are available

        xmean=np.mean(x) #Mean value of x data
        xmin=np.min(x) #Minimum value of x data
        xmax=np.max(x) #Maximum value of x data
        xstd=np.std(x,ddof=1) #Standard deviation of x data
        xSE=xstd/np.sqrt(len(x)) #Standard error of x data
        xCV=xstd/xmean #Coefficient of variation of x data
    
        ymean=np.mean(y) #Mean value of y data
        ymin=np.min(y) #Minimum value of y data
        ymax=np.max(y) #Maximum value of y data
        ystd=np.std(y,ddof=1) #Standard deviation of y data
        ySE=ystd/np.sqrt(len(y)) #Standard error of y data
        yCV=ystd/ymean #Coefficient of variation of y data
    
    elif ((np.size(x)>0) and (np.size(y)>0) and (np.size(z)>0)): #x, y and z are available
    
        xmean=np.mean(x) #Mean value of x data
        xmin=np.min(x) #Minimum value of x data
        xmax=np.max(x) #Maximum value of x data
        xstd=np.std(x,ddof=1) #Standard deviation of x data
        xSE=xstd/np.sqrt(len(x)) #Standard error of x data
        xCV=xstd/xmean #Coefficient of variation of x data
    
        ymean=np.mean(y) #Mean value of y data
        ymin=np.min(y) #Minimum value of y data
        ymax=np.max(y) #Maximum value of y data
        ystd=np.std(y,ddof=1) #Standard deviation of y data
        ySE=ystd/np.sqrt(len(y)) #Standard error of y data
        yCV=ystd/ymean #Coefficient of variation of y data
    
        zmean=np.mean(z) #Mean value of z data
        zmin=np.min(z) #Minimum value of z data
        zmax=np.max(z) #Maximum value of z data
        zstd=np.std(z,ddof=1) #Standard deviation of z data
        zSE=zstd/np.sqrt(len(z)) #Standard error of z data
        zCV=zstd/zmean #Coefficient of variation of z data
    
    #--------------------------------------------------------------------------
    #Calculating goodness of fit parameters (quantitavie evaluation of model) 

    if ((np.size(x)>0) and (np.size(y)>0) and (np.size(z)==0)): #x and y are available
    
        #x-y data
        #Pearson's correlation coefficient
        rxy=(np.sum((x-xmean)*(y-ymean)))/(np.sqrt((np.sum((x-xmean)**2))*(np.sum((y-ymean)**2)))) #Pearson's correlation coefficient 
        
        #Coefficient of determination
        R2xy=rxy**2 #Coefficient of determination

    elif ((np.size(x)>0) and (np.size(y)>0) and (np.size(z)>0)): #x, y and z are available
    
        #x-y data
        #Pearson's correlation coefficient
        rxy=(np.sum((x-xmean)*(y-ymean)))/(np.sqrt((np.sum((x-xmean)**2))*(np.sum((y-ymean)**2)))) #Pearson's correlation coefficient 
        
        #Coefficient of determination
        R2xy=rxy**2 #Coefficient of determination
    
        #x-z data
        #Pearson's correlation coefficient
        rxz=(np.sum((x-xmean)*(z-zmean)))/(np.sqrt((np.sum((x-xmean)**2))*(np.sum((z-zmean)**2)))) #Pearson's correlation coefficient 
    
        #Coefficient of determination
        R2xz=rxz**2 #Coefficient of determination
    
        #y-z data
        #Pearson's correlation coefficient
        ryz=(np.sum((y-ymean)*(z-zmean)))/(np.sqrt((np.sum((y-ymean)**2))*(np.sum((z-zmean)**2)))) #Pearson's correlation coefficient 
    
        #Coefficient of determination
        R2yz=ryz**2 #Coefficient of determination
    
    #--------------------------------------------------------------------------
    #Calculating probability density distribution

    if ((np.size(x)>0) and (np.size(y)==0) and (np.size(z)==0)): #Only x is available
    
        #Probability of x data
        binedgex=np.linspace(np.min(x),np.max(x),11) #Bin edges for x data
        binwidthx=np.diff(binedgex) #Bin width for x data, length(binwidth) is one element less than length(binedge) 
        bincenterx=binedgex[0:-1]+binwidthx/2 #Bin center for x data
    
        NPointsx,binedgeoutx=np.histogram(x,binedgex) #Calculating number of data points in each bin for x data
        histareax=NPointsx*binwidthx #Area under each element of histogram for x data
        fdensityx=NPointsx/(np.sum(histareax)) #Calculating probability density distribution for x data
        #Note: np.sum(fdensityx*binwidthx)=1
    
    elif ((np.size(x)>0) and (np.size(y)>0) and (np.size(z)==0)): #x and y are available
    
        #Probability of x data
        binedgex=np.linspace(np.min(x),np.max(x),11) #Bin edges for x data
        binwidthx=np.diff(binedgex) #Bin width for x data, length(binwidth) is one element less than length(binedge) 
        bincenterx=binedgex[0:-1]+binwidthx/2 #Bin center for x data
    
        NPointsx,binedgeoutx=np.histogram(x,binedgex) #Calculating number of data points in each bin for x data
        histareax=NPointsx*binwidthx #Area under each element of histogram for x data
        fdensityx=NPointsx/(np.sum(histareax)) #Calculating probability density distribution for x data
        #Note: np.sum(fdensityx*binwidthx)=1
    
        #Probability of y data
        binedgey=np.linspace(np.min(y),np.max(y),11) #Bin edges for y data
        binwidthy=np.diff(binedgey) #Bin width for y data, length(binwidth) is one element less than length(binedge) 
        bincentery=binedgey[0:-1]+binwidthy/2 #Bin center for y data
        
        NPointsy,binedgeouty=np.histogram(y,binedgey) #Calculating number of data points in each bin for y data
        histareay=NPointsy*binwidthy #Area under each element of histogram for y data
        fdensityy=NPointsy/(np.sum(histareay)) #Calculating probability density distribution for y data
        #Note: np.sum(fdensityy*binwidthy)=1
    
    elif ((np.size(x)>0) and (np.size(y)>0) and (np.size(z)>0)): #x, y and z are available
    
        #Probability of x data
        binedgex=np.linspace(np.min(x),np.max(x),11) #Bin edges for x data
        binwidthx=np.diff(binedgex) #Bin width for x data, length(binwidth) is one element less than length(binedge) 
        bincenterx=binedgex[0:-1]+binwidthx/2 #Bin center for x data
    
        NPointsx,binedgeoutx=np.histogram(x,binedgex) #Calculating number of data points in each bin for x data
        histareax=NPointsx*binwidthx #Area under each element of histogram for x data
        fdensityx=NPointsx/(np.sum(histareax)) #Calculating probability density distribution for x data
        #Note: np.sum(fdensityx*binwidthx)=1
    
        #Probability of y data
        binedgey=np.linspace(np.min(y),np.max(y),11) #Bin edges for y data
        binwidthy=np.diff(binedgey) #Bin width for y data, length(binwidth) is one element less than length(binedge) 
        bincentery=binedgey[0:-1]+binwidthy/2 #Bin center for y data
        
        NPointsy,binedgeouty=np.histogram(y,binedgey) #Calculating number of data points in each bin for y data
        histareay=NPointsy*binwidthy #Area under each element of histogram for y data
        fdensityy=NPointsy/(np.sum(histareay)) #Calculating probability density distribution for y data
        #Note: np.sum(fdensityy*binwidthy)=1
    
        #Probability of z data
        binedgez=np.linspace(np.min(z),np.max(z),11) #Bin edges for z data
        binwidthz=np.diff(binedgez) #Bin width for z data, length(binwidth) is one element less than length(binedge) 
        bincenterz=binedgez[0:-1]+binwidthz/2 #Bin center for z data
        
        NPointsz,binedgeoutz=np.histogram(z,binedgez) #Calculating number of data points in each bin for z data
        histareaz=NPointsz*binwidthz #Area under each element of histogram for z data
        fdensityz=NPointsz/(np.sum(histareaz)) #Calculating probability density distribution for z data
        #Note: np.sum(fdensityz*binwidthz)=1
    
    #--------------------------------------------------------------------------
    #Displaying results

    if ((np.size(x)>0) and (np.size(y)==0) and (np.size(z)==0)): #Only x is available
    
        print('--------------------------------------------------')
        val=[NTotal,N,NNaNInf]
        name=['N Input','N Used','N NaN Inf']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            print('     {0:15}     {1:.0f}'.format(name[i],val[i])) 
    
        print('--------------------------------------------------')
        val=[xmean,xmin,xmax,xstd,xSE,xCV]
        name=['x mean','x min','x max','x st. dev.','x st. err.','x coef. of var.']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            print('     {0:15}     {1:+0.10g}'.format(name[i],val[i])) 
    
        print('--------------------------------------------------')
    
        #Plotting
        plt.subplot(2,1,1)
        plt.plot(x)
    
        plt.xlabel('Data Points')
        plt.ylabel('x Data')
    
        plt.subplot(2,3,5)
        #plt.bar(bincenterx,fdensityx,width=0.8*binwidthx,label='x pdf')
        plt.plot(bincenterx,fdensityx,label='x pdf')
        plt.plot([xmean,xmean],[np.min(fdensityx),np.max(fdensityx)],label='x mean')
        plt.plot([xmean+xstd,xmean+xstd],[np.min(fdensityx),np.max(fdensityx)],label='mean+std')
        plt.plot([xmean-xstd,xmean-xstd],[np.min(fdensityx),np.max(fdensityx)],label='mean-std')
        plt.title('Probability Density Distribution, x Data')
        plt.ylabel('Probability Density, f')
        plt.xlabel('x Data')
        #plt.legend()
    
    elif ((np.size(x)>0) and (np.size(y)>0) and (np.size(z)==0)): #x and y are available
    
        print('--------------------------------------------------')
        val=[NTotal,N,NNaNInf]
        name=['N Input','N Used','N NaN Inf']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            print('     {0:15}     {1:.0f}'.format(name[i],val[i])) 
    
        print('--------------------------------------------------')
        val=[xmean,xmin,xmax,xstd,xSE,xCV]
        name=['x mean','x min','x max','x st. dev.','x st. err.','x coef. of var.']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            print('     {0:15}     {1:+0.10g}'.format(name[i],val[i])) 
    
        print('--------------------------------------------------')
        val=[ymean,ymin,ymax,ystd,ySE,yCV]
        name=['y mean','y min','y max','y st. dev.','y st. err.','y coef. of var.']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            print('     {0:15}     {1:+0.10g}'.format(name[i],val[i])) 
    
        print('--------------------------------------------------')
        val=[R2xy]
        name=['R2 x-y']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            print('     {0:15}     {1:+0.10g}'.format(name[i],val[i])) 
    
        print('--------------------------------------------------')
    
        #Plotting
        plt.subplot(1,3,1)
        plt.scatter(x,y)
    
        plt.xlabel('x Data')
        plt.ylabel('y Data')
    
        plt.subplot(1,3,2)
        #plt.bar(bincenterx,fdensityx,width=0.8*binwidthx,label='x pdf')
        plt.plot(bincenterx,fdensityx,label='x pdf')
        plt.plot([xmean,xmean],[np.min(fdensityx),np.max(fdensityx)],label='x mean')
        plt.plot([xmean+xstd,xmean+xstd],[np.min(fdensityx),np.max(fdensityx)],label='mean+std')
        plt.plot([xmean-xstd,xmean-xstd],[np.min(fdensityx),np.max(fdensityx)],label='mean-std')
        plt.title('Probability Density Distribution, x Data')
        plt.ylabel('Probability Density, f')
        plt.xlabel('x Data')
        #plt.legend()
    
        plt.subplot(1,3,3)
        #plt.bar(bincentery,fdensityy,width=0.8*binwidthy,label='y pdf')
        plt.plot(bincentery,fdensityy,label='y pdf')
        plt.plot([ymean,ymean],[np.min(fdensityy),np.max(fdensityy)],label='y mean')
        plt.plot([ymean+ystd,ymean+ystd],[np.min(fdensityy),np.max(fdensityy)],label='mean+std')
        plt.plot([ymean-ystd,ymean-ystd],[np.min(fdensityy),np.max(fdensityy)],label='mean-std')
        plt.title('Probability Density Distribution, y Data')
        plt.xlabel('y Data')
        plt.ylabel('Probability Density, f')
        #plt.legend()
    
    elif ((np.size(x)>0) and (np.size(y)>0) and (np.size(z)>0)): #x, y and z are available
    
        print('--------------------------------------------------')
        val=[NTotal,N,NNaNInf]
        name=['N Input','N Used','N NaN Inf']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            print('     {0:15}     {1:.0f}'.format(name[i],val[i])) 
    
        print('--------------------------------------------------')
        val=[xmean,xmin,xmax,xstd,xSE,xCV]
        name=['x mean','x min','x max','x st. dev.','x st. err.','x coef. of var.']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            print('     {0:15}     {1:+0.10g}'.format(name[i],val[i])) 
    
        print('--------------------------------------------------')
        val=[ymean,ymin,ymax,ystd,ySE,yCV]
        name=['y mean','y min','y max','y st. dev.','y st. err.','y coef. of var.']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            print('     {0:15}     {1:+0.10g}'.format(name[i],val[i])) 
    
        print('--------------------------------------------------')
        val=[zmean,zmin,zmax,zstd,zSE,zCV]
        name=['z mean','z min','z max','z st. dev.','z st. err.','z coef. of var.']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            print('     {0:15}     {1:+0.10g}'.format(name[i],val[i])) 
    
        print('--------------------------------------------------')
        val=[R2xy,R2xz,R2yz]
        name=['R2 x-y','R2 x-z','R2 y-z']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            print('     {0:15}     {1:+0.10g}'.format(name[i],val[i])) 
    
        print('--------------------------------------------------')
    
        #Plotting
        fig = plt.figure()
        plt.subplot(2,3,1)
        plt.scatter(x,z)
    
        plt.xlabel('x Data')
        plt.ylabel('z Data')
    
        plt.subplot(2,3,2)
        plt.scatter(y,z)
    
        plt.xlabel('y Data')
        plt.ylabel('z Data')
    
        #plt.subplot(2,3,3)
        ax = fig.add_subplot(2,3,3, projection='3d')
        ax.scatter(x,y,z)
    
        ax.set_xlabel('x Data')
        ax.set_ylabel('y Data')
        ax.set_zlabel('z Data')
    
        plt.subplot(2,3,4)
        #plt.bar(bincenterx,fdensityx,width=0.8*binwidthx,label='x pdf')
        plt.plot(bincenterx,fdensityx,label='x pdf')
        plt.plot([xmean,xmean],[np.min(fdensityx),np.max(fdensityx)],label='x mean')
        plt.plot([xmean+xstd,xmean+xstd],[np.min(fdensityx),np.max(fdensityx)],label='mean+std')
        plt.plot([xmean-xstd,xmean-xstd],[np.min(fdensityx),np.max(fdensityx)],label='mean-std')
        plt.title('Probability Density Distribution, x Data')
        plt.ylabel('Probability Density, f')
        plt.xlabel('x Data')
        #plt.legend()
    
        plt.subplot(2,3,5)
        #plt.bar(bincentery,fdensityy,width=0.8*binwidthy,label='y pdf')
        plt.plot(bincentery,fdensityy,label='y pdf')
        plt.plot([ymean,ymean],[np.min(fdensityy),np.max(fdensityy)],label='y mean')
        plt.plot([ymean+ystd,ymean+ystd],[np.min(fdensityy),np.max(fdensityy)],label='mean+std')
        plt.plot([ymean-ystd,ymean-ystd],[np.min(fdensityy),np.max(fdensityy)],label='mean-std')
        plt.title('Probability Density Distribution, y Data')
        plt.xlabel('y Data')
        plt.ylabel('Probability Density, f')
        #plt.legend()
    
        plt.subplot(2,3,6)
        #plt.bar(bincenterz,fdensityz,width=0.8*binwidthz,label='z pdf')
        plt.plot(bincenterz,fdensityz,label='z pdf')
        plt.plot([zmean,zmean],[np.min(fdensityz),np.max(fdensityz)],label='z mean')
        plt.plot([zmean+zstd,zmean+zstd],[np.min(fdensityz),np.max(fdensityz)],label='mean+std')
        plt.plot([zmean-zstd,zmean-zstd],[np.min(fdensityz),np.max(fdensityz)],label='mean-std')
        plt.title('Probability Density Distribution, z Data')
        plt.xlabel('z Data')
        plt.ylabel('Probability Density, f')
        #plt.legend()
    
    #--------------------------------------------------------------------------
