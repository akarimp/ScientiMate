def replacespikediff(x, WindowSize=5, spikedefineMethod='ellipse', nrpeat=1, DespikeScaleFactor=1, interpMethod='linear', dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2017-02-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.replacespikediff
    ============================

    .. code:: python

        xDespiked, Indx = scientimate.replacespikediff(x, WindowSize=5, spikedefineMethod='ellipse', nrpeat=1, DespikeScaleFactor=1, interpMethod='linear', dispout='yes')

    Description
    -----------

    Remove spikes in the time series using a local difference of data respect to a moving average window

    Inputs
    ------

    x
        Input data
    WindowSize=5
        Window size (number of adjacent elements) that is used for smoothing, should be equal or larger than 3
    spikedefineMethod='ellipse'
        | Method to define spike points
        | 'ellipse': use both local difference and its gradient
        | 'circle': use local difference only
    nrpeat=1
        Number of time despiking procedure is repeating
    DespikeScaleFactor=1
        Scaling a despiking threshold
    interpMethod='linear'
        | Interpolation method for replacing spike points:
        | Matlab/Octave: 'linear', 'nearest', 'next', 'previous', 'pchip', 'cubic', 'spline'
        | Python/Scipy : 'linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic'
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    xDespiked
        Dispiked data
    Indx
        Index of despiked points

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np
        import scipy as sp
        from scipy import signal

        fs=128
        t=np.linspace(0,9.5,10*fs)
        x=np.sin(2*np.pi*0.3*t)+0.1*np.sin(2*np.pi*4*t)
        spikeloc=np.arange(10,len(t),100)
        x[np.int64(spikeloc+np.round(2*np.random.randn(len(spikeloc))))]=np.sign(np.random.randn(len(spikeloc)))
        x[220:225]=1.5
        x=x+5
        xDespiked,Indx=sm.replacespikediff(x,5,'ellipse',2,1,'linear','yes')

        fs=2
        t=np.linspace(0,1023.5,1024*fs)
        x=sp.signal.detrend(0.5*np.cos(2*np.pi*0.2*t)+(-0.1+(0.1-(-0.1)))*np.random.rand(1024*fs))
        spikeloc=np.arange(10,len(t),100)
        x[np.int64(spikeloc+np.round(2*np.random.randn(len(spikeloc))))]=np.sign(np.random.randn(len(spikeloc)))
        xDespiked,Indx=sm.replacespikediff(x,5,'ellipse',1,1,'linear','yes')

    References
    ----------

    Goring, D. G., & Nikora, V. I. (2002). 
    Despiking acoustic Doppler velocimeter data. 
    Journal of Hydraulic Engineering, 128(1), 117-126.

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

    #--------------------------------------------------------------------------

    #Preserving an input data
    xInput=x.copy()

    #--------------------------------------------------------------------------

    N=len(x) #Length of the input data
    Nw=WindowSize #Window size
    
    #Checking size of the window
    if Nw>N: Nw=N
    if Nw<=3: Nw=3
    
    #Creating window function in time domain
    
    nwin=np.arange(0,Nw,1) #Counter from 0 to Nw-1
    
    #Calculating window function
    WnRec=np.ones(Nw) #Rectangular window function (moving average filter)
    Wn=WnRec
    WnNorm=Wn/np.sum(Wn) #Normalizing a window function in time domain

    #--------------------------------------------------------------------------
    #Removing spikes

    for i in range(0,nrpeat,1):
        
        #Smoothing input using convolution
        xSmoothed=np.convolve(x,WnNorm,'same') #Smoothing input using convolution
    
        #Removing alias from edge of smooted data
        for n in range (0,int(np.floor(Nw/2)),1):
            WnEdge=WnNorm[(int(np.ceil(Nw/2))-1)-(n):(int(np.ceil(Nw/2))-1)+(n)+1]
            WnEdgeNorm=WnEdge/np.sum(WnEdge) #Normalizing a window function in time domain
            xSmoothed[n]=np.sum(x[0:n+(n)+1]*WnEdgeNorm)
            xSmoothed[N-n-1]=np.sum(x[(N-n-1)-(n):(N+1)]*WnEdgeNorm)
    
        
        #Local difference
        xDiff=np.abs(x-xSmoothed) #Local difference of x respecte to moving average window
    
        #Standard deviation
        xStd=np.std(x) #Standard deviation of x
        xDiffStd=np.std(xDiff) #Standard deviation of xDiff
    
        #Expected absolute maximum based on normal distributaion
        Lambda=np.sqrt(2*np.log(N))
    
        #Axis of ellipse
        a1=Lambda*xStd #Major axis of ellipse for xDiff vs x
        b1=Lambda*xDiffStd #Minor axis of ellipse for xDiff vs x
    
        #Calculating mean values
        xmean=np.mean(x)
        xDiffmean=np.mean(xDiff)
    
        #Calculating distances based on ellipse equation
        distellipse=((x-xmean)**2/a1**2+(xDiff-xDiffmean)**2/b1**2)
    
        #Defining spike points 
        if spikedefineMethod=='ellipse':
            #If a point is locatated outside an ellipse, it is considered as spike
            Indx1=np.int_((np.nonzero(distellipse>1*DespikeScaleFactor))[0])
        elif spikedefineMethod=='circle':
            #If a point is locatated outside a circle, it is considered as spike
            #Considering a Normal distributaion, all points are smaller than 3 times of standard deviation
            Indx1=np.int_((np.nonzero(xDiff>(xDiffmean+3*xDiffStd)*DespikeScaleFactor))[0])
    
        
        #Setting spike points as NaN
        xDespiked=x.copy()
        xDespiked[Indx1]=np.nan
    
        #Locating all spike pints
        Indx2=np.int_((np.nonzero(np.isnan(xDespiked)==1))[0])
    
        #Locating all non-spike pints
        Indx3=np.int_((np.nonzero(np.isnan(xDespiked)!=1))[0])
    
        #Replacing spike points
        samples=np.arange(1,len(x)+1,1)
        if len(Indx2)!=0:
            #xDespiked[Indx2]=np.interp1(samples[Indx2],samples[Indx3],xDespiked[Indx3],interpMethod,'extrap')
            fun=sp.interpolate.interp1d(samples[Indx3],xDespiked[Indx3],kind=interpMethod,fill_value='extrapolated')
            xDespiked[Indx2]=fun(samples[Indx2])
    
        #Assigning a despiked data to x for next repeat
        x=xDespiked.copy()


    #--------------------------------------------------------------------------
    #Defining despiked points

    Indx=np.int_((np.nonzero(xInput!=xDespiked))[0])
    NumberDespikedPoint=len(Indx)
    PercentDespikedPoint=len(Indx)/len(x)*100

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':

        print('Total number of pints despiked = ',str(NumberDespikedPoint))
        print('Percent of despiked points     = ',str(round(PercentDespikedPoint,2)),' %')
        
        plt.subplot(2,1,1)
        plt.plot(samples,xInput)
        #plt.hold='True'
        plt.scatter(samples[Indx],xInput[Indx])
        plt.xlabel('Sample')
        plt.ylabel('Data')
        plt.title('Input Data')
    
        plt.subplot(2,1,2)
        plt.plot(samples,xDespiked)
        #plt.hold='True'
        plt.scatter(samples[Indx],xDespiked[Indx])
        plt.xlabel('Sample')
        plt.ylabel('Data')
        plt.title('Despiked Data')
    

    #--------------------------------------------------------------------------
    #Outputs
    return xDespiked, Indx

    #--------------------------------------------------------------------------
