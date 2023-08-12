def replaceoutlier(x, WindowSize=15, zscore_threshold=2, interpMethod='linear', dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2020-05-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    scientimate.replaceoutlier
    ==========================
    
    .. code:: python
    
        xReplaced, outlier_Indx = scientimate.replaceoutlier(x, WindowSize=15, zscore_threshold=2, interpMethod='linear', dispout='no')
    
    Description
    -----------
    
    Remove outliers in the time series using moving z-score window
    
    Inputs
    ------
    
    x
        Input data
    WindowSize=15
        Window size (number of adjacent elements) that is used for moving window, should be equal or larger than 3
    zscore_threshold=2
        | z-score threshold to define outliers
        | data in range of x < (xmean-zscore_threshold*std) or x > (xmean+zscore_threshold*std) considered outliers
    interpMethod='linear'
        | Interpolation method for replacing spike points:
        | Matlab/Octave: 'linear', 'nearest', 'next', 'previous', 'pchip', 'cubic', 'spline'
        | Python/Scipy : 'linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic'
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)
    
    Outputs
    -------
    
    xReplaced
        Replaced data
    outlier_Indx
        Logical index of replaced points
    
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
        xReplaced,outlier_Indx=sm.replaceoutlier(x,37,2,'linear','yes')
    
        fs=2
        t=np.linspace(0,1023.5,1024*fs)
        x=sp.signal.detrend(0.5*np.cos(2*np.pi*0.2*t)+(-0.1+(0.1-(-0.1)))*np.random.rand(1024*fs))
        spikeloc=np.arange(10,len(t),100)
        x[np.int64(spikeloc+np.round(2*np.random.randn(len(spikeloc))))]=np.sign(np.random.randn(len(spikeloc)))
        xReplaced,outlier_Indx=sm.replaceoutlier(x,21,2,'linear','yes')
    
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
    from scipy import interpolate
    #import pandas as pd
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
    
    #--------------------------------------------------------------------------
    #Removing outliers
    
    #Assigning empty to NaN_Indx
    NaN_Indx=[]
    
    Nw_half=int(Nw/2) #Half of the range width

    #Moving mean
    rolled_mean=np.zeros(N)

    for i in range(0,N,1):
        if ((i>=0) and (i<=Nw_half-1)):
            rolled_mean[i]=np.mean(x[0:i+Nw_half+1])
            
        elif ((i>=Nw_half) and (i<=N-Nw_half-1)):
            rolled_mean[i]=np.mean(x[i-Nw_half:i+Nw_half+1])
            
        elif ((i>=N-Nw_half) and (i<=N)):
            rolled_mean[i]=np.mean(x[i-Nw_half:])

    #Or
    #import pandas as pd
    #df=pd.DataFrame(data=xInput)
    #rolled_mean=df.rolling(Nw,min_periods=1,center=True).mean()

    #Moving std
    rolled_std=np.zeros(N)

    for i in range(0,N,1):
        if ((i>=0) and (i<=Nw_half-1)):
            rolled_std[i]=np.std(x[0:i+Nw_half+1])
            
        elif ((i>=Nw_half) and (i<=N-Nw_half-1)):
            rolled_std[i]=np.std(x[i-Nw_half:i+Nw_half+1])
            
        elif ((i>=N-Nw_half) and (i<=N)):
            rolled_std[i]=np.std(x[i-Nw_half:])

    #Or
    #import pandas as pd
    #df=pd.DataFrame(data=xInput)
    #rolled_std=df.rolling(Nw,min_periods=1,center=True).std()
    
    #Calculate moving z-score
    rolled_z_score=(x-rolled_mean)/rolled_std
    #Or rolled_z_score=(df-rolled_mean)/rolled_std
    
    #Replace outliers with NaN
    xInput[np.absolute(rolled_z_score)>zscore_threshold]=None
    #Or xInput[np.absolute(rolled_z_score.iloc[:,0])>zscore_threshold]=None

    #Find NaN
    NaN_Indx=(np.isnan(xInput)==1)
    Valid_Indx=(np.isnan(xInput)!=1)
    
    #Assigning x to xReplaced
    xReplaced=x.copy()
    
    #Replacing outliers
    samples=np.arange(1,len(x)+1,1)
    #print(samples[NaN_Indx==True].shape)
    if True in NaN_Indx:
        #xReplaced[NaN_Indx]=np.interp(samples[NaN_Indx==True],samples[Valid_Indx==True],xReplaced[Valid_Indx==True])
        fun=sp.interpolate.interp1d(samples[Valid_Indx==True],xReplaced[Valid_Indx==True],kind=interpMethod,bounds_error=False,fill_value='extrapolate')
        xReplaced[NaN_Indx==True]=fun(samples[NaN_Indx==True])
    
    #Assigning NaN_Indx to outlier_Indx
    outlier_Indx=NaN_Indx.copy()
    
    #--------------------------------------------------------------------------
    #Defining replaced points
    
    if True in NaN_Indx:
        NumberReplacedPoint=np.sum(NaN_Indx)
        PercentReplacedPoint=np.sum(NaN_Indx)/np.size(x)*100
    else:
        NumberReplacedPoint=0
        PercentReplacedPoint=0
    
    #--------------------------------------------------------------------------
    #Displaying results
    
    if dispout=='yes':

        print('Total number of pints replaced = ',str(NumberReplacedPoint))
        print('Percent of replaced points     = ',str(round(PercentReplacedPoint,2)),' %')
        
        plt.subplot(2,1,1)
        plt.plot(samples,x)
        plt.xlabel('Sample')
        plt.ylabel('Data')
        plt.title('Input Data')
    
        plt.subplot(2,1,2)
        plt.plot(samples,xReplaced)
        plt.scatter(samples[NaN_Indx],xReplaced[NaN_Indx])
        plt.xlabel('Sample')
        plt.ylabel('Data')
        plt.title('Replaced Data')
    
    #--------------------------------------------------------------------------
    #Outputs
    return xReplaced, outlier_Indx

    #--------------------------------------------------------------------------
    