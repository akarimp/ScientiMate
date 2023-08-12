def movingwindow(x, WindowSize=3, StatisticalMethod='mean', dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2020-08-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.movingwindow
    ========================

    .. code:: python

        moving_statistics = scientimate.movingwindow(x, WindowSize=3, StatisticalMethod='mean', dispout='no')

    Description
    -----------

    Calculate statistics of moving window through 1-d x data

    Inputs
    ------

    x
        Input data
    WindowSize=3
        | Window size (number of adjacent elements) that is used for moving window, should be equal or larger than 3
        | Window size should be an odd integer
    StatisticalMethod='mean'
        | Statistical value of moving window to be reported:
        | 'mean': Moving mean
        | 'std': Moving standard deviation
        | 'min': Moving minimum
        | 'max': Moving maximum
        | 'sum': Moving sum
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    moving_statistics
        Statistical value of moving window

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        fs=128
        t=np.linspace(0,9.5,10*fs)
        x=np.sin(2*np.pi*0.3*t)+0.1*np.sin(2*np.pi*4*t)
        moving_statistics=sm.movingwindow(x,37,'mean','yes')

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

    N=len(x) #Length of the input data
    Nw=WindowSize #Window size
    
    #Checking size of the window
    if Nw>N: Nw=N
    if Nw<=3: Nw=3

    #--------------------------------------------------------------------------

    Nw_half=int(Nw/2) #Half of the range width

    moving_statistics=np.zeros(N)

    #Calculate moving mean
    if StatisticalMethod=='mean':

        for i in range(0,N,1):
            if ((i>=0) and (i<=Nw_half-1)):
                moving_statistics[i]=np.mean(x[0:i+Nw_half+1])
                
            elif ((i>=Nw_half) and (i<=N-Nw_half-1)):
                moving_statistics[i]=np.mean(x[i-Nw_half:i+Nw_half+1])
                
            elif ((i>=N-Nw_half) and (i<=N)):
                moving_statistics[i]=np.mean(x[i-Nw_half:])
                
        #Or
        #import pandas as pd
        #df=pd.DataFrame(data=x)
        #moving_statistics=df.rolling(Nw,min_periods=1,center=True).mean()

    #Calculate moving standard deviation
    elif StatisticalMethod=='std':

        for i in range(0,N,1):
            if ((i>=0) and (i<=Nw_half-1)):
                moving_statistics[i]=np.std(x[0:i+Nw_half+1])
                
            elif ((i>=Nw_half) and (i<=N-Nw_half-1)):
                moving_statistics[i]=np.std(x[i-Nw_half:i+Nw_half+1])
                
            elif ((i>=N-Nw_half) and (i<=N)):
                moving_statistics[i]=np.std(x[i-Nw_half:])
                
        #Or
        #import pandas as pd
        #df=pd.DataFrame(data=x)
        #moving_statistics=df.rolling(Nw,min_periods=1,center=True).std()

    #Calculate moving minimum
    elif StatisticalMethod=='min':

        for i in range(0,N,1):
            if ((i>=0) and (i<=Nw_half-1)):
                moving_statistics[i]=np.min(x[0:i+Nw_half+1])
                
            elif ((i>=Nw_half) and (i<=N-Nw_half-1)):
                moving_statistics[i]=np.min(x[i-Nw_half:i+Nw_half+1])
                
            elif ((i>=N-Nw_half) and (i<=N)):
                moving_statistics[i]=np.min(x[i-Nw_half:])

        #Or
        #import pandas as pd
        #df=pd.DataFrame(data=x)
        #moving_statistics=df.rolling(Nw,min_periods=1,center=True).min()

    #Calculate moving maximum
    elif StatisticalMethod=='max':

        for i in range(0,N,1):
            if ((i>=0) and (i<=Nw_half-1)):
                moving_statistics[i]=np.max(x[0:i+Nw_half+1])
                
            elif ((i>=Nw_half) and (i<=N-Nw_half-1)):
                moving_statistics[i]=np.max(x[i-Nw_half:i+Nw_half+1])
                
            elif ((i>=N-Nw_half) and (i<=N)):
                moving_statistics[i]=np.max(x[i-Nw_half:])

        #Or
        #import pandas as pd
        #df=pd.DataFrame(data=x)
        #moving_statistics=df.rolling(Nw,min_periods=1,center=True).max()

    #Calculate moving sum
    elif StatisticalMethod=='sum':

        for i in range(0,N,1):
            if ((i>=0) and (i<=Nw_half-1)):
                moving_statistics[i]=np.sum(x[0:i+Nw_half+1])
                
            elif ((i>=Nw_half) and (i<=N-Nw_half-1)):
                moving_statistics[i]=np.sum(x[i-Nw_half:i+Nw_half+1])
                
            elif ((i>=N-Nw_half) and (i<=N)):
                moving_statistics[i]=np.sum(x[i-Nw_half:])

        #Or
        #import pandas as pd
        #df=pd.DataFrame(data=x)
        #moving_statistics=df.rolling(Nw,min_periods=1,center=True).sum()


    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
        
        plt.plot(x, label='Input')
        plt.plot(moving_statistics, label='Moving Statistics')
        
        plt.xlabel('Sample')
        plt.ylabel('Data')
        plt.legend()
        

    #--------------------------------------------------------------------------
    #Outputs
    return moving_statistics

    #--------------------------------------------------------------------------
