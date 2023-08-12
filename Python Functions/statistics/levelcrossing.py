def levelcrossing(x, CrossingLevel=None, fs=1, dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2017-04-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.levelcrossing
    =========================

    .. code:: python

        UpCrossIndx, UpCrossTime, UpCrossValue, CrestIndx, CrestTime, CrestValue, TroughIndx, TroughTime, TroughValue, t = scientimate.levelcrossing(x, CrossingLevel=None, fs=1, dispout='no')

    Description
    -----------

    Calculate crossing point for a given level by using an upward zero crossing method

    Inputs
    ------

    x
        Input time series (oscillatory) data
    CrossingLevel=mean(x)
        | Level that crossing reported for
        | If x is deterended wave (mean=0), then set CrossingLevel=0
    fs=1
        | Sampling frequency that data collected at in (Hz)
        | If fs is not given, then default fs is fs=1;
        | If fs=1, then index of data points represents time as well
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    UpCrossIndx
        Index of up-crossing location
    UpCrossTime
        Time of up-crossing location (s)
    UpCrossValue
        Value of up-crossing
    CrestIndx
        Index of crest location
    CrestTime
        Time of wave crest location (s)
    CrestValue
        Value of wave crest
    TroughIndx
        Index of trough location
    TroughTime
        Time of wave trough location (s)
    TroughValue
        Value of wave trough
    t
        Time (s)

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np
        import scipy as sp
        from scipy import signal

        fs=2 #Sampling frequency
        duration=1024 #Duration of the data
        N=fs*duration #Total number of points
        df=fs/N #Frequency difference 
        dt=1/fs #Time difference, dt=1/fs
        t=np.linspace(0,duration-dt,N) #Time
        x=sp.signal.detrend(0.5*np.cos(2*np.pi*0.2*t)+(-0.1+(0.1-(-0.1)))*np.random.rand(N))
        UpCrossIndx,UpCrossTime,UpCrossValue,CrestIndx,CrestTime,CrestValue,TroughIndx,TroughTime,TroughValue,t=sm.levelcrossing(x,np.mean(x),fs,'yes')

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
    #Assign default values

    if CrossingLevel is None: CrossingLevel=np.mean(x)

    #--------------------------------------------------------------------------
    #Checking if a crossing level has any intersection with input data

    if (CrossingLevel>=np.max(x)) or (CrossingLevel<=np.min(x)):
        import sys
        sys.exit('No intersection between input data and given level, try another CrossingLevel value')

    #--------------------------------------------------------------------------

    #Deterending input data
    xDetrended=x-CrossingLevel

    #--------------------------------------------------------------------------

    #Generating time vector
    N=len(x) #Length of a time series, total number of points in input file, N=fs*duration
    duration=N/fs #Duration time that data are collected (second), duration=N/fs
    dt=1/fs #Time difference between consecutive samples, dt=1/fs=duration/N
    # t(:,1)=linspace(dt,N/fs,N); #Time from 0 to T-dt, equally spaced at dt
    t=np.arange(0,duration,dt) #Time from 0 to T-dt, equally spaced at dt

    #--------------------------------------------------------------------------
    #Detecting first and last wave (first and last complete crest-trough cycle)

    #Detecting the start point of the first wave (first complete crest-trough cycle)
    if ((xDetrended[0]==0) and (xDetrended[1]>0)): 
        len3=1

    if (((xDetrended[0]==0) and (xDetrended[1]<0)) or (xDetrended[0]!=0)):
        for i in range(0,N-1,1):
            if ((xDetrended[i]<0) and (xDetrended[i+1]>0)): 
                len3=i
                break

    #Detecting the end point of last wave (last complete crest-trough cycle)
    if ((xDetrended[-1]==0) and (xDetrended[-2]<0)):
        len4=N

    if (((xDetrended[-1]==0) and (xDetrended[-2]>0)) or (xDetrended[-1]!=0)):
        for i in range(N-2,0,-1):
            if ((xDetrended[i]<0) and (xDetrended[i+1]>0)):
                len4=i
                break

    #--------------------------------------------------------------------------
    #Detecting zero crossing points

    m=-1
    n=-1

    xupcross=np.array([])
    yupcross=np.array([])
    xupcrossIndx=np.array([])
    
    xdowncross=np.array([])
    ydowncross=np.array([])
    xdowncrossIndx=np.array([])

    for i in range(len3,len4+1,1): 
        
        #Detecting up-crossing zero-crossing points
        if i==1:
            m=m+1
            xupcross=np.append(xupcross,dt)
            yupcross=np.append(yupcross,0)
            xupcrossIndx=np.append(xupcrossIndx,i)
        
        if ((i>1) and (i<N)):
            if ((xDetrended[i]<0) and (xDetrended[i+1]>0)):
                m=m+1
                xupcross=np.append(xupcross,t[i]-(t[i+1]-t[i])/(xDetrended[i+1]-xDetrended[i])*xDetrended[i])
                yupcross=np.append(yupcross,0)
                xupcrossIndx=np.append(xupcrossIndx,i)
        
        if i==N:
            m=m+1
            xupcross=np.append(xupcross,t[-1])
            yupcross=np.append(yupcross,0)
            xupcrossIndx=np.append(xupcrossIndx,i)
        
        #Detecting down-crossing zero-crossing points
        if ((i>1) and (i<N)):
            if ((xDetrended[i]>0) and (xDetrended[i+1]<0)):
                n=n+1
                xdowncross=np.append(xdowncross,t[i]-(t[i+1]-t[i])/(xDetrended[i+1]-xDetrended[i])*xDetrended[i])
                ydowncross=np.append(ydowncross,0)
                xdowncrossIndx=np.append(xdowncrossIndx,i)
                if xdowncrossIndx[n]>=N:
                    xdowncrossIndx[n]=N-1

    #Converting to int
    xupcrossIndx=np.int64(xupcrossIndx)
    xdowncrossIndx=np.int64(xdowncrossIndx)

    #--------------------------------------------------------------------------
    #Detecting crest and trough

    m=-1
    n=-1

    xmax=np.array([])
    ymax=np.array([])
    xmaxIndx=np.array([])
    
    xmin=np.array([])
    ymin=np.array([])
    xminIndx=np.array([])

    for i in range(xupcrossIndx[0],xupcrossIndx[-1]+1,1):
        
        #Detecting crest
        if ((i>1) and (i<N)):
            
            if ((xDetrended[i]>xDetrended[i-1]) and (xDetrended[i]>xDetrended[i+1])):
                
                #Check if after last crest, a trough is detected or not 
                #m==n means it is detected and the new y(i,1) is the crest of the next wave
                if xDetrended[i]>0:
                    
                    if m==n:
                        
                        m=m+1
                        xmax=np.append(xmax,t[i])
                        ymax=np.append(ymax,xDetrended[i])
                        xmaxIndx=np.append(xmaxIndx,i)

                    elif m!=0:
                            
                        #Replacingthe old crest location with new one if the new one is larger
                        if ((xDetrended[i]>ymax[m]) and (m!=0)):
                            
                            xmax=np.append(xmax,t[i])
                            ymax=np.append(ymax,xDetrended[i])
                            xmaxIndx=np.append(xmaxIndx,i)
                                
        #Detecting trough
        if ((i>1) and (i<N)):
            
            if ((xDetrended[i]<xDetrended[i-1]) and (xDetrended[i]<xDetrended[i+1])):
                
                if xDetrended[i]<0:
                    
                    #Check if after last trough, a crest is detected or not 
                    #n==m-1 means it is detected and the new y(i,1) is the trough of next wave
                    if n==m-1:
                        
                        n=n+1
                        xmin=np.append(xmin,t[i])
                        ymin=np.append(ymin,xDetrended[i])
                        xminIndx=np.append(xminIndx,i)
                        
                    elif n!=0:
                            
                        #Replacingthe old crest location with new one if the new one is smaller
                        if xDetrended[i]<ymin[n]:
                            
                            xmin=np.append(xmin,t[i])
                            ymin=np.append(ymin,xDetrended[i])
                            xminIndx=np.append(xminIndx,i)
                                
    #x=xDetrended.copy() #Water surface elevation time series

    #Converting to int
    xmaxIndx=np.int64(xmaxIndx)
    xminIndx=np.int64(xminIndx)

    #--------------------------------------------------------------------------
    #Assigning output

    UpCrossIndx=xupcrossIndx.copy() #Up crossing Index
    UpCrossTime=xupcross.copy() #Up crossing time
    UpCrossValue=yupcross+CrossingLevel #Up crossing value
    
    CrestIndx=xmaxIndx.copy() #Crest index
    CrestTime=xmax.copy() #Crest time
    CrestValue=ymax+CrossingLevel #Crest value
    
    TroughIndx=xminIndx.copy() #Trough index
    TroughTime=xmin.copy() #Trough time
    TroughValue=ymin+CrossingLevel #Trough value

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
    
        #Plotting data time series
        plt.plot(t,x,label='Time Series')
        plt.plot([t[0],t[-1]],[CrossingLevel,CrossingLevel],label='Crossing Level')
        
        plt.scatter(xmax,ymax+CrossingLevel,marker='*',label='Crest')
        plt.scatter(xmin,ymin+CrossingLevel,marker='o',label='Trough')
        
        plt.scatter(xupcross,yupcross+CrossingLevel,marker='^',label='Up Crossing')
        plt.scatter(xdowncross,ydowncross+CrossingLevel,marker='v',label='Down Crossing')
    
        plt.xlabel('Time (s)')
        plt.ylabel('Wave')
        plt.legend()
    
    #--------------------------------------------------------------------------
    #Outputs
    return UpCrossIndx, UpCrossTime, UpCrossValue, CrestIndx, CrestTime, CrestValue, TroughIndx, TroughTime, TroughValue, t

    #--------------------------------------------------------------------------
