def wavefromsurfaceelevzcross(Eta, fs, dispout='no'):
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

    scientimate.wavefromsurfaceelevzcross
    =====================================

    .. code:: python

        Hs, Ts, Hz, Tz, Hrms, H, T = scientimate.wavefromsurfaceelevzcross(Eta, fs, dispout='no')

    Description
    -----------

    Calculate wave properties from water surface elevation by using an upward zero crossing method

    Inputs
    ------

    Eta
        Water surface elevation time series data in (m)
    fs
        Sampling frequency that data collected at in (Hz)
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    Hs
        Significant Wave Height (m)
    Ts
        Significant Wave Period (second)
    Hz
        Zero Crossing Mean Wave Height (m)
    Tz
        Zero Crossing Mean Wave Period (second)
    Hrms
        Root Mean Square Wave Height (m)
    H
        Wave Height Data Series array (m)
    T
        Wave Period Data Series array (second)

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
        Eta=sp.signal.detrend(0.5*np.cos(2*np.pi*0.2*t)+(-0.1+(0.1-(-0.1)))*np.random.rand(N))
        Hs,Ts,Hz,Tz,Hrms,H,T=sm.wavefromsurfaceelevzcross(Eta,fs,'yes')

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
    from scipy import signal
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
    
    Eta=type2numpy(Eta)

    #--------------------------------------------------------------------------

    #Deterending input data
    EtaDetrended=sp.signal.detrend(Eta)

    #--------------------------------------------------------------------------

    #Generating time vector
    N=len(Eta) #Length of a time series, total number of points in input file, N=fs*duration
    duration=N/fs #Duration time that data are collected (second), duration=N/fs
    dt=1/fs #Time difference between consecutive samples, dt=1/fs=duration/N
    # t(:,1)=linspace(dt,N/fs,N); #Time from 0 to T-dt, equally spaced at dt
    t=np.arange(0,duration,dt) #Time from 0 to T-dt, equally spaced at dt

    #--------------------------------------------------------------------------
    #Detecting first and last wave (first and last complete crest-trough cycle)

    #Detecting the start point of the first wave (first complete crest-trough cycle)
    if ((EtaDetrended[0]==0) and (EtaDetrended[1]>0)): 
        len3=1

    if (((EtaDetrended[0]==0) and (EtaDetrended[1]<0)) or (EtaDetrended[0]!=0)):
        for i in range(0,N-1,1):
            if ((EtaDetrended[i]<0) and (EtaDetrended[i+1]>0)): 
                len3=i
                break

    #Detecting the end point of last wave (last complete crest-trough cycle)
    if ((EtaDetrended[-1]==0) and (EtaDetrended[-2]<0)):
        len4=N

    if (((EtaDetrended[-1]==0) and (EtaDetrended[-2]>0)) or (EtaDetrended[-1]!=0)):
        for i in range(N-2,0,-1):
            if ((EtaDetrended[i]<0) and (EtaDetrended[i+1]>0)):
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
            if ((EtaDetrended[i]<0) and (EtaDetrended[i+1]>0)):
                m=m+1
                xupcross=np.append(xupcross,t[i]-(t[i+1]-t[i])/(EtaDetrended[i+1]-EtaDetrended[i])*EtaDetrended[i])
                yupcross=np.append(yupcross,0)
                xupcrossIndx=np.append(xupcrossIndx,i)
        
        if i==N:
            m=m+1
            xupcross=np.append(xupcross,t[-1])
            yupcross=np.append(yupcross,0)
            xupcrossIndx=np.append(xupcrossIndx,i)
        
        #Detecting down-crossing zero-crossing points
        if ((i>1) and (i<N)):
            if ((EtaDetrended[i]>0) and (EtaDetrended[i+1]<0)):
                n=n+1
                xdowncross=np.append(xdowncross,t[i]-(t[i+1]-t[i])/(EtaDetrended[i+1]-EtaDetrended[i])*EtaDetrended[i])
                ydowncross=np.append(ydowncross,0)
                xdowncrossIndx=np.append(xdowncrossIndx,i)
                if xdowncrossIndx[n]>=N:
                    xdowncrossIndx[n]=N-1

    #Converting to int
    xupcrossIndx=np.int_(xupcrossIndx)
    xdowncrossIndx=np.int_(xdowncrossIndx)

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
            
            if ((EtaDetrended[i]>EtaDetrended[i-1]) and (EtaDetrended[i]>EtaDetrended[i+1])):
                
                #Check if after last crest, a trough is detected or not 
                #m==n means it is detected and the new y(i,1) is the crest of the next wave
                if EtaDetrended[i]>0:
                    
                    if m==n:
                        
                        m=m+1
                        xmax=np.append(xmax,t[i])
                        ymax=np.append(ymax,EtaDetrended[i])
                        xmaxIndx=np.append(xmaxIndx,i)

                    elif m!=-1:
                            
                        #Replacingthe old crest location with new one if the new one is larger
                        if ((EtaDetrended[i]>ymax[m]) and (m!=-1)):
                            
                            xmax[m]=t[i]
                            ymax[m]=EtaDetrended[i]
                            xmaxIndx[m]=i
                                
        #Detecting trough
        if ((i>1) and (i<N)):
            
            if ((EtaDetrended[i]<EtaDetrended[i-1]) and (EtaDetrended[i]<EtaDetrended[i+1])):
                
                if EtaDetrended[i]<0:
                    
                    #Check if after last trough, a crest is detected or not 
                    #n==m-1 means it is detected and the new y(i,1) is the trough of next wave
                    if n==m-1:
                        
                        n=n+1
                        xmin=np.append(xmin,t[i])
                        ymin=np.append(ymin,EtaDetrended[i])
                        xminIndx=np.append(xminIndx,i)
                        
                    elif n!=-1:
                            
                        #Replacingthe old crest location with new one if the new one is smaller
                        if EtaDetrended[i]<ymin[n]:
                            
                            xmin[n]=t[i]
                            ymin[n]=EtaDetrended[i]
                            xminIndx[n]=i
                                
    #Eta=EtaDetrended.copy() #Water surface elevation time series

    #Converting to int
    xmaxIndx=np.int_(xmaxIndx)
    xminIndx=np.int_(xminIndx)

    #--------------------------------------------------------------------------
    #Calculating wave height and surface elevation of the wave crest and trough

    len1=len(xmax)
    len2=len(xmin)
    
    if len1!=len2:
        import sys
        sys.exit('Numbers of crests and troughs are not the same')
    

    H=np.zeros(len1)
    xmean=np.zeros(len1)
    Etac=np.zeros(len1)
    Etat=np.zeros(len1)
    
    for i in range(0,len1,1):
        
        H[i]=ymax[i]-ymin[i] #Wave height
        xmean[i]=(xmax[i]+xmin[i])/2
        Etac[i]=ymax[i] #Water surface elevation of the wave crest
        Etat[i]=ymin[i] #Water surface elevation of the wave trough
        
    #--------------------------------------------------------------------------
    #Calculating wave period 

    T=np.zeros(len(H))
    
    for i in range(0,len(H),1):
        for j in range (0,len(xupcross)-1,1):
            
            if ((xupcross[j]<xmean[i]) and (xupcross[j+1]>xmean[i])):
                T[i]=xupcross[j+1]-xupcross[j]

    #--------------------------------------------------------------------------
    #Calculating wave properties

    numberofwaves=len(H)
    Etarms=np.sqrt(np.sum(EtaDetrended**2)/N)
    
    HsortIndex=np.argsort(H)[::-1] #Sorting in descending order
    Hsort=np.sort(H)[::-1] #Sorting in descending order
    HTop3rdIndx=int(np.round(1/3*len(Hsort))) #Index where top 3rd wave height starts
    Hs=np.mean(Hsort[0:HTop3rdIndx+1]) #Zero-crossing significant wave height
    
    Tsort=T[HsortIndex]
    TTop3rdIndx=int(np.round(1/3*len(Tsort))) #Index where top 3rd wave period starts
    Ts=np.mean(Tsort[0:TTop3rdIndx+1]) #Zero-crossing significant wave period
    
    Hrms=np.sqrt(np.sum(H**2)/len(H)) #Root Mean Square Wave Height
    Hz=np.mean(H) #Zero-crossing mean wave height
    Tz=np.mean(T) #Zero-crossing mean wave period
    
    #Hs=np.sqrt(2)*np.sqrt(np.sum(H**2)/len(H)) #Significant wave height, Hs=(2^0.5)*Hrms
    #Hz=np.sum(H)/numberofwaves #Zero-crossing mean wave height
    #Tz=duration/numberofwaves #Zero-crossing mean wave period

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
    
        val=[Hs,Ts,Hz,Tz,Hrms]
        name=['Hs','Ts','Hz','Tz','Hrms']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            print('{0:10}= {1:0.10f}'.format(name[i],val[i])) 
    
        #Plotting data time series
        plt.plot(t,EtaDetrended,label='Time Series')
        
        plt.scatter(xmax,ymax,marker='*',label='Crest')
        plt.scatter(xmin,ymin,marker='o',label='Trough')
        
        plt.scatter(xupcross,yupcross,marker='^',label='Up Crossing')
        plt.scatter(xdowncross,ydowncross,marker='v',label='Down Crossing')
    
        for i in range (0,len(H),1):
            
            if i==0:
                plt.plot([xmean[i],xmean[i]],[(ymax[i]+ymin[i])/2-H[i]/2,(ymax[i]+ymin[i])/2+H[i]/2],':m',label='Wave Height')
            else:
                plt.plot([xmean[i],xmean[i]],[(ymax[i]+ymin[i])/2-H[i]/2,(ymax[i]+ymin[i])/2+H[i]/2],':m')
            
        plt.xlabel('Time (s)')
        plt.ylabel('Water Surface Elevatio (m)')
        plt.legend()
    
    #--------------------------------------------------------------------------
    #Outputs
    return Hs, Ts, Hz, Tz, Hrms, H, T

    #--------------------------------------------------------------------------
