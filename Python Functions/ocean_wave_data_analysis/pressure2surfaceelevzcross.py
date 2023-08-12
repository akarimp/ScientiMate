def pressure2surfaceelevzcross(P, fs, h, heightfrombed=0, Kpmin=0.15, kCalcMethod='beji', Rho=1000, dispout='no'):
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

    scientimate.pressure2surfaceelevzcross
    ======================================

    .. code:: python

        Eta, t = scientimate.pressure2surfaceelevzcross(P, fs, h, heightfrombed=0, Kpmin=0.15, kCalcMethod='beji', Rho=1000, dispout='no')

    Description
    -----------

    Calculate water surface elevation time series from water pressure time series by using an upward zero crossing method

    Inputs
    ------

    P
        Water pressure time series data in (N/m^2)
    fs
        Sampling frequency that data collected at in (Hz)
    h
        Water depth in (m)
    heightfrombed=0
        Height from bed that data collected at in (m)
    Kpmin=0.15
        | Minimum acceptable value for a pressure response factor
        | If Kpmin=0.15, it avoid wave amplification larger than 6 times (1/0.15)
    kCalcMethod='beji'
        | Wave number calculation method 
        | 'hunt': Hunt (1979), 'beji': Beji (2013), 'vatankhah': Vatankhah and Aghashariatmadari (2013) 
        | 'goda': Goda (2010), 'exact': calculate exact value 
    Rho=1000
        Water density (kg/m^3)
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    Eta
        Water surface elevation time series in (m)
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
        Eta=sp.signal.detrend(0.5*np.cos(2*np.pi*0.2*t)+(-0.1+(0.1-(-0.1)))*np.random.rand(N))
        hfrombed=4
        h=5
        k=0.2
        P=Eta*9.81*1000*(np.cosh(k*hfrombed)/np.cosh(k*h))
        Eta,t=sm.pressure2surfaceelevzcross(P,fs,5,4,0.15,'beji',1025,'yes')

    References
    ----------

    Beji, S. (2013). 
    Improved explicit approximation of linear dispersion relationship for gravity waves. 
    Coastal Engineering, 73, 11-12.

    Goda, Y. (2010). 
    Random seas and design of maritime structures. 
    World scientific.

    Hunt, J. N. (1979). 
    Direct solution of wave dispersion equation. 
    Journal of the Waterway Port Coastal and Ocean Division, 105(4), 457-459.

    Vatankhah, A. R., & Aghashariatmadari, Z. (2013). 
    Improved explicit approximation of linear dispersion relationship for gravity waves: A discussion. 
    Coastal engineering, 78, 21-22.

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
    from scipy import optimize
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
    
    P=type2numpy(P)
    h=type2numpy(h)

    #--------------------------------------------------------------------------

    #Deterending input data
    PDetrended=sp.signal.detrend(P)
    EtaDetrended=PDetrended/(Rho*9.81) #Converting pressure data to water elevation (depth) data, P=Rho.g.h
    #if (h<=0): h=0.001

    #--------------------------------------------------------------------------

    #Generating time vector
    N=len(P) #Length of a time series, total number of points in input file, N=fs*duration
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
    xmaxIndx=np.int64(xmaxIndx)
    xminIndx=np.int64(xminIndx)

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
    #==========================================================================
    #--------------------------------------------------------------------------
    #Calculating wave number (k)

    f=1/T #Wave frequency
    w=2*np.pi*f #Wave angular frequency, w=2*pi/T=2*pi*f
    
    #Deep water
    k0=(w**2)/9.81 #Deep water wave number
    k0h=k0*h
    
    #Estimation of wave number (k) from Hunt (1979)
    if kCalcMethod=='hunt':
        # c2gh=(k0h+(1.0+0.6522*k0h+0.4622*k0h**2+0.0864*k0h**4+0.0675*k0h**5)**(-1))**(-1) #Calculating wave number from Hunt (1979)
        # k=w/(np.sqrt(c2gh*9.81*h))
        d1=0.6666666667; d2=0.3555555556; d3=0.1608465608; d4=0.0632098765; d5=0.0217540484; d6=0.0065407983;
        kh=np.sqrt(k0h**2+(k0h)/(1+d1*k0h+d2*k0h**2+d3*k0h**3+d4*k0h**4+d5*k0h**5+d6*k0h**6)) #Calculating wave number from Hunt (1979)
        k=kh/h #Calculating wave number from Hunt (1979)
    
    
    #Estimation of wave number (k) from Beji (2013)
    elif kCalcMethod=='beji':
        kh=k0h*(1+k0h**1.09*np.exp(-(1.55+1.3*k0h+0.216*k0h**2)))/np.sqrt(np.tanh(k0h)) #Calculating wave number from Beji (2013)
        k=kh/h #Calculating wave number from Beji (2013)
        k[w==0]=0
    
    
    #Estimation of wave number (k) from Vatankhah and Aghashariatmadari (2013)
    elif kCalcMethod=='vatankhah':
        kh=(k0h+k0h**2*np.exp(-(3.2+k0h**1.65)))/np.sqrt(np.tanh(k0h))+k0h*(1-np.exp(-(k0h**0.132)))**(5.0532+2.158*k0h**1.505) #Calculating wave number from Vatankhah and Aghashariatmadari (2013)
        k=kh/h #Calculating wave number from Vatankhah and Aghashariatmadari (2013)
        k[w==0]=0


    #Estimation of wave number (k) from Goad (2010)
    elif kCalcMethod=='goda':
        kh=np.zeros(len(k0h))
        kh[k0h>=1]=k0h[k0h>=1]
        kh[k0h<1]=(k0h[k0h<1])**0.5
        for i in range(0,3,1):
            kh=kh-((kh-k0h*(np.tanh(kh))**-1)/(1+k0h*((np.tanh(kh))**(-2)-1))) #Calculating wave number from Goda (2010)

        k=kh/h #Calculating wave number from Goda (2010)
        k[w==0]=0
    
    
    #Calculating exact wave number (k)
    elif kCalcMethod=='exact':
    
        #Estimation of wave number (k) from Beji (2013)
        kh=k0h*(1+k0h**1.09*np.exp(-(1.55+1.3*k0h+0.216*k0h**2)))/np.sqrt(np.tanh(k0h)) #Calculating wave number from Beji (2013)
        kini=kh/h #Initial value for k (Wave number from Beji (2013))
        kini[w==0]=0
    
        #Calculating exact wave number (k)
        k=np.zeros(len(f)) #Pre-assigning array to make program run faster
        for i in range(0,len(f),1):
            if len(h)==1:
                fun = lambda x : w[i]**2-(9.81*x*np.tanh(x*h)) #w**2=g*k*tanh(kh)
            else:
                fun = lambda x : w[i]**2-(9.81*x*np.tanh(x*h[i])) #w**2=g*k*tanh(kh)
    
            k[i]=sp.optimize.fsolve(fun,kini[i]) #Wave number

    #--------------------------------------------------------------------------
    #Calculating the pressure response factor, Kp
    Kp=np.cosh(k*heightfrombed)/np.cosh(k*h)

    #Calculating KmaxL based on linear wave theory
    kmaxL=np.pi/(h-heightfrombed) #Wave number associated with fmaxpcorrL
    
    #Comparing Kp with KpminL calculated based on linear wave theory
    KpminL=np.cosh(kmaxL*heightfrombed)/np.cosh(kmaxL*h) #Minimum Limit for Kp calculated based on linear wave theory
    Kp[Kp<KpminL]=KpminL #Check to avoid large amplification, Kp should be larger than minimum Kp calculated based on linear wave theory

    #Check to avoid wave amplification larger than defined by Kpmin
    Kp[Kp<Kpmin]=Kpmin 

    #--------------------------------------------------------------------------
    #Pressure attenuation correction

    EtawithKp=np.zeros(N)

    for i in range(0,len(H),1):
        for j in range(xupcrossIndx[i],xupcrossIndx[i+1]+1,1):
            EtawithKp[j]=EtaDetrended[j]/Kp[i]

    #Data before very first wave, first Kp is used due lack of enough information
    for j in range(0,xupcrossIndx[0]+1,1):
        EtawithKp[j]=EtaDetrended[j]/Kp[0] 
    
    #Data after very last wave, last Kp is used due lack of enough information
    for j in range(xupcrossIndx[-1],len(EtaDetrended),1):
        EtawithKp[j]=EtaDetrended[j]/Kp[-1] 
    
    Eta=EtawithKp.copy() #Water surface elevation time series

    #--------------------------------------------------------------------------
    #==========================================================================
    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':

        #Plotting data time series
        plt.plot(t,Eta)
        #plt.plot(t,EtaDetrended,':')
    
        plt.xlabel('Time (s)')
        plt.ylabel('Water Surface Elevatio (m)')

    #--------------------------------------------------------------------------
    #Outputs
    return Eta, t

    #--------------------------------------------------------------------------
