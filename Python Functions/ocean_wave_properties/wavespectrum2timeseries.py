def wavespectrum2timeseries(f, Sxx, fs=2, dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2018-05-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.wavespectrum2timeseries
    ===================================

    .. code:: python

        Eta, t, Hm0, fp, fEta, SxxEta, a, w, Phi = scientimate.wavespectrum2timeseries(f, Sxx, fs=2, dispout='no')

    Description
    -----------

    | Generate random water wave data from a given water wave spectrum using wave superposition
    | For more options use psd2timeseries

    Inputs
    ------

    f
        Frequency (Hz)
    Sxx
        Wave power spectral density (m^2s)
    fs=2
        Sampling frequency that data are collected at in (Hz)
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    Eta
        Water surface level time series in (m)
    t
        Time in (s)
    Hm0
        Zero moment wave height (m)
    fp
        Peak wave frequency (Hz), fp=1/Tp (Tp: Peak wave period (s))
    fEta
        Frequency from generated time series(Hz)
    SxxEta
        Power spectral density from generated time series (m^2s)
    a
        Wave amplitude for for one-sided spectrum (0<fEta<fs/2) from generated time series (m)
    w
        Wave angular frequency for for one-sided spectrum (0<fEta<fs/2) from generated time series (rad/s)
    Phi
        Wave random phase for for one-sided spectrum (0<fEta<fs/2) from generated time series (rad)

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        N=2**11
        fs=8
        df=fs/N #Frequency difference 
        f=np.arange(0,fs/2+df,df) #Frequency vector 
        f[0]=f[1]/2 #Assign temporarily non-zero value to fisrt element of f to prevent division by zero
        Sxx=0.016*9.81**2/((2*np.pi)**4*(f**5))*np.exp(-1.25*(0.33/f)**4) #Calculating Spectrum 
        f[0]=0
        Sxx[0]=0
        Eta,t,Hm0,fp,fEta,SxxEta,a,w,Phi=sm.wavespectrum2timeseries(f,Sxx,fs,'yes')

    References
    ----------

    Branlard, E. (2010).
    Generation of time series from a spectrum.
    Technical University Denmark. National Laboratory for Sustainable Energy.

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
    from numpy import random
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
    
    f=type2numpy(f)
    Sxx=type2numpy(Sxx)

    #--------------------------------------------------------------------------
    #projecting f to f1(0:fs/2,1) and Sxx to Sxx1(0:fs/2,1)
    df=f[1]-f[0] #Delta f
    
    lenf=len(f)

    f1=np.empty(0)
    Sxx1=np.empty(0)

    if f[0]>df:
        f11=np.arange(f[0]-df,0,-df)
        f11=np.flipud(f11)
        lenf11=len(f11)
        f1=np.concatenate((f1,f11))
        Sxx1[0:lenf11]=0
    else:
        lenf11=0
     
    f1=np.concatenate((f1,f))
    Sxx1=np.concatenate((Sxx1,Sxx))
    
    if f[-1]<fs/2:
        f12=np.arange(f[-1]+df,fs/2+df,df)
        lenf12=len(f12)
        f1=np.concatenate((f1,f12))
        Sxx1=np.concatenate((Sxx1,np.zeros(lenf12)))
    else:
        lenf12=0
     
    fEta=f1
    sample=len(f1) #Number of sample for 0<f<fs/2 which is equal to int(N/2+1), i.e [0:N/2]
    N=2*(sample-1) #Total number of points between 0<f<fs is N+1 where int(N/2+1)=sample is a total number of points between 0<f<fs/2
    dt=1/fs #Time difference between consecutive samples, dt=1/fs
    duration=N*dt #Total time of time series
    t=np.arange(0,N*dt,dt) #Time from 0 to T-dt, equally spaced at dt

    #--------------------------------------------------------------------------
    #Calculating random time series, using wave superposition

    #----------------------------------------------------------------------
    #Calculating random phase
    mu=0 #mean=0
    sigma=1 #standard deviation=1
    rng = np.random.default_rng()
    RandNum=sigma*rng.standard_normal(int(N/2+1))+mu  # Random number with mean=0 and standard deviation=1 (normal distribution)
    Phi=2*np.pi*RandNum # Random phase
    
    #----------------------------------------------------------------------
    #Calculating wave properties
    #Note: Hm0=4*sqrt(Sxx*df), Hrms=Hm0/sqrt(2)=4/sqrt(2)*sqrt(Sxx*df), a=Hrms/2=4/(2*sqrt(2))*sqrt(Sxx*df)=sqrt(2*Sxx*df)

    Hm01=4*(Sxx1*df)**0.5 #Wave height for each deltaf from one-sided spectrum
    H=Hm01/np.sqrt(2) #Wave height: Hrms=Hm0/sqrt(2)
    w=2*np.pi*f1 #Wave angular frequency
    
    #----------------------------------------------------------------------
    #Calculating wave amplitude for one-sided spectrum
    #Note: Hm0=4*sqrt(Sxx*df), Hrms=Hm0/sqrt(2)=4/sqrt(2)*sqrt(Sxx*df), a=Hrms/2=4/(2*sqrt(2))*sqrt(Sxx*df)=sqrt(2*Sxx*df)
    #a=np.sqrt(2*Sxx1*df) #Wave amplitude for 0<f<fs/2
    a=H[0:int(N/2+1)]/2 #Wave amplitude: a=H/2 for 0<f<fs/2

    #----------------------------------------------------------------------
    #Calculating random time series using wave superposition
        
    Etaij=np.zeros((len(t),int(N/2+1))) #Pre-assigning array to make program run faster
    for i in range(0,int(N/2+1)):
        Etaij[:,i]=a[i]*np.cos(-w[i]*t+Phi[i])

    Eta=np.sum(Etaij,axis=1)

    #----------------------------------------------------------------------
    #Calculating spectrum
    SxxEta=a**2/2/df # Wave energy spectrum Sxx=1/2*a^2/deltaf
    
    #--------------------------------------------------------------------------
    #Calculating wave properties

    #Calculating spectrum using Welch method
    #[Sxxwelch,fwelch] = pwelch(Eta,hanning(256),[],N,fs); #Wave power spectrum and Frequency
    #[SxxPG,fPG] = periodogram(Eta,hamming(length(Eta(:,1))),length(Eta(:,1)),fs); #Wave power spectrum and Frequency
    #SxxPGSmooth=smooth(SxxPG,0.01,'lowess');
    
    #Calculating spectral moments
    m0=np.sum(SxxEta*fEta**0*df)
    m1=np.sum(SxxEta*fEta**1*df)
    m2=np.sum(SxxEta*fEta**2*df)
    
    #Calculating wave properties
    Hm0=4*np.sqrt(m0) #Zero-Moment wave height
    Tm01=m0/m1 #mean period
    Tm02=(m0/m2)**0.5 #zero crossing period
    
    #Calculating peak period
    SxxmaxIndx=np.argmax(Sxx) #Locating an index of the spectrum peak
    fp=f[SxxmaxIndx] #peak frequency
    Tp=1/fp #peak period

    #Calculating peak frequency from weighted integral (Young, 1995)
    #fp=(np.sum(Sxx**5*f**1*df))/(np.sum(Sxx**5*f**0*df)) #Peak frequency

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
        
        val=[Hm0,fp,Tp,Tm01,Tm02,m0,m1,m2]
        name=['Hm0','fp','Tp','Tm01','Tm02','m0','m1','m2']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            print('{0:10}= {1:0.10f}'.format(name[i],val[i])) 
        
        #plotting
        plt.subplot(2,1,1)
        plt.plot(f1[f1!=0],Sxx1[f1!=0],label='Input PSD')
        #plt.hold='True'
        plt.plot(fEta[fEta!=0],SxxEta[fEta!=0],'--',label='Output PSD')
        # plt.plot(fwelch(fwelch~=0),Sxxwelch(fwelch~=0),label='Welch Sxx')
        
        plt.title('Power Spectral Density')
        plt.xlabel('Frequency(Hz)')
        plt.ylabel('Spectral Density(m^2/Hz)')
        
        #plt.legend(loc='upper right',ncol=1,fontsize=16,scatterpoints=1)
        plt.legend(loc='upper right')
        
        plt.subplot(2,1,2)
        plt.plot(t,Eta)
        #plt.hold='True'
        
        plt.title('Water Level')
        plt.xlabel('Time(s)')
        plt.ylabel('Water Level(m)')
        
    #--------------------------------------------------------------------------
    #Outputs
    return Eta, t, Hm0, fp, fEta, SxxEta, a, w, Phi

    #--------------------------------------------------------------------------
