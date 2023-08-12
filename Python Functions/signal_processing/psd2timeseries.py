def psd2timeseries(f, Sxx, fs=2, CalcMethod='fft', dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2017-01-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.psd2timeseries
    ==========================

    .. code:: python

        Eta, t, Hm0, fp, fEta, SxxEta, a = scientimate.psd2timeseries(f, Sxx, fs=2, CalcMethod='fft', dispout='no')

    Description
    -----------

    Generate random wave data from a given spectrum

    Inputs
    ------

    f
        Frequency (Hz)
    Sxx
        Power spectral density (m^2s)
    fs=8
        Sampling frequency that data collected at in (Hz)
    CalcMethod='fft'
        | Method for Calculating random time series, 
        | 'fft': using Fast Fourier Transform, 'sp': using wave superposition
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    Eta
        Water Surface Level Time Series in (m)
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
        Eta,t,Hm0,fp,fEta,SxxEta,a=sm.psd2timeseries(f,Sxx,fs,'fft','yes')

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
    #%%
    #projecting f to f1(0:fs/2,1) and Sxx to Sxx1(0:fs/2,1)
    df=f[1]-f[0] #Delta f
    
    lenf=len(f)

    f1=np.empty(0)
    Sxx1=np.empty(0)

    if f[0]>df:
        f11=np.arange(f[0]-df,0,-df)
        f11=np.flipud(f11)
        lenf11=len(f11)
        #f1=np.zeros(lenf11)f1[0:lenf11]=f11
        f1=np.concatenate((f1,f11))
        Sxx1[0:lenf11]=0
    else:
        lenf11=0
     
    #f1[lenf11:lenf11+lenf]=f
    #Sxx1[lenf11+1:lenf11+lenf]=Sxx
    f1=np.concatenate((f1,f))
    Sxx1=np.concatenate((Sxx1,Sxx))
    
    if f[-1]<fs/2:
        f12=np.arange(f[-1]+df,fs/2+df,df)
        lenf12=len(f12)
        #f1[lenf11+lenf+1:lenf11+lenf+lenf12]=f12
        f1=np.concatenate((f1,f12))
        #Sxx1[lenf11+lenf+1:lenf11+lenf+lenf12]=0
        Sxx1=np.concatenate((Sxx1,np.zeros(lenf12)))
    else:
        lenf12=0
     
    fEta=f1
    sample=len(f1) #Number of sample for 0<f<fs/2 which is equal to int(N/2+1), i.e [0:N/2]
    N=2*(sample-1) #Total number of points between 0<f<fs is N+1 where int(N/2+1)=sample is a total number of points between 0<f<fs/2
    dt=1/fs #Time difference between consecutive samples, dt=1/fs
    #duration=N/fs #Total time of time series
    duration=N*dt #Total time of time series
    #t=np.linspace(dt,N/fs,N) #Time from 0 to T-dt, equally spaced at dt
    #t=np.arange(0,duration,dt) #Time from 0 to T-dt, equally spaced at dt
    t=np.arange(0,N*dt,dt) #Time from 0 to T-dt, equally spaced at dt
    
    #--------------------------------------------------------------------------
    # converts one-sided spectrum into two-sided
    Sxx2=Sxx1.copy() 
    Sxx2[1:-1]=Sxx1[1:-1]/2 #division by 2 converts one-sided spectrum into two-sided, 
                             #by reducing the energy in two-sides to half to make it equal with one-sided spectrum
    #Sxx3=np.zeros(N+1)
    #Sxx3[0:int(N/2+1)]=Sxx2 #make Sxx symetric around fs/2
    #Sxx3[int(N/2+1):N+2]=np.flipud(Sxx2[0:-1]) #make Sxx symetric around fs/2
    Sxx3=np.concatenate((Sxx2,np.flipud(Sxx2[1:-1]))) #make Sxx symetric around fs/2

    #--------------------------------------------------------------------------
    #%%
    # Calculating random time series, First method, using FFT

    if CalcMethod=='fft':
        
        #----------------------------------------------------------------------
        #Calculating FFT from Spectrum
        
        # Note: Hm0=4*sqrt(Sxx*df), Hrms=Hm0/sqrt(2)=4/sqrt(2)*sqrt(Sxx*df), a=Hrms/2=4/(2*sqrt(2))*sqrt(Sxx*df)=sqrt(2*Sxx*df)
        # Note: total number of points=N
        
        # method 1, 1st approach
        EtaFFT1=np.sqrt(2*df*Sxx3[0:N+1]) #Wave amplitude from two-sided spectrum, Eta=Sigma(dn*cos(nwt+phi)), Sxx.df=1/(2)*(dn)^2, dn=sqrt(2.Sxx.df)
        EtaFFT2=EtaFFT1*(N) # in Matlab Y=fft(y,N)/length(y)
        
        # method 1, 2nd approach
        # EtaFFT2=sqrt(2*(N+1)*fs*Sxx3(1:N+1,1)); #Wave amplitude from two-sided spectrum, Sxx=2*(1/(fs*N))*abs(fft)^2
        
        #----------------------------------------------------------------------
        #Calculating wave amplitude for one-sided spectrum
        dn=np.sqrt(2*df*Sxx1) #Wave amplitude from one-sided spectrum, Eta=Sigma(dn*cos(nwt+phi)), Sxx.df=1/(2)*(dn)^2, dn=sqrt(2.Sxx.df)
        a=dn #Wave amplitude: a=dn=H/2 for 0<f<fs/2
        
        #----------------------------------------------------------------------
        #Calculating white noise
        mu=0 #mean=0
        sigma=1 #standard deviation=1
        RandNum=sigma*np.random.randn(int(N/2+1))+mu  #Random number with mean=0 and standard deviation=1 (normal distribution)
        #RandNum=(-1)+(1-(-1))*np.random.rand(int(N/2+1)) #Random number between -1 and 1 (uniform distribution)
        Phi=2*np.pi*RandNum # Random phase
        WhiteNoise1=sigma*np.exp(1j*Phi) # Generating white noise
        
        #WhiteNoise2=np.flipud(WhiteNoise1[0:-1])
        #WhiteNoise=WhiteNoise1.copy()
        #WhiteNoise[int(N/2+1):N+1]=WhiteNoise2 #make White Noise symetric around fs/2

        WhiteNoise=np.concatenate((WhiteNoise1,np.flipud(WhiteNoise1[1:-1]))) #make White Noise symetric around fs/2
                
        #----------------------------------------------------------------------
        #Calculating random time series
        
        #Adding white noise to spectrum (double-sided) in frequency-domain
        EtaFFT=EtaFFT2*WhiteNoise # (meters)
        
        #Calculating time series, method 1
        Eta=np.real(np.fft.ifft(EtaFFT[0:N+1],N))	# corected water surface levels time series, total number of points=N+1
        Eta=Eta[0:N] #change the size Eta equal to size of t
        
        #Calculating random time series, method 1, 3rd approach 
        # w(:,1)=2*pi*f1; # Angular frequency
        
        # for i=1:N/2
        #     Etaij(:,i)=(EtaFFT1(i,1)/sqrt(2))*cos(w(i,1)*t+Phi(i,1)); #Hrms=Hm0/sqrt(2)
        # end
        
        # for i=int(N/2+1):N
        #     Etaij(:,i)=Etaij(:,N/2-(i-(int(N/2+1)))); #make Wij symetric around fs/2
        # end
        
        # for i=1:N
        #     Eta1(i,1)=sum(Etaij(i,1:end));
        # end
        
        #----------------------------------------------------------------------
        #Calculating spectrum
        
        # First method: Sxx=2*(1/(dt^2/duration))*abs(fft)**2 or Sxx=2*(1/(fs*N))*abs(fft)^2
        Y=np.fft.fft(Eta,N) #calculating Fast Fourier transform
        YMagnitude=np.abs(Y) #magnitude of complex fft Y=sqrt(an^2+bn^2)
        # psd=(dt^2/duration)*(YMagnitude)**2 #calculating two-side power density spectrum, #Sxx=T/dt^2*(cn)^2
        psd=(1/(N*fs))*(YMagnitude)**2 #Calculating psd using fs
        SxxEta=np.zeros(int(N/2+1)) #number of sample for 0<f<fs/2 is equal to int(N/2+1), i.e [0:N/2]
        SxxEta[0]=psd[0] #assigning one-sided power density [0]
        SxxEta[int(N/2)]=psd[int(N/2)] #assigning one-sided power density [N/2]
        SxxEta[1:int(N/2)]=2*psd[1:int(N/2)] #calculating one-sided power density [1:N/2-1]
                                               #(multiplying by 2 converts two sides spectrum into one side, 
                                               #by double the energy in one-sided spectrum to compansete for energy from eleminated side of two-sided spectrum)
        
        #Creating Hamming window function
        Nw=int((len(SxxEta))/32)+1 #Number of elements in a window function
        if (Nw<3): Nw=3 #Window size should be larger than 3
        nwin=np.arange(0,Nw,1) #Counter from 0 to Nw-1
        WnHamm=0.54-0.46*np.cos(2*np.pi*nwin/(Nw-1)) #Hamming window function

        WnNorm=WnHamm/np.sum(WnHamm) #Normalizing a filter
        
        # smoothing data
        SxxEta=np.convolve(SxxEta,WnNorm,'same') #Smoothing a power density spectra
        Hm0_Smooth=4*np.sqrt(np.sum(SxxEta*fEta**0*df))
        
 
    #--------------------------------------------------------------------------
    #%%
    # Calculating random time series, Second method, using wave superposition

    if CalcMethod=='sp':
        #----------------------------------------------------------------------
        #Calculating random phase
        mu=0 #mean=0
        sigma=1 #standard deviation=1
        RandNum=sigma*np.random.randn(int(N/2+1))+mu  # Random number with mean=0 and standard deviation=1 (normal distribution)
        #RandNum=(-1)+(1-(-1))*rand(N/2,1); #Random number between -1 and 1 (uniform distribution)
        Phi=2*np.pi*RandNum # Random phase
        
        # Phi1=Phi;
        # Phi1=flipud(Phi1);
        # Phi(N/2+1:2*N/2,1)=Phi1; #make Random phase symetric around fs/2
        
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

        #Etaij=np.empty([len(t),0])
        #for i in range(0,int(N/2+1)):
        #    Etaij=np.column_stack((Etaij,(a[i])*np.cos(-w[i]*t+Phi[i])))
         
        #print(np.shape(Etaij))
        #Eta=np.empty(len(t))
        #for i in range(0,len(t)):
        #    Eta[i]=np.sum(Etaij[i,0:])

        #----------------------------------------------------------------------
        #Calculating spectrum
        SxxEta=a**2/2/df # Wave energy spectrum Sxx=1/2*a^2/deltaf
        
    #--------------------------------------------------------------------------
    #%%
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
        plt.plot(fEta[fEta!=0],SxxEta[fEta!=0],label='Output PSD')
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
    return Eta, t, Hm0, fp, fEta, SxxEta, a

    #--------------------------------------------------------------------------
