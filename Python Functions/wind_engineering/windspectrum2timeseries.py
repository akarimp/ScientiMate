def windspectrum2timeseries(f, Sxx, fs=None, dispout='no'):
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
    
    windspectrum2timeseries
    =======================
    
    .. code:: python
    
        [Eta, t] = windspectrum2timeseries(f, Sxx, fs, dispout)
    
    Description
    -----------
    
    | Generate zero-mean wind velocity time series from a given spectrum
    | Add mean of the wind velocity to generated time series to get desired time series
    
    Inputs
    ------
    
    f
        Frequency (Hz)
    Sxx
        | Power spectral density in ((m/s)^2/Hz)
        | Length of Sxx and f should be odd number
    fs=2*max(f)
        Sampling frequency that data collected at in (Hz)
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)
    
    Outputs
    -------
    
    U
        Wind velocity time series in (m/s)
    t
        Time in (s)
    
    Examples
    --------
    
    .. code:: python
    
        import scientimate as sm
        import numpy as np

        N=2048+1 #Total number of points
        fs=2 #Sampling frequency
        df=fs/N #Frequency difference 
        f=np.arange(0.002,fs/2+df,df) #Frequency vector 
        Sxx=1.0**2*(6.868*100/10)/(1+10.32*f*100/10)**(5/3) #Calculate Spectrum
        [U, t]=sm.windspectrum2timeseries(f, Sxx, fs, 'yes')
    
    References
    ----------
    
    Branlard, E. (2010).
    Generation of time series from a spectrum.
    Technical University Denmark. National Laboratory for Sustainable Energy.
    
    Rose, S., & Apt, J. (2012). 
    Generating wind time series as a hybrid of measured and simulated data. 
    Wind Energy, 15(5), 699-715.
    
    Shinozuka, M., & Jan, C. M. (1972). 
    Digital simulation of random processes and its applications. 
    Journal of sound and vibration, 25(1), 111-128.
    
    Veers, P. (1984). 
    Modeling stochastic wind loads on vertical axis wind turbines. 
    In 25th Structures, Structural Dynamics and Materials Conference (p. 910).
    
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
    
    f=type2numpy(f)
    Sxx=type2numpy(Sxx)

    #--------------------------------------------------------------------------
    #Initial values

    if fs is None: fs=2*np.max(f)

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
    
    f_U=f1
    #sample=len(f1) #Number of sample for 0<f<fs/2 which is equal to int(N/2+1), i.e [0:N/2]
    sample=len(f) #Number of sample for 0<f<fs/2 which is equal to int(N/2+1), i.e [0:N/2]
    N=2*(sample-1) #Total number of points between 0<f<fs is N+1 where int(N/2+1)=sample is a total number of points between 0<f<fs/2
    dt=1/fs #Time difference between consecutive samples, dt=1/fs
    duration=N*dt #Total time of time series
    t=np.arange(0,N*dt,dt) #Time from 0 to T-dt, equally spaced at dt

    #--------------------------------------------------------------------------
    #Calculating random time series, using wave superposition
    
    #----------------------------------------------------------------------
    #Calculate random phase
    mu=0 #mean=0
    sigma=1 #standard deviation=1
    rng = np.random.default_rng()
    #RandNum=sigma*rng.standard_normal(int(N/2+1))+mu  # Random number with mean=0 and standard deviation=1 (normal distribution)
    RandNum=rng.random(int(N/2+1))  #Random number (uniform distribution)
    Phi=2*np.pi*RandNum # Random phase
    
    #----------------------------------------------------------------------
    #Calculate wave properties
    
    dw=2*np.pi*df
    an=np.sqrt(1/2*Sxx*df)*np.sin(Phi) #Veers (1984)
    bn=np.sqrt(1/2*Sxx*df)*np.cos(Phi) #Veers (1984)
    cn=np.sqrt(Sxx*df) #Shinozuka & Jan (1972)
    #cn[0]=0
    w=2*np.pi*f #Wave angular frequency
    
    #----------------------------------------------------------------------
    #Calculating random time series using wave superposition
    
    #Shinozuka, M., & Jan, C. M. (1972). 
    #Digital simulation of random processes and its applications. 
    #Journal of sound and vibration, 25(1), 111-128.
    #Eq(15)
    
    Uij=np.zeros((len(t),int(N/2+1))) #Pre-assigning array to make program run faster
    for i in range(0,int(N/2+1),1):
        #Uij[:,i]=an[i]*np.sin(w[i]*t)+bn[i]*np.cos(w[i]*t) #Veers (1984)
        Uij[:,i]=np.sqrt(2)*cn[i]*np.cos(w[i]*t+Phi[i]) #Shinozuka & Jan (1972)

    U=np.sum(Uij,axis=1)
    #U=2*U #Veers (1984)
    
    #----------------------------------------------------------------------
    #Calculate spectrum
    
    #fwelch,Sxxwelch=sp.signal.welch(U, fs=fs, nfft=N) #Wave power spectrum and Frequency
    #print(np.trapz(Sxx,f)) #Input spectrum
    #print(np.trapz(Sxxwelch,fwelch)) #Spectrum from generated time series
    
    #plt.plot(np.sqrt(Sxxwelch/Sxx))
    #plt.figure(2)
    
    #--------------------------------------------------------------------------
    #Displaying results
    
    if dispout=='yes':
        
        plt.subplot(2,1,1)
        plt.plot(f[f!=0],Sxx[f!=0], label='Input PSD')
        #plt.plot(fwelch[fwelch!=0],Sxxwelch[fwelch!=0], label='Output PSD')
        
        plt.title('Power Spectral Density')
        plt.xlabel('Frequency(Hz)')
        plt.ylabel('Spectral Density ((m/s)^2/Hz)')
        
        #plt.legend()
        
        plt.subplot(2,1,2)
        plt.plot(t,U)
        
        #plt.title('Wind Velocity')
        plt.xlabel('Time(s)')
        plt.ylabel('Wind Velocity (m/s)')
        
    
    #--------------------------------------------------------------------------
    #Outputs
    return U, t

    #--------------------------------------------------------------------------
    