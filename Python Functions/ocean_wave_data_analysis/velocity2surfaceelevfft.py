def velocity2surfaceelevfft(U, fs, duration, h, heightfrombed=0, fKuvmin=None, fcL=0, fcH=None, KuvafterfKuvmin='constant', kCalcMethod='beji', dispout='no'):
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

    scientimate.velocity2surfaceelevfft
    ===================================

    .. code:: python

        Eta, t = scientimate.velocity2surfaceelevfft(U, fs, duration, h, heightfrombed=0, fKuvmin=None, fcL=0, fcH=None, KuvafterfKuvmin='constant', kCalcMethod='beji', dispout='no')

    Description
    -----------

    Calculate water surface elevation time series from wave orbital velocity time series by using Fast Fourier Transform

    Inputs
    ------

    U
        Wave horizontal orbital velocity data in (m/s)
    fs
        Sampling frequency that data collected at in (Hz)
    duration
        Duration time that data are collected (second)
    h
        Water depth in (m)
    heightfrombed=0
        Height from bed that data collected at in (m)
    fKuvmin=fs/2
        Frequency that a velocity conversion factor (Kuv) at that frequency is considered as a minimum limit for Kuv
    fcL=0
        Low cut-off frequency, between 0*fs to 0.5*fs (Hz)
    fcH=fs/2
        High cut-off frequency, between 0*fs to 0.5*fs (Hz)
    KuvafterfKuvmin='constant'
        | Define conversion factor, Kuv, value for frequency larger than fKuvmin
        | 'nochange': Kuv is not changed for frequency larger than fKuvmin 
        | 'one': Kuv=1 for frequency larger than fKuvmin 
        | 'constant': Kuv for f larger than fKuvmin stays equal to Kuv at fKuvmin (constant)
    kCalcMethod='beji'
        | Wave number calculation method 
        | 'hunt': Hunt (1979), 'beji': Beji (2013), 'vatankhah': Vatankhah and Aghashariatmadari (2013) 
        | 'goda': Goda (2010), 'exact': calculate exact value 
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
        U=(np.pi/5)*(2*Eta)*(np.cosh(k*hfrombed)/np.sinh(k*h)) 
        Eta1,t=sm.velocity2surfaceelevfft(U,fs,duration,5,4,0.6,0,fs/2,'constant','beji','yes')

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

    Wiberg, P. L., & Sherwood, C. R. (2008). 
    Calculating wave-generated bottom orbital velocities from surface-wave parameters. 
    Computers & Geosciences, 34(10), 1243-1262.

    Welch, P. (1967). 
    The use of fast Fourier transform for the estimation of power spectra: a method based on time averaging over short, modified periodograms. 
    IEEE Transactions on audio and electroacoustics, 15(2), 70-73.

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
    
    U=type2numpy(U)
    h=type2numpy(h)

    #--------------------------------------------------------------------------
    #Assign default values

    if fKuvmin is None: fKuvmin=fs/2
    if fcH is None: fcH=fs/2

    #--------------------------------------------------------------------------
    #Preventing error in scipy signal by making sure it is a single column array
    #U=U[:,0]

    #--------------------------------------------------------------------------
    #Checking number of points in input file to be even

    #if (len(U)%2)==1:
    #    U=U[0:-1] #Removing the last element of the array to make array size even

    #--------------------------------------------------------------------------

    #Deterending input data
    UDetrended=sp.signal.detrend(U)
    #if (h<=0): h=0.001
    if (fKuvmin>fs/2): fKuvmin=int(fs/2)
    if (fcL<0): fcL=0
    if (fcH>fs/2): fcH=int(fs/2)
    #nfft=int(2**(np.ceil(np.log2(len(UxDetrended)))))

    #--------------------------------------------------------------------------

    #Generating frequency vector
    N=len(U) #Length of a time series, total number of points in input file, N=fs*duration
    df=fs/N #Frequency difference between consecutive samples, df=fs/N
    f=np.arange(0,fs/2+df,df) #Frequency vector with (N/2+1) elements from 0 Hz to fNy=fs/2 Hz, equally spaced at df

    #Generating time vector
    dt=1/fs #Time difference between consecutive samples, dt=1/fs=duration/N
    # t(:,1)=linspace(dt,N/fs,N); #Time from 0 to T-dt, equally spaced at dt
    t=np.arange(0,duration,dt) #Time from 0 to T-dt, equally spaced at dt

    #--------------------------------------------------------------------------
    #Calculating the Fast Fourier Transform of the orbital velocity 

    UFFT=np.fft.fft(UDetrended,n=N) #FFT of the orbital velocity

    #--------------------------------------------------------------------------
    #Calculating wave number (k)

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
    #Calculating the velocity response factor, Kuv, (Wiberg, 2008)
    Kuv=(2*np.pi*f)*np.cosh(k*heightfrombed)/np.sinh(k*h)
    Kuv[0]=Kuv[1]
    
    #Defining a value of Kuv for frequency larger than fKuvmin
    if KuvafterfKuvmin=='nochange':
    
        Kuv
    
    elif KuvafterfKuvmin=='one':
    
        Kuv[f>fKuvmin]=1 #Kuv for f larger than fKuvmin should be 1 (no correction)
    
        #Linear change of Kuv to 1 for f larger than fKuvmin with bandwith of 0.1 Hz
        Indx1=int(np.max((np.nonzero(f<=fKuvmin-0.05))[0]))
        Indx2=int(np.max((np.nonzero(f<=fKuvmin+0.05))[0]))
        if (Indx2>len(f)): Indx2=len(f)
        for i in range(Indx1,Indx2+1,1):
            Kuv[i]=(Kuv[Indx2]-Kuv[Indx1])/(Indx2-Indx1)*(i-Indx1)+Kuv[Indx1]

    elif KuvafterfKuvmin=='constant': #Apply minimum velocity conversion factor after fKuvmin
        #Kuv[Kuv<1]=1 #Conversion factor, Kuv, should be larger than 1
        fKuvminIndx=int(np.max((np.nonzero(f<=fKuvmin))[0]))
        if (fKuvminIndx>len(f)): fKuvminIndx=len(f)
        Kuv[f>fKuvmin]=Kuv[fKuvminIndx] #Kuv for f larger than fKuvmin stays equal to Kuv at fKuvmin (constant)

    #Calculating KmaxL based on linear wave theory
    kmaxL=np.pi/(h-heightfrombed) #Wave number associated with fmaxuvcorrL
    
    #Calculating fmaxuvcorrL based on linear wave theory
    fmaxuvcorrL=1/(2*np.pi)*np.sqrt(9.81*np.pi/(h-heightfrombed)*np.tanh(np.pi*h/(h-heightfrombed))) #fmaxuvcorrL based on linear wave theory
    
    #Comparing Kuv with KuvminL calculated based on linear wave theory
    KuvminL=(2*np.pi*fmaxuvcorrL)*np.cosh(kmaxL*heightfrombed)/np.sinh(kmaxL*h) #Minimum Limit for Kuv calculated based on linear wave theory
    Kuv[Kuv<KuvminL]=KuvminL #Check to avoid large amplification, Kuv should be larger than minimum Kuv calculated based on linear wave theory

    #Make Kuv length equal to N and symetric around fs/2
    Kuv=np.concatenate((Kuv,np.flipud(Kuv[1:-1]))) #Make Kuv symetric around fs/2

    #--------------------------------------------------------------------------
    #Replacing FFT values by zero based on fcL and fcH

    #import matplotlib.pylab as plt
    #plt.plot(np.arange(0,N,1)*df,np.abs(UFFT)) #Plot before applying fcL and fcH

    if ((fcL>f[0]) and (fcL<fs/2)):
        UFFT_1sthalf_1=UFFT[0:len(f)].copy()
        UFFT_1sthalf_1[f<fcL]=0
        UFFT=np.concatenate((UFFT_1sthalf_1,np.flipud(UFFT_1sthalf_1[1:-1]))) #Make UFFT symetric around fs/2
        #UFFT[f<fcL]=0

    if ((fcH>f[0]) and (fcH<fs/2)):
        UFFT_1sthalf_2=UFFT[0:len(f)].copy()
        UFFT_1sthalf_2[f>fcH]=0
        UFFT=np.concatenate((UFFT_1sthalf_2,np.flipud(UFFT_1sthalf_2[1:-1]))) #Make UFFT symetric around fs/2
        #UFFT[f>fcH]=0

    #plt.plot(np.arange(0,N,1)*df,np.abs(UFFT),'--') #Plot after applying fcL and fcH

    #--------------------------------------------------------------------------
    #Calculating FFT of the water surface elevation

    EtaFFT= UFFT/Kuv #Calculating water elevation FFT from orbital velocity FFT
    Eta=np.real(np.fft.ifft(EtaFFT,n=N))	#Water surface elevation time series

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':

        plt.plot(t,Eta)
        plt.xlabel('Time (s)')
        plt.ylabel('Water Surface Elevation (m)')

    #--------------------------------------------------------------------------
    #Outputs
    return Eta, t

    #--------------------------------------------------------------------------
