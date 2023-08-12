def pressure2surfaceelevfft(P, fs, duration, h, heightfrombed=0, fmaxpcorr=None, fminpcorr=0, fcL=0, fcH=None, fmaxpcorrCalcMethod='auto', Kpafterfmaxpcorr='constant', kCalcMethod='beji', Rho=1000, dispout='no'):
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

    scientimate.pressure2surfaceelevfft
    ===================================

    .. code:: python

        Eta, t = scientimate.pressure2surfaceelevfft(P, fs, duration, h, heightfrombed=0, fmaxpcorr=None, fminpcorr=0, fcL=0, fcH=None, fmaxpcorrCalcMethod='auto', Kpafterfmaxpcorr='constant', kCalcMethod='beji', Rho=1000, dispout='no')

    Description
    -----------

    Calculate water surface elevation time series from water pressure time series by using Fast Fourier Transform

    Inputs
    ------

    P
        Water pressure time series data in (N/m^2)
    fs
        Sampling frequency that data collected at in (Hz)
    duration
        Duration time that data are collected (second)
    h
        Water depth in (m)
    heightfrombed=0
        Height from bed that data collected at in (m)
    fmaxpcorr=fs/2
        | Maximum frequency that a pressure attenuation factor applies up on that (Hz)
        | If fmaxpcorrCalcMethod='user', then the smaller of calculated and user defined fmaxpcorr will be chosen
    fminpcorr=0
        | Minimum frequency that is used for defining fmaxpcorr if fmaxpcorrCalcMethod='auto' (Hz)
        | fminpcorr should be smaller than fp 
        | If swell energy exists, fminpcorr should be smaller than fp of wind sea (fpsea) and larger than fp of swell (fpswell) if there swell 
    fcL=0
        Low cut-off frequency, between 0*fs to 0.5*fs (Hz)
    fcH=fs/2
        High cut-off frequency, between 0*fs to 0.5*fs (Hz)
    fmaxpcorrCalcMethod='auto'
        | Define if to calculate fmaxpcorr and ftail or to use user defined
        | 'user': use user defined value for fmaxpcorr
        | 'auto': automatically define value for fmaxpcorr
    Kpafterfmaxpcorr='constant'
        | Define a pressure response factor, Kp, value for frequency larger than fmaxpcorr
        | 'nochange': Kp is not changed for frequency larger than fmaxpcorr 
        | 'one': Kp=1 for frequency larger than fmaxpcorr 
        | 'constant': Kp for f larger than fmaxpcorr stays equal to Kp at fmaxpcorr (constant)
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
        Eta,t=sm.pressure2surfaceelevfft(P,fs,duration,5,4,0.7,0.15,0,fs/2,'auto','constant','beji',1025,'yes')

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
    
    P=type2numpy(P)
    h=type2numpy(h)

    #--------------------------------------------------------------------------
    #Assign default values

    if fmaxpcorr is None: fmaxpcorr=fs/2
    if fcH is None: fcH=fs/2

    #--------------------------------------------------------------------------
    #Preventing error in scipy signal by making sure it is a single column array
    #P=P[:,0]

    #--------------------------------------------------------------------------
    #Checking number of points in input file to be even

    #if (len(P)%2)==1:
    #    P=P[0:-1] #Removing the last element of the array to make array size even

    #--------------------------------------------------------------------------

    #Deterending input data
    PDetrended=sp.signal.detrend(P)
    EtaDetrended=PDetrended/(Rho*9.81) #Converting pressure data to water elevation (depth) data, P=Rho.g.h
    #if (h<=0): h=0.001
    if (fmaxpcorr>fs/2): fmaxpcorr=int(fs/2)
    if (fcL<0): fcL=0
    if (fcH>fs/2): fcH=int(fs/2)
    #nfft=int(2**(np.ceil(np.log2(len(EtaDetrended)))))

    #--------------------------------------------------------------------------

    #Generating frequency vector
    N=len(P) #Length of a time series, total number of points in input file, N=fs*duration
    df=fs/N #Frequency difference between consecutive samples, df=fs/N
    f=np.arange(0,fs/2+df,df) #Frequency vector with (N/2+1) elements from 0 Hz to fNy=fs/2 Hz, equally spaced at df

    #Generating time vector
    dt=1/fs #Time difference between consecutive samples, dt=1/fs=duration/N
    # t(:,1)=linspace(dt,N/fs,N); #Time from 0 to T-dt, equally spaced at dt
    t=np.arange(0,duration,dt) #Time from 0 to T-dt, equally spaced at dt

    #--------------------------------------------------------------------------
    #Calculating the Fast Fourier Transform of the water surface elevation before applying Kp 

    EtaFFTBeforeKp=np.fft.fft(EtaDetrended,n=N) #FFT of the water surface elevation before applying Kp

    #--------------------------------------------------------------------------
    #Calculating water surface elevation power spectral density and frequency from Welch method

    f1,Syy=sp.signal.welch(EtaDetrended,fs,nfft=N) 

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
    #Calculating the pressure response factor, Kp
    Kp=np.cosh(k*heightfrombed)/np.cosh(k*h)

    #Calculating KmaxL based on linear wave theory
    kmaxL=np.pi/(h-heightfrombed) #Wave number associated with fmaxpcorrL
    
    #Automatically estimating fmaxpcorr and ftailcorrection
    if fmaxpcorrCalcMethod=='auto':
    
        fminpcorrIndx=int(np.max((np.where(f<=fminpcorr))[0])) #Locating the index of fminpcorr (fmaxpcorr should be larger than fminpcorr)
        SyymaxIndx=np.argmax(Syy[fminpcorrIndx:]) #Locating the peak frequency, fp, of original spectrum before applying Kp
    
        fmaxpcorrL=1/(2*np.pi)*np.sqrt(9.81*kmaxL*np.tanh(kmaxL*h)) #Maximum frequency that Kp can be applied, calculated from linear wave theory
        fmaxpcorrLIndx=int(np.max((np.where(f<=fmaxpcorrL))[0])) #Locating the index of fmaxpcorrL
        if (fmaxpcorrLIndx<fminpcorrIndx+(SyymaxIndx)): fmaxpcorrLIndx=fminpcorrIndx+(SyymaxIndx) #Check if fmaxpcorrLIndx locataed after fp
    
        Syy1=Syy/(Kp**2) #Applying Kp on spectrum causing an increase to infinity at the tail of spectrum
        SyyminIndx=np.argmin(Syy1[fminpcorrIndx+(SyymaxIndx):fmaxpcorrLIndx+1]) #Locating the index of minimum value for Syy between fp and fmaxpcorrL
    
        fmaxpcorr1=f[fminpcorrIndx+(SyymaxIndx)+(SyyminIndx)] #Asigning the frequency of the minimum value of Syy between fp and fmaxpcorrL
        if (fmaxpcorr1>fmaxpcorrL): fmaxpcorr1=fmaxpcorrL #Check fmaxpcorr1 be smaller than fmaxpcorrL
        if ((fmaxpcorr1==f[fminpcorrIndx+(SyymaxIndx)]) and (fmaxpcorrL>f[fminpcorrIndx+(SyymaxIndx)])): fmaxpcorr1=fmaxpcorrL #If fmaxpcorrL>fp then fmaxpcorr1 should not be equal to fp
        if (fmaxpcorr>fmaxpcorr1): fmaxpcorr=fmaxpcorr1
    
        #ftail=f[fminpcorrIndx+(SyymaxIndx)+(SyyminIndx)] #Asigning the frequency of the minimum value of Syy between fp and fmaxpcorrL for ftail
        #if (ftail>fmaxpcorrL): ftail=fmaxpcorrL
        #if (ftail>fmaxpcorr1): ftail=fmaxpcorr1
        ftail=fmaxpcorr
    
    elif fmaxpcorrCalcMethod=='user':
    
        ftail=fmaxpcorr
    
    
    #Defining a value of Kp for frequency larger than fmaxpcorr
    if Kpafterfmaxpcorr=='nochange':

        Kp

    elif Kpafterfmaxpcorr=='one':
    
        Kp[f>fmaxpcorr]=1 #Kp for f larger than fmaxpcorr should be 1 (no correction)
    
        #Linear increase of Kp to 1 for f larger than fmaxpcorr with bandwith of 0.1 Hz
        Indx1=int(np.max((np.where(f<=fmaxpcorr-0.05))[0]))
        Indx2=int(np.max((np.where(f<=fmaxpcorr+0.05))[0]))
        if (Indx2>len(f)): Indx2=len(f)
        for i in range(Indx1,Indx2+1,1):
            Kp[i]=(Kp[Indx2]-Kp[Indx1])/(Indx2-Indx1)*(i-Indx1)+Kp[Indx1]

    
    elif Kpafterfmaxpcorr=='constant':
        fmaxpcorrIndx=int(np.max((np.where(f<=fmaxpcorr))[0]))
        if (fmaxpcorrIndx>len(f)): fmaxpcorrIndx=len(f)
        Kp[f>fmaxpcorr]=Kp[fmaxpcorrIndx] #Kp for f larger than fmaxpcorr stays equal to Kp at fmaxpcorr (constant)

    #Comparing Kp with KpminL calculated based on linear wave theory
    KpminL=np.cosh(kmaxL*heightfrombed)/np.cosh(kmaxL*h) #Minimum Limit for Kp calculated based on linear wave theory
    Kp[Kp<KpminL]=KpminL #Check to avoid large amplification, Kp should be larger than minimum Kp calculated based on linear wave theory

    #Make Kp length equal to N and symetric around fs/2
    Kp=np.concatenate((Kp,np.flipud(Kp[1:-1]))) #Make Kuv symetric around fs/2

    #--------------------------------------------------------------------------
    #Replacing FFT values by zero based on fcL and fcH

    #import matplotlib.pylab as plt
    #plt.plot(np.arange(0,N,1)*df,np.abs(EtaFFTBeforeKp)) #Plot before applying fcL and fcH

    if ((fcL>f[0]) and (fcL<fs/2)):
        EtaFFTBeforeKp_1sthalf_1=EtaFFTBeforeKp[0:len(f)].copy()
        EtaFFTBeforeKp_1sthalf_1[f<fcL]=0
        EtaFFTBeforeKp=np.concatenate((EtaFFTBeforeKp_1sthalf_1,np.flipud(EtaFFTBeforeKp_1sthalf_1[1:-1]))) #Make EtaFFTBeforeKp symetric around fs/2
        #EtaFFTBeforeKp[f<fcL]=0

    if ((fcH>f[0]) and (fcH<fs/2)):
        EtaFFTBeforeKp_1sthalf_2=EtaFFTBeforeKp[0:len(f)].copy()
        EtaFFTBeforeKp_1sthalf_2[f>fcH]=0
        EtaFFTBeforeKp=np.concatenate((EtaFFTBeforeKp_1sthalf_2,np.flipud(EtaFFTBeforeKp_1sthalf_2[1:-1]))) #Make EtaFFTBeforeKp symetric around fs/2
        #EtaFFTBeforeKp[f>fcH]=0

    #plt.plot(np.arange(0,N,1)*df,np.abs(EtaFFTBeforeKp),'--') #Plot after applying fcL and fcH

    #--------------------------------------------------------------------------
    #Calculating FFT of the water surface elevation

    EtaFFT=EtaFFTBeforeKp/Kp #Calculating water surface elevation FFT from surface elevation FFT before Kp
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
