def wavedirpuv(P, Ux, Uy, fs, h, Pheightfrombed=0, UVheightfrombed=0, dirCalcMethod='puv1', coordinatesys='xyz', fmaxpcorr=None, fminpcorr=0, fKuvmin=None, fcL=0, fcH=None, fmaxpcorrCalcMethod='auto', Kpafterfmaxpcorr='constant', KuvafterfKuvmin='constant', kCalcMethod='beji', Rho=1000, nfft=None, SegmentSize=256, OverlapSize=128, dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2017-05-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.wavedirpuv
    ======================

    .. code:: python

        Wavedir, theta1, theta2, f = scientimate.wavedirpuv(P, Ux, Uy, fs, h, Pheightfrombed=0, UVheightfrombed=0, dirCalcMethod='puv1', coordinatesys='xyz', fmaxpcorr=None, fminpcorr=0, fKuvmin=None, fcL=0, fcH=None, fmaxpcorrCalcMethod='auto', Kpafterfmaxpcorr='constant', KuvafterfKuvmin='constant', kCalcMethod='beji', Rho=1000, nfft=None, SegmentSize=256, OverlapSize=128, dispout='no')

    Description
    -----------

    Calculate wave direction using pressure and horizontal orbital velocity

    Inputs
    ------

    P
        Water pressure time series data in (N/m^2)
    Ux
        Wave horizontal orbital velocity data in x direction in (m/s)
    Uy
        Wave horizontal orbital velocity data in y direction in (m/s)
    fs
        Sampling frequency that data collected at in (Hz)
    h
        Water depth in (m)
    Pheightfrombed=0
        Pressure sensor height from bed that data collected at in (m)
    UVheightfrombed=0
        Velocity sensor height from bed that data collected at in (m)
    dirCalcMethod='puv1'
        | Wave number calculation method 
        | 'puv1': PUV Method 1, 'puv2': PUV Method 2, 'puv3': PUV Method 3 
    coordinatesys='xyz'
        | Define the coordinate system 
        | 'xyz': XYZ coordinate system, 'enu': ENU (East North Up) coordinate system 
        | If coordinatesys='enu', then x is East and y is North  
        | If coordinatesys='enu', results are reported with respect to true north  
        | In true north coordinate system, wave comes from as:
        | 0 degree: from north, 90 degree: from east, 180 degree: from south, 270 degree: from west  
    fmaxpcorr=fs/2
        | Maximum frequency that a pressure attenuation factor applies up on that (Hz)
        | If fmaxpcorrCalcMethod='user', then the smaller of calculated and user defined fmaxpcorr will be chosen
    fminpcorr=0
        | Minimum frequency that is used for defining fmaxpcorr if fmaxpcorrCalcMethod='auto' (Hz)
        | fminpcorr should be smaller than fp 
        | If swell energy exists, fminpcorr should be smaller than fp of wind sea (fpsea) and larger than fp of swell (fpswell) if there swell 
    fKuvmin=fs/2
        Frequency that a velocity conversion factor (Kuv) at that frequency is considered as a minimum limit for Kuv
    fcL=0
        Low cut-off frequency, between 0*fs to 0.5*fs (Hz)
    fcH=fs/2
        High cut-off frequency, between 0*fs to 0.5*fs (Hz)
    fmaxpcorrCalcMethod='auto'
        | Define if to calculate fmaxpcorr and ftail or to use user defined
        | 'user': use user defined value for fmaxpcorr
        | 'auto': automatically define value for fmaxpcorr
    Kpafterfmaxpcorr='constant'
        | Define a apressure response factor, Kp, value for frequency larger than fmaxpcorr
        | 'nochange': Kp is not changed for frequency larger than fKuvmin 
        | 'one': Kp=1 for frequency larger than fmaxpcorr 
        | 'constant': Kp for f larger than fmaxpcorr stays equal to Kp at fmaxpcorr (constant)
    KuvafterfKuvmin='constant'
        | Define conversion factor, Kuv, value for frequency larger than fKuvmin
        | 'nochange': Kuv is not changed for frequency larger than fKuvmin 
        | 'one': Kuv=1 for frequency larger than fKuvmin 
        | 'constant': Kuv for f larger than fKuvmin stays equal to Kuv at fKuvmin (constant)
    kCalcMethod='beji'
        | Wave number calculation method 
        | 'hunt': Hunt (1979), 'beji': Beji (2013), 'vatankhah': Vatankhah and Aghashariatmadari (2013) 
        | 'goda': Goda (2010), 'exact': calculate exact value 
    Rho=1000
        Water density (kg/m^3)
    nfft=length(P)
        Total number of points between 0 and fs that spectrum reports at is (nfft+1)
    SegmentSize=256
        Segment size, data are divided into the segments each has a total element equal to SegmentSize
    OverlapSize=128
        Number of data points that are overlaped with data in previous segments 
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    Wavedir
        Mean wave direction (Degree)
    theta1
        Mean wave direction as a function of frequency (Degree)
    theta2
        Principal wave direction as a function of frequency (Degree)
    f
        Frequency (Hz)

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
        Ux=(np.pi/5)*(2*Eta)*(np.cosh(k*hfrombed)/np.sinh(k*h)) 
        Uy=0.2*Ux
        Wavedir,theta1,theta2,f=sm.wavedirpuv(P,Ux,Uy,fs,h,4,4,'puv1','xyz',0.7,0,0.7,0,fs/2,'auto','constant','constant','beji',1025,N,256,128,'yes')

    References
    ----------

    Beji, S. (2013). 
    Improved explicit approximation of linear dispersion relationship for gravity waves. 
    Coastal Engineering, 73, 11-12.

    Deo, M. C., Gondane, D. S., & Sanil Kumar, V. (2002). 
    Analysis of wave directional spreading using neural networks. 
    Journal of waterway, port, coastal, and ocean engineering, 128(1), 30-37.

    Earle, M. D., McGehee, D., & Tubman, M. (1995). 
    Field Wave Gaging Program, Wave Data Analysis Standard (No. WES/IR/CERC-95-2). 
    ARMY ENGINEER WATERWAYS EXPERIMENT STATION VICKSBURG MS.

    Ewans, K. C. (1998). 
    Observations of the directional spectrum of fetch-limited waves. 
    Journal of Physical Oceanography, 28(3), 495-512.

    Goda, Y. (2010). 
    Random seas and design of maritime structures. 
    World scientific.

    Grosskopf, W., Aubrey, D., Mattie, M., & Mathiesen, M. (1983). 
    Field intercomparison of nearshore directional wave sensors. 
    IEEE Journal of Oceanic Engineering, 8(4), 254-271.

    Herbers, T. H. C., Elgar, S., & Guza, R. T. (1999). 
    Directional spreading of waves in the nearshore. 
    Journal of Geophysical Research: Oceans, 104(C4), 7683-7693.

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
    Ux=type2numpy(Ux)
    Uy=type2numpy(Uy)
    h=type2numpy(h)

    #--------------------------------------------------------------------------
    #Assign default values

    if fmaxpcorr is None: fmaxpcorr=fs/2
    if fKuvmin is None: fKuvmin=fs/2
    if fcH is None: fcH=fs/2
    if nfft is None: nfft=len(P)

    #--------------------------------------------------------------------------

    #Deterending input data
    PDetrended=sp.signal.detrend(P)
    EtaDetrended=PDetrended/(Rho*9.81) #Converting pressure data to water elevation (depth) data, P=Rho.g.h
    UxDetrended=sp.signal.detrend(Ux)
    UyDetrended=sp.signal.detrend(Uy)
    
    #if (h<=0): h=0.001
    if (fmaxpcorr>fs/2): fmaxpcorr=int(fs/2)
    if (fKuvmin>fs/2): fKuvmin=int(fs/2)
    if (fcL<0): fcL=0
    if (fcH>fs/2): fcH=int(fs/2)
    #nfft=int(2**(np.ceil(np.log2(len(EtaDetrended)))))

    #--------------------------------------------------------------------------
    #Calculating power spectral density

    #Power spectral density and frequency from Welch method
    f,Syy=sp.signal.welch(EtaDetrended,fs=fs,window='hamming',nperseg=SegmentSize,noverlap=OverlapSize,nfft=nfft) 
    f,Suu=sp.signal.welch(UxDetrended,fs=fs,window='hamming',nperseg=SegmentSize,noverlap=OverlapSize,nfft=nfft) 
    f,Svv=sp.signal.welch(UyDetrended,fs=fs,window='hamming',nperseg=SegmentSize,noverlap=OverlapSize,nfft=nfft) 

    #Cross power spectral density and frequency from Welch method
    f,Syu=sp.signal.csd(EtaDetrended,UxDetrended,fs=fs,window='hamming',nperseg=SegmentSize,noverlap=OverlapSize,nfft=nfft)
    f,Syv=sp.signal.csd(EtaDetrended,UyDetrended,fs=fs,window='hamming',nperseg=SegmentSize,noverlap=OverlapSize,nfft=nfft)
    f,Suv=sp.signal.csd(UxDetrended,UyDetrended,fs=fs,window='hamming',nperseg=SegmentSize,noverlap=OverlapSize,nfft=nfft)

    df=f[1]-f[0] #Frequency interval 

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
    Kp=np.cosh(k*Pheightfrombed)/np.cosh(k*h)

    #Calculating KmaxL based on linear wave theory
    kmaxL=np.pi/(h-Pheightfrombed) #Wave number associated with fmaxpcorrL
    
    #Automatically estimating fmaxpcorr and ftailcorrection
    if fmaxpcorrCalcMethod=='auto':
    
        fminpcorrIndx=int(np.max((np.nonzero(f<=fminpcorr))[0])) #Locating the index of fminpcorr (fmaxpcorr should be larger than fminpcorr)
        SyymaxIndx=np.argmax(Syy[fminpcorrIndx:]) #Locating the peak frequency, fp, of original spectrum before applying Kp
    
        fmaxpcorrL=1/(2*np.pi)*np.sqrt(9.81*kmaxL*np.tanh(kmaxL*h)) #Maximum frequency that Kp can be applied, calculated from linear wave theory
        fmaxpcorrLIndx=int(np.max((np.nonzero(f<=fmaxpcorrL))[0])) #Locating the index of fmaxpcorrL
        if (fmaxpcorrLIndx<fminpcorrIndx+(SyymaxIndx)): fmaxpcorrLIndx=fminpcorrIndx+(SyymaxIndx-1) #Check if fmaxpcorrLIndx locataed after fp
    
        Syy1=Syy/(Kp**2) #Applying Kp on spectrum causing an increase to infinity at the tail of spectrum
        SyyminIndx=np.argmin(Syy1[fminpcorrIndx+(SyymaxIndx):fmaxpcorrLIndx+1]) #Locating the index of minimum value for Syy between fp and fmaxpcorrL
    
        fmaxpcorr1=f[fminpcorrIndx+(SyymaxIndx)+(SyyminIndx)] #Asigning the frequency of the minimum value of Syy between fp and fmaxpcorrL
        if (fmaxpcorr1>fmaxpcorrL): fmaxpcorr1=fmaxpcorrL #Check fmaxpcorr1 be smaller than fmaxpcorrL
        if ((fmaxpcorr1==f[fminpcorrIndx+(SyymaxIndx)]) and (fmaxpcorrL>f[fminpcorrIndx+(SyymaxIndx)])): fmaxpcorr1=fmaxpcorrL #If fmaxpcorrL>fp then fmaxpcorr1 should not be equal to fp
        if (fmaxpcorr>fmaxpcorr1): fmaxpcorr=fmaxpcorr1
    
        #ftail=f[fminpcorrIndx+(SyymaxIndx-1)+(SyyminIndx-1)] #Asigning the frequency of the minimum value of Syy between fp and fmaxpcorrL for ftail
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
        Indx1=int(np.max((np.nonzero(f<=fmaxpcorr-0.05))[0]))
        Indx2=int(np.max((np.nonzero(f<=fmaxpcorr+0.05))[0]))
        if (Indx2>len(f)): Indx2=len(f)
        for i in range(Indx1,Indx2+1,1):
            Kp[i]=(Kp[Indx2]-Kp[Indx1])/(Indx2-Indx1)*(i-Indx1)+Kp[Indx1]

    
    elif Kpafterfmaxpcorr=='constant':
        fmaxpcorrIndx=int(np.max((np.nonzero(f<=fmaxpcorr))[0]))
        if (fmaxpcorrIndx>len(f)): fmaxpcorrIndx=len(f)
        Kp[f>fmaxpcorr]=Kp[fmaxpcorrIndx] #Kp for f larger than fmaxpcorr stays equal to Kp at fmaxpcorr (constant)

    #Comparing Kp with KpminL calculated based on linear wave theory
    KpminL=np.cosh(kmaxL*Pheightfrombed)/np.cosh(kmaxL*h) #Minimum Limit for Kp calculated based on linear wave theory
    Kp[Kp<KpminL]=KpminL #Check to avoid large amplification, Kp should be larger than minimum Kp calculated based on linear wave theory

    #--------------------------------------------------------------------------
    #Pressure attenuation correction

    SyyBeforeKp=Syy.copy() #Power spectral density before applying Kp
    
    #Applying Kp on power spectral density
    Syy= Syy/(Kp**2) #Applies pressure response factor, Kp

    #--------------------------------------------------------------------------
    #Calculating the velocity response factor, Kuv, (Wiberg, 2008)
    Kuv=(2*np.pi*f)*np.cosh(k*UVheightfrombed)/np.sinh(k*h)
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
    kmaxL=np.pi/(h-UVheightfrombed) #Wave number associated with fmaxuvcorrL
    
    #Calculating fmaxuvcorrL based on linear wave theory
    fmaxuvcorrL=1/(2*np.pi)*np.sqrt(9.81*np.pi/(h-UVheightfrombed)*np.tanh(np.pi*h/(h-UVheightfrombed))) #fmaxuvcorrL based on linear wave theory
    
    #Comparing Kuv with KuvminL calculated based on linear wave theory
    KuvminL=(2*np.pi*fmaxuvcorrL)*np.cosh(kmaxL*UVheightfrombed)/np.sinh(kmaxL*h) #Minimum Limit for Kuv calculated based on linear wave theory
    Kuv[Kuv<KuvminL]=KuvminL #Check to avoid large amplification, Kuv should be larger than minimum Kuv calculated based on linear wave theory

    #--------------------------------------------------------------------------
    #Applying Kuv and Kp on power spectral density

    Suu=Suu/(Kuv**2)
    Svv=Svv/(Kuv**2)
    Suv=Suv/(Kuv**2)
    Syu=Syu/(Kp*Kuv)
    Syv=Syv/(Kp*Kuv)

    #--------------------------------------------------------------------------
    #Checking power spectral density for negative values

    Syy[Syy<0]=0 #Syy can not be negative
    Suu[Syy<0]=0 #Suu can not be negative
    Svv[Syy<0]=0 #Svv can not be negative
    Syu[Syy<0]=0 #Syu can not be negative
    Syv[Syy<0]=0 #Syv can not be negative
    Suv[Syy<0]=0 #Suv can not be negative

    #--------------------------------------------------------------------------
    #Cut off spectrum based on fcL and fcH

    if ((fcL>f[0]) and (fcL<fs/2)):
        Indx=np.int_((np.nonzero(f<fcL))[0])
        Syy=np.delete(Syy,Indx)
        Suu=np.delete(Suu,Indx)
        Svv=np.delete(Svv,Indx)
        Syu=np.delete(Syu,Indx)
        Syv=np.delete(Syv,Indx)
        Suv=np.delete(Suv,Indx)
        f=np.delete(f,Indx)

    if ((fcH>f[0]) and (fcH<fs/2)):
        Indx=np.int_((np.nonzero(f>fcH))[0])
        Syy=np.delete(Syy,Indx)
        Suu=np.delete(Suu,Indx)
        Svv=np.delete(Svv,Indx)
        Syu=np.delete(Syu,Indx)
        Syv=np.delete(Syv,Indx)
        Suv=np.delete(Suv,Indx)
        f=np.delete(f,Indx)

    #--------------------------------------------------------------------------
    #Calculating wave properties

    #Calculating spectral moments
    m0=np.sum(Syy*f**0*df)
    m1=np.sum(Syy*f**1*df)
    m2=np.sum(Syy*f**2*df)

    #Calculating wave properties
    Hm0=4*np.sqrt(m0) #Zero-Moment wave height
    Tm01=m0/m1 #mean period
    Tm02=(m0/m2)**0.5 #zero crossing period

    #Calculating peak period
    SyymaxIndx=np.argmax(Syy) #Locating an index of the spectrum peak
    fp=f[SyymaxIndx] #peak frequency
    Tp=1/fp #peak period

    #Calculating peak frequency from weighted integral (Young, 1995)
    #fp=(np.sum(Syy**5*f**1*df))/(np.sum(Syy**5*f**0*df)) #Peak frequency

    #--------------------------------------------------------------------------
    #Calculating Fourier coefficients
    #(e.g. Grosskopf (1983); Earle (1995); Ewans (1998); Herbers et al. (1999); Deo et al. (2002))

    if dirCalcMethod=='puv1':
    
        #Note: power spectra are already devided by Kp and Kuv
        #Note: Kuv already multiplied by (2*pi*f)
        a0=np.real(Syy/(2*np.pi))
        a1=np.real(Syu/(Syy*(Suu+Svv))**0.5)
        b1=np.real(Syv/(Syy*(Suu+Svv))**0.5)
        a2=np.real((Suu-Svv)/(Suu+Svv))
        b2=np.real(2*Suv/(Suu+Svv))
    
    elif dirCalcMethod=='puv2':
    
        #Note: Power spectra are already devided by Kp and Kuv
        #Note: Kuv already multiplied by (2*pi*f)
        a0=np.real(Syy/(2*np.pi))  
        a1=np.real(Syu/(np.pi))  
        b1=np.real(Syv/(np.pi))
        a2=np.real((Suu-Svv)/(np.pi))
        b2=np.real(2*Suv/(np.pi))
    
    
    #Checking for NaN and Inf
    if ((dirCalcMethod=='puv1') or (dirCalcMethod=='puv2')): 
    
        a0[((np.isnan(a0)==1) | (np.isinf(a0)==1))]=0
        a1[((np.isnan(a1)==1) | (np.isinf(a1)==1))]=0
        b1[((np.isnan(b1)==1) | (np.isinf(b1)==1))]=0
        a2[((np.isnan(a2)==1) | (np.isinf(a2)==1))]=0
        b2[((np.isnan(b2)==1) | (np.isinf(b2)==1))]=0
    
    #--------------------------------------------------------------------------
    #Calculating mean wave direction

    if ((dirCalcMethod=='puv1') or (dirCalcMethod=='puv2')): 
    
        theta1=np.arctan2(np.real(b1),np.real(a1)) #Mean wave direction as a function of frequency
        theta2=0.5*np.arctan2(np.real(b2),np.real(a2)) #Principal wave direction as a function of frequency
        #theta_m_f=np.arctan2(np.real(b1),np.real(a1)) #Mean wave direction as a function of frequency
        #theta_m_rad=np.arctan2(np.sum(np.sin(theta_m_f)*(Syy)*df),np.sum(np.cos(theta_m_f)*(Syy)*df)) #Mean wave direction in radian, Note: Syy already devided by (Kp)
        theta_m_rad=np.arctan2(np.sum(np.sin(theta1)*(Syy)*df),np.sum(np.cos(theta1)*(Syy)*df)) #Mean wave direction in radian, Note: Syy already devided by (Kp)
        Wavedir=np.rad2deg(theta_m_rad) #Mean wave direction in degree
    
    elif dirCalcMethod=='puv3': #Calculating mean wave direction using FFT
    
        #Calculating Fast Fourier transform
        EtaFFT=np.fft.fft(EtaDetrended)
        UxFFT=np.fft.fft(UxDetrended)
        UyFFT=np.fft.fft(UyDetrended)
    
        #Generate symmetrical Kp and Kuv respect with fs/2
        KpSym=np.concatenate((Kp,np.flipud(Kp[1:-1])))
        KuvSym=np.concatenate((Kuv,np.flipud(Kuv[1:-1])))
    
        #Applying Kp and Kuv
        EtaFFT=EtaFFT/(KpSym)
        UxFFT=UxFFT/KuvSym
        UyFFT=UyFFT/KuvSym
    
        #Inverse Fast Fourier transform
        EtaIFFT=np.real(np.fft.ifft(EtaFFT))	
        UxIFFT=np.real(np.fft.ifft(UxFFT))	
        UyIFFT=np.real(np.fft.ifft(UyFFT))
    
        #Choosing data related to wave motion direction
        UxIFFT=UxIFFT[EtaDetrended>0]
        UyIFFT=UyIFFT[EtaDetrended>0]
        EtaIFFT=EtaIFFT[EtaDetrended>0]
    
        #Calculating mean wave direction
        theta1=np.arctan2(np.abs(UyFFT[0:len(f)]),np.abs(UxFFT[0:len(f)])) #Mean wave direction as a function of frequency
        theta2=np.nan #Principal wave direction as a function of frequency
        theta_rad=np.arctan2(UyIFFT,UxIFFT) #Time series of wave direction in radian
        theta_m_rad=(np.arctan2(np.sum(np.sin(theta_rad)*EtaIFFT**1),np.sum(np.cos(theta_rad)*EtaIFFT**1))) #Mean wave direction in radian
        Wavedir=np.rad2deg(theta_m_rad) #Mean wave direction in degree
    

    #--------------------------------------------------------------------------
    #Converting mathematical angle direction to compass direction with respect to true north

    if coordinatesys=='enu':
        theta1=270-theta1 #Converting mathematical direction to compass direction
        theta1[theta1<0]=theta1[theta1<0]+360
    
        theta2=270-theta2 #Converting mathematical direction to compass direction
        theta2[theta2<0]=theta2[theta2<0]+360
    
        Wavedir=270-Wavedir #Converting mathematical direction to compass direction
        if (Wavedir<0): Wavedir=Wavedir+360
        if (Wavedir>360): Wavedir=Wavedir-360

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':

        val=[Hm0,fp,Tp,Tm01,Tm02,m0,m1,m2]
        name=['Hm0','fp','Tp','Tm01','Tm02','m0','m1','m2']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            print('{0:10}= {1:0.10f}'.format(name[i],val[i])) 

        #Plotting
        def cart2pol(x,y):
            rho=np.sqrt(x**2+y**2)
            phi = np.arctan2(y,x)
            return(phi,rho)

        def pol2cart(phi,rho):
            x=rho*np.cos(phi)
            y=rho*np.sin(phi)
            return(x,y)

        ax=plt.subplot(111, projection='polar')
        #Wavedirx,Wavediry=pol2cart(np.deg2rad(Wavedir),1)
        #arr1=ax.plot([np.deg2rad(Wavedir),np.deg2rad(Wavedir)],[0,1])
        arr1=ax.arrow(np.deg2rad(Wavedir),0,0,1)
    
        plt.title('Mean Wave Direction')


    #--------------------------------------------------------------------------
    #Outputs
    return Wavedir, theta1, theta2, f

    #--------------------------------------------------------------------------
