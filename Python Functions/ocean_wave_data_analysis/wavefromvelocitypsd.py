def wavefromvelocitypsd(Ux, Uy, fs, h, heightfrombed=0, fKuvmin=None, fcL=0, fcH=None, KuvafterfKuvmin='constant', kCalcMethod='beji', nfft=None, SegmentSize=256, OverlapSize=128, dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2017-04-0                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.wavefromvelocitypsd
    ===============================

    .. code:: python

        Hm0, fp, Tp, Tm01, Tm02, f, Syy, m0 = scientimate.wavefromvelocitypsd(Ux, Uy, fs, h, heightfrombed=0, fKuvmin=None, fcL=0, fcH=None, KuvafterfKuvmin='constant', kCalcMethod='beji', nfft=None, SegmentSize=256, OverlapSize=128, dispout='no')

    Description
    -----------

    Calculate wave properties from wave orbital velocity by converting it to water surface elevation power spectral density

    Inputs
    ------

    Ux
        Wave horizontal orbital velocity data in x direction in (m/s)
    Uy
        Wave horizontal orbital velocity data in y direction in (m/s)
    fs
        Sampling frequency that data collected at in (Hz)
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
    nfft=length(Eta)
        Total number of points between 0 and fs that spectrum reports at is (nfft+1)
    SegmentSize=256
        Segment size, data are divided into the segments each has a total element equal to SegmentSize
    OverlapSize=128
        Number of data points that are overlapped with data in previous segments 
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    Hm0
        Zero-Moment Wave Height (m)
    fp
        Peak wave frequency (Hz)
    Tp
        Peak wave period (second)
    Tm01
        Wave Period from m01 (second), Mean Wave Period
    Tm02
        Wave Period from m02 (second), Mean Zero Crossing Period
    f
        Frequency (Hz)
    Syy
        Power spectral density (m^2/Hz)
    m0
        Zero-Moment of the power spectral density (m^2)

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
        Ux=(np.pi/5)*(2*Eta)*(np.cosh(k*hfrombed)/np.sinh(k*h)) 
        Uy=0.2*Ux
        Hm0,fp,Tp,Tm01,Tm02,f,Syy,m0=sm.wavefromvelocitypsd(Ux,Uy,fs,5,4,0.6,0,fs/2,'constant','beji',N,256,128,'yes')

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
    
    Ux=type2numpy(Ux)
    Uy=type2numpy(Uy)
    h=type2numpy(h)

    #--------------------------------------------------------------------------
    #Assign default values

    if fKuvmin is None: fKuvmin=fs/2
    if fcH is None: fcH=fs/2
    if nfft is None: nfft=len(P)

    #--------------------------------------------------------------------------
    #Preventing error in scipy signal by making sure it is a single column array
    #Ux=Ux[:,0]
    #Uy=Uy[:,0]

    #--------------------------------------------------------------------------

    #Deterending input data
    UxDetrended=sp.signal.detrend(Ux)
    UyDetrended=sp.signal.detrend(Uy)
    #if (h<=0): h=0.001
    if (fKuvmin>fs/2): fKuvmin=int(fs/2)
    if (fcL<0): fcL=0
    if (fcH>fs/2): fcH=int(fs/2)
    #nfft=int(2**(np.ceil(np.log2(len(UxDetrended)))))

    #--------------------------------------------------------------------------
    #Calculating power spectral density

    #Velocity power spectral density and frequency from Welch method
    f,Suu=sp.signal.welch(UxDetrended,fs,window='hamming',nperseg=SegmentSize,noverlap=OverlapSize,nfft=nfft) 
    f,Svv=sp.signal.welch(UyDetrended,fs,window='hamming',nperseg=SegmentSize,noverlap=OverlapSize,nfft=nfft) 
    
    Suv=Suu+Svv #Combined horizontal spectrum (Wiberg, 2008)

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

    #--------------------------------------------------------------------------
    #Converting velocity power spectral density to water level power spectral density

    Syy=Suv/(Kuv**2) #Calculating water level spectrum from orbital velocity spectrum

    #--------------------------------------------------------------------------
    #Cut off spectrum based on fcL and fcH

    if ((fcL>f[0]) and (fcL<fs/2)):
        #Indx=np.int_((np.nonzero(f<fcL))[0])
        #Syy=np.delete(Syy,Indx)
        #f=np.delete(f,Indx)
        Syy[f<fcL]=0

    if ((fcH>f[0]) and (fcH<fs/2)):
        #Indx=np.int_((np.nonzero(f>fcH))[0])
        #Syy=np.delete(Syy,Indx)
        #f=np.delete(f,Indx)
        Syy[f>fcH]=0

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
    #Displaying results

    if dispout=='yes':

        val=[Hm0,fp,Tp,Tm01,Tm02,m0,m1,m2]
        name=['Hm0','fp','Tp','Tm01','Tm02','m0','m1','m2']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            print('{0:10}= {1:0.10f}'.format(name[i],val[i])) 

        #Plotting
        plt.loglog(f[f!=0],Syy[f!=0])
    
        plt.title('Power Spectral Density')
        plt.xlabel('Frequency(Hz)')
        plt.ylabel('Spectral Density(m^2/Hz)')

    #--------------------------------------------------------------------------
    #Outputs
    return Hm0, fp, Tp, Tm01, Tm02, f, Syy, m0

    #--------------------------------------------------------------------------
