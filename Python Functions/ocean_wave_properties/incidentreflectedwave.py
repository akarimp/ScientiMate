def incidentreflectedwave(Eta1, Eta2, dx, h, fs=2, fmin=0, fmax=None, SepMethod='goda', kCalcMethod='beji', dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2017-02-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.incidentreflectedwave
    =================================

    .. code:: python

    Kr, EtaInc, EtaRef, aInc, aRef, t, f = scientimate.incidentreflectedwave(Eta1, Eta2, dx, h, fs=2, fmin=0, fmax=None, SepMethod='goda', kCalcMethod='beji', dispout='no')

    Description
    -----------

    Separate incident and reflected waves

    Inputs
    ------

    Eta1
        Wave signal at position x1 in (m), should have same size as Eta2  
    Eta2
        Wave signal at position x2 in (m), should have same size as Eta1 
    dx
        Distance between x1 and x2 in (m), dx=x2-x1 
    h
        Mean water depth in (m) 
    fs=2
        Sampling frequency that data collected at in (Hz), if fs=1 then output is equivalent to normalized filter
    fmin=0
        Minimum frequency to be considered in (Hz)
    fmax=fs/2
        Maximum frequency to be considered in (Hz)
    SepMethod='goda'
        | Incident and reflected waves Separation method 
        | 'goda': Goda and Suzuki (1977) 
        | 'ma': Ma et al. (2010)
        | 'frigaard': Frigaard and Brorsen (1995) 
    kCalcMethod='beji'
        | Wave number calculation method 
        | 'hunt': Hunt (1979), 'beji': Beji (2013), 'vatankhah': Vatankhah and Aghashariatmadari (2013) 
        | 'goda': Goda (2010), 'exact': calculate exact value 
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    Kr
        Reflection coefficient
    EtaInc
        Incident wave (m)
    EtaRef
        Reflected wave (m)
    aInc
        Amplitude of incident wave (m)
    aRef
        Amplitude of reflected wave (m)
    t
        Time (s)
    f
        Frequency (Hz)

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        h=1
        fs=2
        dt=1/fs
        duration=1024
        t=np.linspace(0,duration-dt,duration*fs)
        x1=1 
        x2=3.7
        dx=x2-x1
        W1=0.5*np.cos(0.412*x1-2*np.pi*0.2*t)
        W2=0.5*np.cos(0.739*x1-2*np.pi*0.34*t)
        EtaIncGauge1=(W1+W2)
        W3=0.5*np.cos(0.412*x2-2*np.pi*0.2*t)
        W4=0.5*np.cos(0.739*x2-2*np.pi*0.34*t)
        EtaIncGauge2=(W3+W4)
        W5=0.1*np.cos(0.412*x1+2*np.pi*0.2*t)
        W6=0.1*np.cos(0.739*x1+2*np.pi*0.34*t)
        EtaRefGauge1=(W5+W6)
        W7=0.1*np.cos(0.412*x2+2*np.pi*0.2*t)
        W8=0.1*np.cos(0.739*x2+2*np.pi*0.34*t)
        EtaRefGauge2=(W7+W8)
        Eta1=EtaIncGauge1+EtaRefGauge1
        Eta2=EtaIncGauge2+EtaRefGauge2

        Kr,EtaInc,EtaRef,aInc,aRef,t,f=sm.incidentreflectedwave(Eta1,Eta2,dx,h,fs,0,fs/2,'goda','beji','yes')

    References
    ----------

    Beji, S. (2013). 
    Improved explicit approximation of linear dispersion relationship for gravity waves. 
    Coastal Engineering, 73, 11-12.

    Baldock, T. E., & Simmonds, D. J. (1999). Separation of incident and reflected waves over sloping bathymetry. 
    Coastal Engineering, 38(3), 167-176.

    Frigaard, P., Brorsen, M., 1995. A time domain method for separating incident and reflected irregular waves.
    Coastal Eng. 24, 205–215

    Goda, Y., & Suzuki, Y. (1977). Estimation of incident and reflected waves in random wave experiments. 
    In Coastal Engineering 1976 (pp. 828-845).

    Goda, Y. (2010). 
    Random seas and design of maritime structures. 
    World scientific.

    Hunt, J. N. (1979). 
    Direct solution of wave dispersion equation. 
    Journal of the Waterway Port Coastal and Ocean Division, 105(4), 457-459.

    Ma, Y., Dong, G., Ma, X., & Wang, G. (2010). 
    A new method for separation of 2D incident and reflected waves by the Morlet wavelet transform. 
    Coastal Engineering, 57(6), 597-603.

    Mansard, E. P., & Funke, E. R. (1980). 
    The measurement of incident and reflected spectra using a least squares method. 
    In Coastal Engineering 1980 (pp. 154-172).

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
    
    h=type2numpy(h)
    Eta1=type2numpy(Eta1)
    Eta2=type2numpy(Eta2)

    #--------------------------------------------------------------------------
    #Assign default values

    if fmax is None: fmax=fs/2

    #--------------------------------------------------------------------------

    #Detrending data
    Eta1=sp.signal.detrend(Eta1)
    Eta2=sp.signal.detrend(Eta2)

    #Generating frequency vector
    N=np.min((len(Eta1),len(Eta2))) #Length of a time series

    #Make Eta1 and Eta2 the same size
    Eta1=Eta1[0:N]
    Eta2=Eta2[0:N]

    df=fs/N #Frequency difference between consecutive samples, df=fs/N
    f=np.arange(0,fs/2+df,df) #Frequency vector with (N/2+1) elements from 0 Hz to fNy=fs/2 Hz, equally spaced at df
    # Note: in general f(:,1)=[0:N-1]*df from 0 Hz to (fs-df) Hz, equally spaced at df

    dt=1/fs
    t=np.linspace(0,(N-1)*dt,N) #Time from 0 to T-dt, equally spaced at dt

    #Check the frequency be between fmain and fmax and not ne zero
    if fmin<0: fmin=0
    if fmax>fs/2: fmax=fs/2

    #--------------------------------------------------------------------------
    #Calculating wave number (k)

    w=2*np.pi*f #Wave angular frequency, w=2.pi/Tp
    
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
    #Calculating Incident and Reflected waves

    #Method 2: Goda, Y., & Suzuki, Y. (1977)
    if SepMethod=='goda':

        #Fast Fourier Transform of the first signal
        #Eta1FFT=np.fft.fft(Eta1) #Fourier transform of the first signal
        Eta1FFT=np.fft.fft(Eta1)*(2/N) #Scaling fft 
        Eta1FFT=Eta1FFT[0:len(f)] #(N/2+1)elements if N is even and (N+1)/2 elements if N is odd
        Eta1amp=np.abs(Eta1FFT) #Amplitude, is equal to sqrt((real(fft))**2+(imag(fft))**2)
        A1=np.real(Eta1FFT) #Real pert of the Fourier transform
        B1=np.imag(Eta1FFT) #Imaginary pert of the Fourier transform
        Eta1Phi=np.angle(Eta1FFT) #Wave phase

        #Fast Fourier Transform of the second signal
        #Eta2FFT=np.fft.fft(Eta2) #Fourier transform of the second signal
        Eta2FFT=np.fft.fft(Eta2)*(2/N) #Scaling fft 
        Eta2FFT=Eta2FFT[0:len(f)] #(N/2+1)elements if N is even and (N+1)/2 elements if N is odd
        Eta2amp=np.abs(Eta2FFT) #Amplitude, is equal to sqrt((real(fft))**2+(imag(fft))**2)
        A2=np.real(Eta2FFT) #Real pert of the Fourier transform
        B2=np.imag(Eta2FFT) #Imaginary pert of the Fourier transform
        Eta2Phi=np.angle(Eta2FFT) #Wave phase

        aInc=np.zeros(len(f)) #Pre-assigning array to make program run faster
        aRef=np.zeros(len(f)) #Pre-assigning array to make program run faster
        C=1/(2*np.abs(np.sin(k*dx)))
        C[np.abs(C)>1]=1
        C[0]=1
        #aInc=C*np.sqrt((A2-A1*np.cos(k*dx)-B1*np.sin(k*dx))**2+(B2+A1*np.sin(k*dx)-B1*np.cos(k*dx))**2) #Incident amplitude
        #aRef=C*np.sqrt((A2-A1*np.cos(k*dx)+B1*np.sin(k*dx))**2+(B2-A1*np.sin(k*dx)-B1*np.cos(k*dx))**2) #Reflected amplitude
        aInc=C*np.sqrt((A2-A1*np.cos(k*dx)+B1*np.sin(k*dx))**2+(B2-A1*np.sin(k*dx)-B1*np.cos(k*dx))**2) #Incident amplitude
        aRef=C*np.sqrt((A2-A1*np.cos(k*dx)-B1*np.sin(k*dx))**2+(B2+A1*np.sin(k*dx)-B1*np.cos(k*dx))**2) #Reflected amplitude

        C=1/(2*np.sqrt(-1+0j)*np.sin(k*dx))
        C[np.abs(C)>1]=1
        C[0]=1
        #aIncFFT=C*((Eta1FFT)*np.exp(np.sqrt(-1+0j)*k*dx)-(Eta2FFT)) #Incident amplitude
        #aRefFFT=C*((Eta1FFT)*np.exp(-np.sqrt(-1+0j)*k*dx)-(Eta2FFT)) #Reflected amplitude
        aIncFFT=C*(np.abs(Eta1FFT)*np.exp(np.sqrt(-1+0j)*(np.angle(Eta1FFT))+np.sqrt(-1+0j)*(k*dx))-(Eta2FFT)) #Incident amplitude
        aRefFFT=C*(np.abs(Eta1FFT)*np.exp(np.sqrt(-1+0j)*(np.angle(Eta1FFT))-np.sqrt(-1+0j)*(k*dx))-(Eta2FFT)) #Reflected amplitude
        aIncPhi=-np.angle(aIncFFT) #Wave phase
        aRefPhi=-np.angle(aRefFFT)+np.pi #Wave phase

    #Method 1: Ma, Y., Dong, G., Ma, X., & Wang, G. (2010)
    elif SepMethod=='ma':

        #Fast Fourier Transform of the first signal
        #Eta1FFT=np.fft.fft(Eta1) #Fourier transform of the first signal
        Eta1FFT=np.fft.fft(Eta1)*(2/N) #Scaling fft 
        Eta1FFT=Eta1FFT[0:len(f)] #(N/2+1)elements if N is even and (N+1)/2 elements if N is odd

        #Fast Fourier Transform of the second signal
        #Eta2FFT=np.fft.fft(Eta2) #Fourier transform of the second signal
        Eta2FFT=np.fft.fft(Eta2)*(2/N) #Scaling fft 
        Eta2FFT=Eta2FFT[0:len(f)] #(N/2+1)elements if N is even and (N+1)/2 elements if N is odd

        aInc=np.zeros(len(f)) #Pre-assigning array to make program run faster
        aRef=np.zeros(len(f)) #Pre-assigning array to make program run faster
        C=1/(2*np.sqrt(-1+0j)*np.sin(k*dx))
        C[np.abs(C)>1]=1
        C[0]=1
        #aIncFFT=C*((Eta1FFT)*np.exp(np.sqrt(-1+0j)*k*dx)-(Eta2FFT)) #Incident amplitude
        #aRefFFT=C*((Eta1FFT)*np.exp(-np.sqrt(-1+0j)*k*dx)-(Eta2FFT)) #Reflected amplitude
        aIncFFT=C*(np.abs(Eta1FFT)*np.exp(np.sqrt(-1+0j)*(np.angle(Eta1FFT))+np.sqrt(-1+0j)*(k*dx))-(Eta2FFT)) #Incident amplitude
        aRefFFT=C*(np.abs(Eta1FFT)*np.exp(np.sqrt(-1+0j)*(np.angle(Eta1FFT))-np.sqrt(-1+0j)*(k*dx))-(Eta2FFT)) #Reflected amplitude
        aInc=np.abs(aIncFFT)
        aRef=np.abs(aRefFFT)
        aIncPhi=-np.angle(aIncFFT) #Wave phase
        aRefPhi=-np.angle(aRefFFT)+np.pi #Wave phase

    #Method 3: frigaardaard, P., Brorsen, M. (1995)
    elif SepMethod=='frigaard':
        
        #Fast Fourier Transform of the first signal
        Eta1FFT=np.fft.fft(Eta1) #Fourier transform of the first signal
        #Eta1FFT=np.fft.fft(Eta1)*(2/N) #Scaling fft 
        Eta1FFT=Eta1FFT[0:len(f)] #(N/2+1)elements if N is even and (N+1)/2 elements if N is odd

        #Fast Fourier Transform of the second signal
        Eta2FFT=np.fft.fft(Eta2) #Fourier transform of the second signal
        #Eta2FFT=np.fft.fft(Eta2)*(2/N) #Scaling fft 
        Eta2FFT=Eta2FFT[0:len(f)] #(N/2+1)elements if N is even and (N+1)/2 elements if N is odd

        #n(:,1)=[0:1:length(f(:,1))-1] #Counter from 0 to N/2+1
        #m(:,1)=[0:1:length(f(:,1))-1] #Counter from 0 to N/2+1
        n=0
        m=0
        Phi1=k*dx+np.pi/2+m*np.pi+n*2*np.pi
        Phi2=-np.pi/2-m*np.pi+n*2*np.pi
        C=1/(2*np.cos(-k*dx-np.pi/2-m*np.pi)) #Amplification factor
        C[np.abs(C)>1]=1
        C[0]=1
        
        #Scaling FFT and adding Phi to phase
        #Eta1FFTStar=C*(Eta1FFT*np.exp(np.sqrt(-1+0j)*Phi1))
        #Eta2FFTStar=C*(Eta2FFT*np.exp(np.sqrt(-1+0j)*Phi2))
        Eta1FFTStar=C*(np.abs(Eta1FFT)*np.exp(np.sqrt(-1+0j)*(np.angle(Eta1FFT)+Phi1))) 
        Eta2FFTStar=C*(np.abs(Eta2FFT)*np.exp(np.sqrt(-1+0j)*(np.angle(Eta2FFT)+Phi2)))

        EtaFFTInc=Eta1FFTStar+Eta2FFTStar
        EtaFFTInc[(f<fmin) | (f>fmax)]=0
        aInc=np.abs(EtaFFTInc)*(2/N) #Incident amplitude

        EtaFFTRef=Eta1FFT-EtaFFTInc
        EtaFFTRef[(f<fmin) | (f>fmax)]=0
        aRef=np.abs(EtaFFTRef)*(2/N) #Reflected amplitude

        aIncPhi=-np.angle(EtaFFTInc) #Wave phase
        aRefPhi=-np.angle(EtaFFTRef) #Wave phase

        EtaFFTInc=np.concatenate((EtaFFTInc,np.flipud(EtaFFTInc[1:-1]))) #Make it symmetrical with respect to fs/2
        EtaInc1=np.real(np.fft.ifft(EtaFFTInc)) #Incident water surface

        EtaFFTRef=np.concatenate((EtaFFTRef,np.flipud(EtaFFTRef[1:-1]))) #Make it symmetrical with respect to fs/2
        EtaRef1=np.real(np.fft.ifft(EtaFFTRef)) #Reflected water surface

    #--------------------------------------------------------------------------
    #Generating an incident and reflected waves

    aInc[0]=0
    aRef[0]=0

    #Remove data that are not between fmain and fmax
    aInc[(f<fmin) | (f>fmax)]=0
    aRef[(f<fmin) | (f>fmax)]=0

    EtaijInc=np.zeros((len(t),len(aInc))) #Pre-assigning array to make program run faster
    for i in range(0,len(aInc),1):
        EtaijInc[:,i]=(aInc[i])*np.cos(-2*np.pi*f[i]*t+aIncPhi[i])

    EtaijRef=np.zeros((len(t),len(aRef))) #Pre-assigning array to make program run faster
    for i in range(0,len(aRef),1):
        EtaijRef[:,i]=(aRef[i])*np.cos(-2*np.pi*f[i]*t+aRefPhi[i])

    EtaInc=np.sum(EtaijInc,axis=1) #Incident water surface
    EtaRef=np.sum(EtaijRef,axis=1) #Reflected water surface

    #--------------------------------------------------------------------------
    #Calculating power spectral density

    #Fast Fourier Transform (FFT) of the water surface elevation from f=0 Hz to f=fs Hz
    EtaIncFFT=np.fft.fft(EtaInc)
    EtaRefFFT=np.fft.fft(EtaRef)

    #Fast Fourier Transform (FFT) of the water surface elevation from f=0 Hz to fNy=fs/2 Hz
    EtaIncFFT=EtaIncFFT[0:len(f)] #(N/2+1)elements if N is even and (N+1)/2 elements if N is odd
    EtaRefFFT=EtaRefFFT[0:len(f)] #(N/2+1)elements if N is even and (N+1)/2 elements if N is odd

    #Half of the two-sided power spectral density (psd) from f=0 Hz to fNy=fs/2 Hz
    SxxEtaInc_2Sided=(1/(N*fs))*(np.abs(EtaIncFFT))**2 #Calculating psd using fs in (m**2/Hz)
    SxxEtaRef_2Sided=(1/(N*fs))*(np.abs(EtaRefFFT))**2 #Calculating psd using fs in (m**2/Hz)

    #One-sided power spectral density (psd) from f=0 Hz to fNy=fs/2 Hz
    SxxEtaInc_1Sided=SxxEtaInc_2Sided
    SxxEtaRef_1Sided=SxxEtaRef_2Sided

    SxxEtaInc_1Sided[1:-1]=2*SxxEtaInc_1Sided[1:-1] #one-side-spectrum=2*two-sided-spectrum in (m**2/Hz)
    SxxEtaRef_1Sided[1:-1]=2*SxxEtaRef_1Sided[1:-1] #one-side-spectrum=2*two-sided-spectrum in (m**2/Hz)

    SxxEtaInc=SxxEtaInc_1Sided #One-sided power spectral density (psd)
    SxxEtaRef=SxxEtaRef_1Sided #One-sided power spectral density (psd)

    #--------------------------------------------------------------------------
    #Calculating reflection coefficient

    #Calculating spectral moments
    m0EtaInc=np.sum(SxxEtaInc*f**0*df)
    m0EtaRef=np.sum(SxxEtaRef*f**0*df)

    #Calculating wave properties
    Hm0EtaInc=4*np.sqrt(m0EtaInc) #Zero-Moment wave height
    Hm0EtaRef=4*np.sqrt(m0EtaRef) #Zero-Moment wave height

    #Calculating reflection coefficient
    Kr=Hm0EtaRef/Hm0EtaInc #Reflection coefficient

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':

        plt.subplot(3,2,1)
        plt.plot(f,aInc)
        plt.xlabel('f (Hz)')
        plt.ylabel('Amplitude (m)')
        plt.title('Incident Wave')

        plt.subplot(3,2,2)
        plt.plot(f,aRef)
        plt.xlabel('f (Hz)')
        plt.ylabel('Amplitude (m)')
        plt.title('Reflected Wave')

        plt.subplot(3,1,2)
        plt.plot(t,EtaInc)
        plt.xlabel('Time (s)')
        plt.ylabel('Water Surface (m)')
        plt.title('Incident Wave')

        plt.subplot(3,1,3)
        plt.plot(t,EtaRef)
        plt.xlabel('Time (s)')
        plt.ylabel('Water Surface (m)')
        plt.title('Reflected Wave')

    #--------------------------------------------------------------------------
    #Outputs
    return Kr, EtaInc, EtaRef, aInc, aRef, t, f

    #--------------------------------------------------------------------------
