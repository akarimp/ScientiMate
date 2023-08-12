def spectrogrampsd(x, fs=2, SegmentSize=256, OverlapSize=0, WindowName='hamming', nfft=256, outtype='psd', OutputSmoothSize=0, dispout='no'):
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

    scientimate.spectrogrampsd
    =======================

    .. code:: python

        f, t, Sxx = scientimate.spectrogrampsd(x, fs=2, SegmentSize=256, OverlapSize=0, WindowName='hamming', nfft=256, outtype='psd', OutputSmoothSize=0, dispout='no')

    Description
    -----------

    Calculate spectrogram following Welch's method without averaging

    Inputs
    ------

    x
        Input data
    fs=2
        Sampling frequency that data collected at in (Hz)
    SegmentSize=256
        Segment size, data are divided into the segments each has a total element equal to SegmentSize
    OverlapSize=0
        | Number of data points that are overlaped with data in previous segments 
        | OverlapSize is recomneded to be half of the SegmentSize
    WindowName='hamming'
        | Window name, define if multiplying input data by window function or not ('none': not multiplying)
        | 'none','rectangular','triangular','welch','hanning','hamming','gaussian','blackman','nuttall','blackmanharris'
    nfft=length(x)
        Total number of points between 0 and fs that spectrum reports at is (nfft+1)
    outtype='psd'
        | Define output type
        | 'psd': power spectral density, 'db': decibel   
    OutputSmoothSize=0
        Window size for smoothing calculated spectrum (0, 1 or 2: not smoothing, reports original Welch spectrum)
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    f
        Frequency in (Hz)
    Sxx
        Spectrogram (m^2/Hz) or (dB/Hz)
    t
        Time at midpoint of each section (s)

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np
        import scipy as sp
        from scipy import signal

        x=sp.signal.detrend(0.5*np.cos(2*np.pi*0.2*np.arange(0,1024,1/2))+(-0.1+(0.1-(-0.1)))*np.random.rand(1024*2))
        f,t,Sxx=sm.spectrogrampsd(x,2,256,128,'hamming',2048,'db',0,'yes')

        x=sp.signal.chirp(np.arange(0,4,0.001),50,2,150,'quadratic')
        f,t,Sxx=sm.spectrogrampsd(x,1000,128,64,'hamming',128,'db',0,'yes')

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
    from scipy import interpolate
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
    
    x=type2numpy(x)

    #--------------------------------------------------------------------------

    #Detrending the input data
    x=sp.signal.detrend(x)
    
    N=len(x) #Length of the input data

    #Checking size of the window
    Nseg=SegmentSize #Segment size
    if (Nseg>N): Nseg=N
    if (Nseg<=0): Nseg=N

    #Checking size of the overlap window
    Nlap=OverlapSize #Overlap size
    if (Nlap>=Nseg): Nlap=Nseg-1
    if (Nlap<0): Nlap=Nseg-1

    #--------------------------------------------------------------------------
    #Creating window function in time domain

    if WindowName!='none':
        
        #Creating window function in time domain
        Nwseg=Nseg #Number of elements in a window function
        nwinseg=np.arange(0,Nwseg,1) #Counter from 0 to Nw-1
        
        #Calculating window function
        WnRec=np.ones(Nwseg) #Rectangular window function (moving average filter)
        WnTri=1-abs((nwinseg-((Nwseg-1)/2))/((Nwseg-1)/2)) #Triangular or Bartlett window function
        WnWelch=1-((nwinseg-((Nwseg-1)/2))/((Nwseg-1)/2))**2 #Welch window function
        WnHann=0.5-0.5*np.cos(2*np.pi*nwinseg/(Nwseg-1)) #Hanning window function
        WnHamm=0.54-0.46*np.cos(2*np.pi*nwinseg/(Nwseg-1)) #Hamming window function
        WnGauss=np.exp(-0.5*((nwinseg-((Nwseg-1)/2))/(0.4*(Nwseg-1)/2))**2) #Gaussian window function
        WnBlack=0.42659-0.49656*np.cos(2*np.pi*nwinseg/(Nwseg-1))+0.076849*np.cos(4*np.pi*nwinseg/(Nwseg-1)) #Blackman window function
        WnNutt=0.355768-0.487396*np.cos(2*np.pi*nwinseg/(Nwseg-1))+0.144232*np.cos(4*np.pi*nwinseg/(Nwseg-1))-0.012604*np.cos(6*np.pi*nwinseg/(Nwseg-1)) #Nuttall window function
        WnBlHarr=0.35875-0.48829*np.cos(2*np.pi*nwinseg/(Nwseg-1))+0.14128*np.cos(4*np.pi*nwinseg/(Nwseg-1))-0.01168*np.cos(6*np.pi*nwinseg/(Nwseg-1)) #Blackman-Harris window function
        
        #Assigning window function
        if WindowName=='rectangular':
            Wn=WnRec
        elif WindowName=='triangular':
            Wn=WnTri
        elif WindowName=='welch':
            Wn=WnWelch
        elif WindowName=='hanning':
            Wn=WnHann
        elif WindowName=='hamming':
            Wn=WnHamm
        elif WindowName=='gaussian':
            Wn=WnGauss
        elif WindowName=='blackman':
            Wn=WnBlack
        elif WindowName=='nuttall':
            Wn=WnNutt
        elif WindowName=='blackmanharris':
            Wn=WnBlHarr

        Wnseg=Wn #Window Function
        
        WnNormseg=Wnseg/np.sum(Wnseg) #Normalizing a window function in time domain
        WnSqrSumseg=np.sum((Wnseg)**2) #Sum of the window function in time domain
        U=WnSqrSumseg/Nwseg #Energy loss from multiplying data by window function in time domain is compansated by 1/U
        
    #--------------------------------------------------------------------------
    #Calculating power spectral density using Welch method

    #Generating frequency vector
    dfsegWithZeroPadding=fs/N #Frequency difference between consecutive samples, df=fs/N
    fsegWithZeroPadding=np.arange(0,fs/2+dfsegWithZeroPadding,dfsegWithZeroPadding) #Frequency vector with (N/2+1) elements from 0 Hz to fNy=fs/2 Hz, equally spaced at df

    dt=1/fs #Time difference between consecutive samples, dt=1/fs
    duration=N/fs #total time of time series
    tTotal=np.arange(0,duration,dt) #Time from 0 to T-dt, equally spaced at dt
    
    #Calculating total number segments
    K=0 #Total number of segments
    i=0 #Counter
    StartIndx=0 #Initating the first index of the fisrt data segment
    EndIndx=Nseg-1 #Initating last index of the fisrt data segment
    while ((EndIndx-Nlap+1)+(Nseg-1)<=N-1 or i==0):
        if i==0:
            StartIndx=i*Nseg #Segment start index
            EndIndx=(i+1)*Nseg-1 #Segment end index
        else:
            StartIndx=(EndIndx-Nlap)+1 #Segment start index
            EndIndx=StartIndx+Nseg-1 #Segment end index
        i=i+1 #Counter
        K=K+1 #Total number of segments
    
    #Calculating power spectral density
    Sxx_Segment=np.zeros((int(N/2+1),K)) #Initilizing Sxx_Segment
    FFT_Segment=np.zeros((int(N/2+1),K),dtype=complex) #Initilizing FFT_Segment
    t_Segment=np.zeros((int(Nseg),K)) #Initilizing t_Segment
    t=np.zeros(K) #Initilizing t
    for i in range(0,K,1):
        
        if i==0:
            StartIndx=i*Nseg #Segment start index
            EndIndx=(i+1)*Nseg-1 #Segment end index
        else:
            StartIndx=(EndIndx-Nlap)+1 #Segment start index
            EndIndx=StartIndx+Nseg-1 #Segment end index
            
        xseg=x[StartIndx:EndIndx+1] #Time seris segment with length of Nseg data points
        
        #Appling window function on each segment
        if WindowName!='none':
            xseg=Wnseg*xseg #Multiplying time data by window function
        
        #Padding zero to the end of x to make its length equal to N
        xsegWithZeroPadding=np.zeros(N)
        xsegWithZeroPadding[0:Nseg]=xseg
        
        #Fast Fourier Transform (FFT) of the water surface elevation from f=0 Hz to f=fs Hz
        xFFT=np.fft.fft(xsegWithZeroPadding)
        
        #Fast Fourier Transform (FFT) of the water surface elevation from f=0 Hz to fNy=fs/2 Hz
        xFFT=xFFT[0:len(fsegWithZeroPadding)] #(N/2+1)elements if N is even and (N+1)/2 elements if N is odd
        
        #Half of the two-sided power spectral density (psd) from f=0 Hz to fNy=fs/2 Hz
        Sxx_2Sided=(1/(Nseg*fs))*(np.abs(xFFT))**2 #Calculating psd using fs in (m^2/Hz)
        # Sxx_2Sided=(dt**2/(N*dt))*(np.abs(xFFT))**2 #Calculating psd using dt
        
        #One-sided power spectral density (psd) from f=0 Hz to fNy=fs/2 Hz
        Sxx_1Sided=Sxx_2Sided
        Sxx_1Sided[1:-1]=2*Sxx_1Sided[1:-1] #one-side-spectrum=2*two-sided-spectrum in (m^2/Hz)
        
        #Compensating energy loss by multiplying by 1/U
        if WindowName!='none':
            Sxx_Segment[:,i]=Sxx_1Sided/U #Multiplying by 1/U to compensate energy loss from multiplying time data by window function
        else:
            Sxx_Segment[:,i]=Sxx_1Sided             
        
        FFT_Segment[:,i]=xFFT #FFT of each segment

        t_Segment[:,i]=tTotal[StartIndx:EndIndx+1] #Time segment with length of Nseg data points
        t[i]=np.mean(t_Segment[:,i]) #Time at midpoint of each section
        
    SxxAvg=np.mean(Sxx_Segment,1) #Averaging value for each frequency to get mean power spectral density

    #--------------------------------------------------------------------------
    # Converting length of Syy from SegmentSize to nfft elements

    #Checking size of the nfft
    if (nfft>N): nfft=N
    if (nfft<=0): nfft=N
    
    if nfft==N:
        f=fsegWithZeroPadding.copy()
        Sxx=Sxx_Segment.copy()
        FFT=FFT_Segment.copy()
    else:
        #Generating frequency vector with nfft elements
        dfnfft=fs/nfft #Frequency difference between consecutive samples, df=fs/N
        fnfft=np.arange(0,fs/2+dfnfft,dfnfft) #Frequency vector with (N/2+1) elements from 0 Hz to fNy=fs/2 Hz, equally spaced at df
        f=fnfft.copy()

        Sxx=np.zeros((len(f),K)) #Initilizing Sxx with size (int(nfft/2+1),K)
        FFT=np.zeros((len(f),K),dtype=complex) #Initilizing FFT with size (int(nfft/2+1),K)
        for i in range(0,K,1):
            Sxx[:,i]=np.interp(fnfft,fsegWithZeroPadding,Sxx_Segment[:,i]) #Linear interpolation of (fwin,Syywin) to (f,Syy)
            #FFT[:,i]=np.interp(fnfft,fsegWithZeroPadding,FFT_Segment[:,i]) #Linear interpolation of (fwin,FFTwin) to (f,Syy)
            fun=sp.interpolate.interp1d(fsegWithZeroPadding,FFT_Segment[:,i],fill_value='extrapolate') #Linear and complex interpolation
            FFT[:,i]=fun(fnfft) #Linear interpolation of (fwin,FFTwin) to (f,Syy)

    #--------------------------------------------------------------------------
    #Set output type

    if outtype=='psd':
        Sxx=Sxx #Power spectral density (m^2/Hz)
    elif outtype=='db':
        #Sxx=10.*log10(abs(Sxx)./max(abs(Sxx(:)))); #Power spectral density, Converting to decibel, (dB/Hz)
        Sxx=20*np.log10(np.abs(FFT)/np.max(np.abs(FFT))) #Power spectral density, Converting to decibel, (dB/Hz)

    #--------------------------------------------------------------------------
    #Smoothing a power spectral density using convolution

    if OutputSmoothSize>2: #Window size at least should be 3
        
        #Checking size of the window
        if (OutputSmoothSize>N): OutputSmoothSize=N
        
        #Creating window function in time domain
        Nw=OutputSmoothSize #Number of elements in a window function
        nwin=np.arange(0,Nw,1) #Counter from 0 to Nw-1
        
        #Calculating window function
        WnRec=np.ones(Nw) #Rectangular window function (moving average filter)
        WnTri=1-abs((nwin-((Nw-1)/2))/((Nw-1)/2)) #Triangular or Bartlett window function
        WnWelch=1-((nwin-((Nw-1)/2))/((Nw-1)/2))**2 #Welch window function
        WnHann=0.5-0.5*np.cos(2*np.pi*nwin/(Nw-1)) #Hanning window function
        WnHamm=0.54-0.46*np.cos(2*np.pi*nwin/(Nw-1)) #Hamming window function
        WnGauss=np.exp(-0.5*((nwin-((Nw-1)/2))/(0.4*(Nw-1)/2))**2) #Gaussian window function
        WnBlack=0.42659-0.49656*np.cos(2*np.pi*nwin/(Nw-1))+0.076849*np.cos(4*np.pi*nwin/(Nw-1)) #Blackman window function
        WnNutt=0.355768-0.487396*np.cos(2*np.pi*nwin/(Nw-1))+0.144232*np.cos(4*np.pi*nwin/(Nw-1))-0.012604*np.cos(6*np.pi*nwin/(Nw-1)) #Nuttall window function
        WnBlHarr=0.35875-0.48829*np.cos(2*np.pi*nwin/(Nw-1))+0.14128*np.cos(4*np.pi*nwin/(Nw-1))-0.01168*np.cos(6*np.pi*nwin/(Nw-1)) #Blackman-Harris window function
        
        #Assigning window function
        if WindowName=='rectangular':
            Wn=WnRec
        elif WindowName=='triangular':
            Wn=WnTri
        elif WindowName=='welch':
            Wn=WnWelch
        elif WindowName=='hanning':
            Wn=WnHann
        elif WindowName=='hamming':
            Wn=WnHamm
        elif WindowName=='gaussian':
            Wn=WnGauss
        elif WindowName=='blackman':
            Wn=WnBlack
        elif WindowName=='nuttall':
            Wn=WnNutt
        elif WindowName=='blackmanharris':
            Wn=WnBlHarr
        elif WindowName=='none':
            Wn=WnHamm

        WnNorm=Wn/np.sum(Wn) #Normalizing a window function in time domain
        
        for i in range(0,K,1):
            Sxx[:,i]=np.convolve(Sxx[:,i],WnNorm,'same') #Smoothing a power spectral density
        
    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
        
        if len(Sxx[1,:])>1:
            #X,Y=np.meshgrid(t,f)
            #plt.contourf(X,Y,Sxx)
            
            #plt.pcolormesh(t,f,Sxx)

            plt.imshow(Sxx,extent=[t[0],t[-1],f[0],f[-1]],aspect='auto',origin='lower')
            
            if outtype=='psd':
                plt.title('Power Spectral Density (m^2/Hz)')
            elif outtype=='db':
                plt.title('Power Spectral Density (dB/Hz)')

            plt.xlabel('Time (s)')
            plt.ylabel('Frequency (Hz)')
            plt.colorbar()

    #--------------------------------------------------------------------------
    #Outputs
    return f, t, Sxx

    #--------------------------------------------------------------------------
