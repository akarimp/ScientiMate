def bartlettpsd(x, fs=2, SegmentSize=256, WindowName='hamming', nfft=256, OutputSmoothSize=0, dispout='no'):
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

    scientimate.bartlettpsd
    =======================

    .. code:: python

        f, Sxx = scientimate.bartlettpsd(x, fs=2, SegmentSize=256, WindowName='hamming', nfft=256, OutputSmoothSize=0, dispout='no')

    Description
    -----------

    Calculate power spectral density using Bartlett's method

    Inputs
    ------

    x
        Input data
    fs=2
        Sampling frequency that data collected at in (Hz)
    SegmentSize=256
        Segment size, data are divided into the segments each has a total element equal to SegmentSize
    WindowName='hamming'
        | Window name, define if multiplying input data by window function or not ('none': not multiplying)
        | 'none','rectangular','triangular','welch','hanning','hamming','gaussian','blackman','nuttall','blackmanharris'
    nfft=length(x)
        Total number of points between 0 and fs that spectrum reports at is (nfft+1)
    OutputSmoothSize=0
        | Window size for smoothing calculated spectrum (0, 1 or 2: not smoothing, reports original periodogram)
        | if WindowName='none' and OutputSmoothSize>2, then WindowName='hamming'
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    f
        Frequency in (Hz)
    Sxx
        Power spectral density using Bartlett's method (m^2/Hz)

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np
        import scipy as sp
        from scipy import signal

        x=sp.signal.detrend(0.5*np.cos(2*np.pi*0.2*np.arange(0,1024,1/2))+(-0.1+(0.1-(-0.1)))*np.random.rand(1024*2))
        f,Sxx=sm.bartlettpsd(x,2,256,'hamming',2048,0,'yes')

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
    #Calculating power spectral density using Bartlett method

    K=int(np.fix(N/Nseg)) #Total number of segments

    #Generating frequency vector
    dfsegWithZeroPadding=fs/N #Frequency difference between consecutive samples, df=fs/N
    fsegWithZeroPadding=np.arange(0,fs/2+dfsegWithZeroPadding,dfsegWithZeroPadding) #Frequency vector with (N/2+1) elements from 0 Hz to fNy=fs/2 Hz, equally spaced at df

    Sxx_Segment=np.zeros((int(N/2+1),K)) #Initilizing Sxx_Segment
    for i in range(0,K,1):
        
        StartIndx=i*Nseg #Segment start index
        EndIndx=(i+1)*Nseg-1 #Segment end index
        xseg=x[StartIndx:EndIndx+1] #Time seris segment with length of Nseg data points
        
        #Multiplying time data by window function
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
        
        Sxx_Segment[:,i]=Sxx_1Sided
        
    
    SxxAvg=np.mean(Sxx_Segment,1) #Averaging value for each frequency to get mean power spectral density
    
    #Compensating energy loss by multiplying time data by window function
    if WindowName!='none':
        SxxAvg=SxxAvg/U #Compensating energy loss by multiplying by 1/U

    #--------------------------------------------------------------------------
    # Converting length of Syy from SegmentSize to nfft elements

    #Checking size of the nfft
    if (nfft>N): nfft=N
    if (nfft<=0): nfft=N
    
    if nfft==N:
        f=fsegWithZeroPadding.copy()
        Sxx=SxxAvg.copy()
    else:
        #Generating frequency vector with nfft elements
        dfnfft=fs/nfft #Frequency difference between consecutive samples, df=fs/N
        fnfft=np.arange(0,fs/2+dfnfft,dfnfft) #Frequency vector with (N/2+1) elements from 0 Hz to fNy=fs/2 Hz, equally spaced at df
        
        f=fnfft.copy()
        #Sxx=interp1(fsegWithZeroPadding,SxxAvg,fnfft); #Linear interpolation of (fwin,Syywin) to (f,Syy)
        Sxx=np.interp(fnfft,fsegWithZeroPadding,SxxAvg) #Linear interpolation of (fwin,Syywin) to (f,Syy)

        
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
        
        Sxx=np.convolve(Sxx,WnNorm,'same') #Smoothing a power spectral density

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
        
        plt.loglog(f[f!=0],Sxx[f!=0])
        
        plt.title('Power Spectral Density')
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Spectral Density (m^2/Hz)')
        

    #--------------------------------------------------------------------------
    #Outputs
    return f, Sxx

    #--------------------------------------------------------------------------
