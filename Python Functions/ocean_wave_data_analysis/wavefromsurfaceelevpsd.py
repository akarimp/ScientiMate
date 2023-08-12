def wavefromsurfaceelevpsd(Eta, fs, fcL=0, fcH=None, nfft=None, SegmentSize=256, OverlapSize=128, dispout='no'):
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

    scientimate.wavefromsurfaceelevpsd
    ==================================

    .. code:: python

        Hm0, fp, Tp, Tm01, Tm02, f, Syy, m0 = scientimate.wavefromsurfaceelevpsd(Eta, fs, fcL=0, fcH=None, nfft=None, SegmentSize=256, OverlapSize=128, dispout='no')

    Description
    -----------

    Calculate wave properties from water surface elevation power spectral density

    Inputs
    ------

    Eta
        Water surface elevation time series data in (m)
    fs
        Sampling frequency that data collected at in (Hz)
    fcL=0
        Low cut-off frequency, between 0*fs to 0.5*fs (Hz)
    fcH=fs/2
        High cut-off frequency, between 0*fs to 0.5*fs (Hz)
    nfft=length(Eta)
        Total number of points between 0 and fs that spectrum reports at is (nfft+1)
    SegmentSize=256
        Segment size, data are divided into the segments each has a total element equal to SegmentSize
    OverlapSize=128
        Number of data points that are overlaped with data in previous segments 
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
        Hm0,fp,Tp,Tm01,Tm02,f,Syy,m0=sm.wavefromsurfaceelevpsd(Eta,fs,0,fs/2,N,256,128,'yes')

    References
    ----------

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
    
    Eta=type2numpy(Eta)

    #--------------------------------------------------------------------------
    #Assign default values

    if fcH is None: fcH=fs/2
    if nfft is None: nfft=len(Eta)

    #--------------------------------------------------------------------------
    #Preventing error in scipy signal by making sure it is a single column array
    #Eta=Eta[:,0]

    #--------------------------------------------------------------------------

    #Deterending input data
    EtaDetrended=sp.signal.detrend(Eta)
    if (fcL<0): fcL=0
    if (fcH>fs/2): fcH=int(fs/2)
    #nfft=int(2**(np.ceil(np.log2(len(EtaDetrended)))))

    #--------------------------------------------------------------------------
    #Calculating power spectral density

    #Water surface elevation power spectral density and frequency from Welch method
    f,Syy=sp.signal.welch(EtaDetrended,fs,window='hamming',nperseg=SegmentSize,noverlap=OverlapSize,nfft=nfft) 
    
    df=f[1]-f[0] #Frequency interval 

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
