def filtertimeseries(x, fs, fcL=0, fcH=None, dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2020-09-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    scientimate.filtertimeseries
    ============================
    
    .. code:: python
    
        xFiltered, t = scientimate.filtertimeseries(x, fs, fcL, fcH, dispout)
    
    Description
    -----------
    
    | Filter time-series to retain signals with frequencies of fcL <= f <= fcH
    | It assumes time-series is stationary
    
    Inputs
    ------
    
    x
        | Time-series
    fs
        | Sampling frequency that data collected at in (Hz)
        | fs=1/dt where dt is time interval between two successive points
    fcL=0
        | Low cut-off frequency, between 0*fs to 0.5*fs (Hz)
        | Signals with frequencies f<fcL are removed
        | Must be between 0 and fs/2
    fcH=fs/2
        | High cut-off frequency, between 0*fs to 0.5*fs (Hz)
        | Signals with frequencies f>fcH are removed
        | Must be between 0 and fs/2
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)
    
    Outputs
    -------
    
    xFiltered
        Filtered time-series
    t
        Time
    
    Examples
    --------
    
    .. code:: python
    
        import scientimate as sm
        import numpy as np

        #Generate time-series from 3 waves with frequencies of 0.5, 2, 4 Hz
        fs=32 #Sampling frequency
        d=20 #Duration
        f1=0.5 #1st wave frequency
        f2=2 #2nd wave frequency
        f3=4 #3rd wave frequency
        a1=1 #1st wave amplitude
        a2=0.2 #2nd wave amplitude
        a3=0.1 #3rd wave amplitude
        dt=1/fs #Time interval
        t=np.linspace(0,d-dt,fs*d) #Time data points
        x=a1*np.sin(2*np.pi*f1*t)+a2*np.sin(2*np.pi*f2*t)+a3*np.sin(2*np.pi*f3*t) #Time-series
    
        #Keep only wave with f1 frequency and remove waves with f2 and f3 frequencies
        fcL=f1-0.2 #Lower cut-off frequency
        fcH=f1+0.2 #Higher cut-off frequency
        xFiltered, time = sm.filtertimeseries(x, fs, fcL, fcH, 'yes')
    
        #Keep waves with f2 and f3 frequencies and remove wave with f1 frequency
        fcL=f2-0.2 #Lower cut-off frequency
        fcH=f3+0.2 #Higher cut-off frequency
        xFiltered, time = sm.filtertimeseries(x, fs, fcL, fcH, 'yes')
    
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
    from numpy import fft
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
    #Assign default values

    if fcH is None: fcH=fs/2

    #--------------------------------------------------------------------------
    
    N=len(x) #Length of the input data
    xDetrended=sp.signal.detrend(x,type='constant') #Deterending input data
    xMean=np.mean(x) #Mean
    if (fcL<0): fcL=0
    if (fcH>fs/2): fcH=int(fs/2)
    dt=1/fs #Time interval
    t=np.arange(0,N,1)*dt #Time data points
    
    #--------------------------------------------------------------------------
    #Filter time-series
    
    #Windowed-sinc filter order (must be even)
    if (N%2)==0:
        nf=N 
    else:
        nf=N+1 #Make number of datapoints even by padding a single zero at the end of FFT
    
    #Band-pass windowed-sinc filter
    b = sp.signal.firwin(nf+1, [fcL,fcH], window='hann', pass_zero='bandpass', fs=fs) #b is impulse response
    f, h = sp.signal.freqz(b, a=1, worN=len(b)-1, whole=True, fs=fs) #h is frequency response
    #plt.plot(f,abs(h)) #Plot frequency response vs frequency

    #FFT of time-series
    x_fft=np.fft.fft(xDetrended,nf)

    #Apply band-pass windowed-sinc filter to data
    x_fft_Filtered = x_fft*np.abs(h)

    #Inverst FFT
    x_ifft_Filtered=np.real(np.fft.ifft(x_fft_Filtered));

    #Filtered time-series
    if (N%2)==0:
        xFiltered=x_ifft_Filtered.copy()
    else:
        xFiltered=x_ifft_Filtered[0:-1]

    #Add mean back to time-series
    xFiltered=xFiltered+xMean

    #--------------------------------------------------------------------------
    #Displaying results
    
    if dispout=='yes':
        
        plt.plot(t,x,label='Input')
        plt.plot(t,xFiltered,label='Filtered')
        
        plt.title('Filtered Data')
        plt.xlabel('Time (s)')
        plt.ylabel('Data')
        plt.legend()
        
    
    #--------------------------------------------------------------------------
    #Outputs
    return xFiltered, t

    #--------------------------------------------------------------------------
    