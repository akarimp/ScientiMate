def smoothsignal(x, WindowSize=33, method='moveavg', dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2020-01-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.smoothsignal
    ========================

    .. code:: python

        xSmoothed = scientimate.smoothsignal(x, WindowSize=33, method='moveavg', dispout='no')

    Description
    -----------

    Smooth input data using a window function

    Inputs
    ------

    x
        Input data
    WindowSize=33
        Window size (number of adjacent elements)that is used for smoothing, should be equal or larger than 3
    method='moveavg'
        | Smoothing method
        | 'moveavg': moving average
        | 'lowpass': Low-Pass filter
        | 'savgol': Savitzky-Golay filter
        | 'butter': Butterworth filter
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    xSmoothed
        Smoothed data

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np
        from numpy import random

        rng = np.random.default_rng()
        x=np.sin(np.linspace(0,6*np.pi,200))+rng.random(200)
        xSmoothed=sm.smoothsignal(x,33,'moveavg','yes')

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

    N=len(x) #Length of the input data
    fs=2 #Sampling frequency
    fc=0.25 #Cut-off frequency (Considered as fc/(fs/2)=0.25)

    #Checking size of the window
    Nw=WindowSize
    if Nw>N: Nw=int(N/8)
    if (Nw%2)==0: Nw=Nw+1 #Window and filter length must be odd
    if Nw<=3: Nw=3 #Window and filter length should be equal or larger than 3


    #--------------------------------------------------------------------------
    #Smoothig data

    if method=='moveavg':
        win = sp.signal.windows.boxcar(Nw)
        win_norm = win/np.sum(win) #Normalize 'boxcar' window function
        xSmoothed = sp.signal.convolve(x, win_norm, mode='same') #Smoothing input using convolution

    elif method=='lowpass':
        impulse_response = sp.signal.firwin(Nw, fc, window='hann', pass_zero='lowpass', fs=fs)
        xSmoothed = sp.signal.convolve(x, impulse_response, mode='same')

    elif method=='savgol':
        poly_order=5 #Polynomial order
        xSmoothed = sp.signal.savgol_filter(x, Nw, poly_order, mode='constant', cval=0)

    elif method=='butter':
        filter_order=4 #Filter order
        b, a = sp.signal.butter(filter_order, fc, btype='lowpass', fs=fs) #Design Butterworth filter
        xSmoothed = sp.signal.filtfilt(b, a, x)


    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
        
        plt.plot(x,label='Input')
        plt.plot(xSmoothed,label='Smoothed')
        
        plt.title('Soothed Data')
        plt.xlabel('Samples')
        plt.ylabel('Data')
        plt.legend(loc='upper right')
        
    #--------------------------------------------------------------------------
    #Outputs
    return xSmoothed

    #--------------------------------------------------------------------------
