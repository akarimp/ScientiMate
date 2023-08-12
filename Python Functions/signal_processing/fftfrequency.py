def fftfrequency(N, fs=2):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2020-05-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    scientimate.fftfrequency
    ========================
    
    .. code:: python
    
        freq = scientimate.fftfrequency(N, fs)
    
    Description
    -----------
    
    Return frequencies for Fast Fourier Transform
    
    Inputs
    ------
    
    N
        Total number of data points
    fs=2
        Sampling frequency (Hz)
    
    Outputs
    -------
    
    freq
        Frequency in (Hz)
    
    Examples
    --------
    
    .. code:: python
    
        import scientimate as sm

        freq = sm.fftfrequency(64, 4)
        freq = sm.fftfrequency(33, 2)
    
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
    
    #x=type2numpy(x)

    #--------------------------------------------------------------------------
    #Calculate frequencies
    #https://www.mathworks.com/matlabcentral/answers/141271-what-are-the-frequencies-when-n-in-fft-x-n-is-odd
    
    dt=1/fs #Time interval
    df=fs/N #Frequency interval
    
    #Frequencies
    if (N%2)==0:
        freq = np.arange(0,N,1)*df #All frequencies
        freq[int(N/2):] = freq[int(N/2):]-fs #Map negative frequencies
    else:
        freq = np.arange(0,N,1)*df #All frequencies
        freq[int((N+1)/2):] = freq[int((N+1)/2):]-fs #Map negative frequencies
    
    #--------------------------------------------------------------------------
    #Outputs
    return freq

    #--------------------------------------------------------------------------
    
#freq = fftfrequency(64, 4)
freq = fftfrequency(33, 2)
