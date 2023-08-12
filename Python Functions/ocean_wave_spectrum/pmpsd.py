def pmpsd(U195=10, fp=0.33, fs=2, N=256, CalSpectralSP='yes', dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2017-08-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.pmpsd
    =================

    .. code:: python

        f, Syy, Hm0, fp, Tp, Tm01, Tm02 = scientimate.pmpsd(U195=10, fp=0.33, fs=2, N=256, CalSpectralSP='yes', dispout='no')

    Description
    -----------

    Calculate Pierson-Moskowitz spectrum (power spectral density), (Pierson and Moskowitz 1964) 

    Inputs
    ------

    U195=10
        Wind velocity at 19.5 meter above surface level in (m/s)
    fp=0.33
        | Peak wave frequency (fp=1/Tp) in (Hz)
        | If CalSpectralSP='yes'; then fp is calculated from U195
    fs=8
        Sampling frequency that data collected at in (Hz)
    N=256
        Total number of points between 0 and fs that spectrum reports at is (N+1)
    CalSpectralSP='yes'
        Define to calculate spectral shape parameters or not ('yes': calculate, 'no': use given parameters by user)
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    f
        Frequency (Hz)
    Syy
        Wave Energy Power Spectrum (m^2/Hz)
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

    Examples
    --------

    .. code:: python

        import scientimate as sm
        f,Syy,Hm0,fp,Tp,Tm01,Tm02=sm.pmpsd(10,0.33,2,256,'yes','yes')

    References
    ----------

    Pierson, W. J., & Moskowitz, L. (1964). 
    A proposed spectral form for fully developed wind seas based on the similarity theory of SA Kitaigorodskii. 
    Journal of geophysical research, 69(24), 5181-5190.

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

    U195=type2numpy(U195)
    fp=type2numpy(fp)

    #--------------------------------------------------------------------------
    #Calculating Pierson-Moskowitz Spectrum
    if (fs<0): fs=-fs
    df=fs/N #Frequency difference between consecutive samples, df=fs/N
    
    f=np.arange(0,fs/2+df,df) #Frequency vector from 0 Hz to fNy=fs/2 Hz, equally spaced at df

    #Calculating spectral shape parameters
    if CalSpectralSP=='yes':
        fp=0.87*9.81/(2*np.pi*U195) #Peak frequency

    alpha=8.1e-3 #Phillip’s constant

    f[0]=f[1]/2 #Assign temporarily non-zero value to fisrt element of f to prevent division by zero
    Syy=alpha*9.81**2/((2*np.pi)**4*(f**5))*np.exp(-1.25*(fp/f)**4) #Calculating Pierson-Moskowitz Spectrum 

    #Assign zero value to fisrt element of f and Syy
    f[0]=0
    Syy[0]=0

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
        plt.plot(f[f!=0],Syy[f!=0])
    
        plt.title('Pierson-Moskowitz Power Spectral Density')
        plt.xlabel('Frequency(Hz)')
        plt.ylabel('Spectral Density(m^2/Hz)')

    #--------------------------------------------------------------------------
    #Outputs
    return f, Syy, Hm0, fp, Tp, Tm01, Tm02

    #--------------------------------------------------------------------------
