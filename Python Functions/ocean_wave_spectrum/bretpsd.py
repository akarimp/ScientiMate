def bretpsd(Hs=1, fp=0.33, fs=2, N=256, dispout='no'):
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

    scientimate.bretpsd
    ===================

    .. code:: python

        f, Syy, Hm0, fp, Tp, Tm01, Tm02 = scientimate.bretpsd(Hs=1, fp=0.33, fs=2, N=256, dispout='no')

    Description
    -----------

    Calculate Bretschneider spectrum (power spectral density), (Bretschneider, 1959), ITTC spectrum 

    Inputs
    ------

    Hs=1
        Significant wave height (m)
    fp=0.33
        Peak wave frequency in (Hz)
    fs=2
        Sampling frequency that data collected at in (Hz)
    N=256
        Total number of points between 0 and fs that spectrum reports at is (N+1)
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
        f,Syy,Hm0,fp,Tp,Tm01,Tm02=sm.bretpsd(1,0.33,2,256,'yes')

    References
    ----------

    Bretschneider, C. L. (1959). 
    Wave variability and wave spectra for wind-generated gravity waves
    (No. TM-118). CORPS OF ENGINEERS WASHINGTON DC BEACH EROSION BOARD.

    Stansberg, C. T., Contento, G., Hong, S. W., Irani, M., Ishida, S., & Mercier, R. (2002). 
    The specialist committee on waves final report and recommendations to the 23rd ITTC. 
    In Proceedings of the 23rd ITTC (Vol. 2, pp. 505-551).

    Zwolan, P., & Czaplewski, K. (2012). 
    Sea waves models used in maritime simulators. 
    Zeszyty Naukowe/Akademia Morska w Szczecinie, (32 (104) z. 2), 186-190.

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

    Hs=type2numpy(Hs)
    fp=type2numpy(fp)

    #--------------------------------------------------------------------------
    #Calculating Bretschneider Spectrum (Bretschneider, 1959), ITTC Spectrum
    if (fs<0): fs=-fs
    df=fs/N #Frequency difference between consecutive samples, df=fs/N
    
    f=np.arange(0,fs/2+df,df) #Frequency vector from 0 Hz to fNy=fs/2 Hz, equally spaced at df

    f[0]=f[1]/2 #Assign temporarily non-zero value to fisrt element of f to prevent division by zero
    # Tz=3.0                         #Mean wave period in (s)
    #Syy=0.0796*(Hs**2*Tz)/(Tz*f)**5*np.exp(-0.318*(1/(Tz*f))**4) #Calculating Bretschneider Spectrum (Bretschneider, 1959)
    Syy=5/16*(Hs**2*fp**4)/(f)**5*np.exp(-5/4*(fp/f)**4) #Calculating Bretschneider Spectrum (Bretschneider, 1959)

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
    
        plt.title('Bretschneider Power Spectral Density')
        plt.xlabel('Frequency(Hz)')
        plt.ylabel('Spectral Density(m^2/Hz)')

    #--------------------------------------------------------------------------
    #Outputs
    return f, Syy, Hm0, fp, Tp, Tm01, Tm02

    #--------------------------------------------------------------------------
