def wavepropfrompsd(Syy, f, fcL=0, fcH=None, dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2017-05-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.wavepropfrompsd
    ===========================

    .. code:: python

        Hm0, fp, Tp, Tm01, Tm02, m0, m1, m2 = scientimate.wavepropfrompsd(Syy, f, fcL=0, fcH=None, dispout='no')

    Description
    -----------

    Calculate wave properties from a power spectral density

    Inputs
    ------

    Syy
        Power spectral density (m^2/Hz)
    f
        Frequency (Hz)
    fcL=0
        Low cut-off frequency, between 0*fs to 0.5*fs (Hz)
    fcH=max(f)
        High cut-off frequency, between 0*fs to 0.5*fs (Hz)
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
    m0
        Zero-Moment of the power spectral density (m^2)
    m1
        First-Moment of the power spectral density (m^2)
    m2
        Second-Moment of the power spectral density (m^2)

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        N=2**11
        fs=8
        df=fs/N #Frequency difference 
        f=np.arange(0,fs/2+df,df) #Frequency vector 
        f[0]=f[1]/2 #Assign temporarily non-zero value to fisrt element of f to prevent division by zero
        Syy=0.016*9.81**2/((2*np.pi)**4*(f**5))*np.exp(-1.25*(0.33/f)**4) #Calculating Spectrum 
        f[0]=0
        Syy[0]=0
        Hm0,fp,Tp,Tm01,Tm02,m0,m1,m2=sm.wavepropfrompsd(Syy,f,0,np.max(f),'yes')

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
    
    Syy=type2numpy(Syy)
    f=type2numpy(f)

    #--------------------------------------------------------------------------
    #Assign default values

    if fcH is None: fcH=np.max(f)

    #--------------------------------------------------------------------------

    if (fcL<0): fcL=0
    if (fcH>np.max(f)): fcH=np.max(f)
    df=f[1]-f[0] #Frequency interval 

    #--------------------------------------------------------------------------
    #Cut off spectrum based on fcL and fcH

    if ((fcL>f[0]) and (fcL<np.max(f))):
        #Indx=np.int_((np.nonzero(f<fcL))[0])
        #Syy=np.delete(Syy,Indx)
        #f=np.delete(f,Indx)
        Syy[f<fcL]=0

    if ((fcH>f[0]) and (fcH<np.max(f))):
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
        plt.plot(f,Syy)
    
        plt.title('Power Spectral Density')
        plt.xlabel('Frequency(Hz)')
        plt.ylabel('Spectral Density(m^2/Hz)')

    #--------------------------------------------------------------------------
    #Outputs
    return Hm0, fp, Tp, Tm01, Tm02, m0, m1, m2

    #--------------------------------------------------------------------------
