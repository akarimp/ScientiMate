def seaswell1d(f, Syy, fsepCalcMethod='hwang', fu=0.5, fmaxswell=0.2, fpminswell=0, Windvel=10, dispout='no'):
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

    scientimate.seaswell1d
    ======================

    .. code:: python

        fseparation, Hm0, Hm0sea, Hm0swell, Tp, Tpsea, Tpswell, m0, m0sea, m0swell, fp, Tm01, Tm02 = scientimate.seaswell1d(f, Syy, fsepCalcMethod='hwang', fu=0.5, fmaxswell=0.2, fpminswell=0, Windvel=10, dispout='no')

    Description
    -----------

    Partition (separate) wind sea from swell in a power spectral density using an one dimensional method

    Inputs
    ------

    f
        Frequency (Hz)
    Syy
        Power spectral density (m^2/Hz)
    fsepCalcMethod='hwang'
        | Wind sea swell separating calculation method 
        | 'celerity': using deep water wave celerity, 'gilhousen': Gilhousen and Hervey (2001), 
        | 'hwang': Hwang et al. (2012), 'exact': calculate exact value 
    fu=0.5
        An upper bound of a spectrum integration frequency (Hz)
    fmaxswell=0.25
        Maximum frequency that swell can have, It is about 0.2 in Gulf of Mexico
    fpminswell=0.1
        A lower bound of a spectrum (minimum frequency) that is used for Tpswell calculation
    Windvel=10
        Wind velocity (m/s), only required for Gilhousen and Hervey (2001) method
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    fseparation
        Wind sea and Swell Separation Frequency (Hz)
    Hm0
        Zero-Moment Wave Height (m)
    Hm0sea
        Sea Zero-Moment Wave Height (m)
    Hm0swell
        Swell Zero-Moment Wave Height (m)
    Tp
        Peak wave period (second)
    Tpsea
        Peak Sea period (second)
    Tpswell
        Peak Swell Period (second)
    m0
        Zero-Moment of the power spectral density (m^2)
    m0sea
        Zero-Moment of the wind sea spectrum (m^2)
    m0swell
        Zero-Moment of the swell spectrum (m^2)
    fp
        Peak wave frequency (Hz)
    Tm01
        Wave Period from m01 (second), Mean Wave Period
    Tm02
        Wave Period from m02 (second), Mean Zero Crossing Period

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        N=2**11 #Total number of points
        fs=2 #Sampling frequency
        df=fs/N #Frequency difference 
        f=np.arange(0,fs/2+df,df) #Frequency vector 
        f[0]=f[1]/2 #Assign temporarily non-zero value to fisrt element of f to prevent division by zero
        SyySea=0.016*9.81**2/((2*np.pi)**4*(f**5))*np.exp(-1.25*(0.33/f)**4) #Calculating Spectrum 
        SyySwell=0.016*9.81**2/((2*np.pi)**4*(f**5))*np.exp(-1.25*(0.15/f)**4)*0.005 #Calculating Spectrum 
        Syy=SyySea+SyySwell
        f[0]=0
        Syy[0]=0
        fseparation,Hm0,Hm0sea,Hm0swell,Tp,Tpsea,Tpswell,m0,m0sea,m0swell,fp,Tm01,Tm02=sm.seaswell1d(f,Syy,'exact',0.5,0.3,0,10,'yes')

    References
    ----------

    Gilhousen, D. B., & Hervey, R. (2002). 
    Improved estimates of swell from moored buoys. 
    In Ocean Wave Measurement and Analysis (2001) (pp. 387-393).

    Hwang, P. A., Ocampo-Torres, F. J., & García-Nava, H. (2012). 
    Wind sea and swell separation of 1D wave spectrum by a spectrum integration method. 
    Journal of Atmospheric and Oceanic Technology, 29(1), 116-128.

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
    
    f=type2numpy(f)
    Syy=type2numpy(Syy)

    #--------------------------------------------------------------------------

    df=f[1]-f[0] #Frequency difference between consecutive samples, df=fs/N

    #--------------------------------------------------------------------------
    #Calculating 1D separation frequency of wind sea from swell

    #Calculating 1D separation frequency of wind sea from swell using wave celerity (phase speed)
    if fsepCalcMethod=='celerity':
    
        # fseparation is associated with a frequency where wave celerity becomes equal to wind velocity 
        fseparation=9.81/(2*np.pi*Windvel)  

    #Calculating 1D separation frequency of wind sea from swell using Gilhousen and Hervey (2001) method
    elif fsepCalcMethod=='gilhousen':
        
        fup=f[f<=fu] #Frequency up to fu
        fstar=np.zeros(len(fup))
        m0fstar=np.zeros(len(fup))
        m2fstar=np.zeros(len(fup))
        for i in range(0,len(fup),1):
            fstar[i]=fup[i]
            m0fstar[i]=np.sum(Syy[i:len(fup)+1]*f[i:len(fup)+1]**(0)*df)
            m2fstar[i]=np.sum(Syy[i:len(fup)+1]*f[i:len(fup)+1]**(2)*df)

        
        alfafstar=(8*np.pi*m2fstar)/(9.81*np.sqrt(m0fstar)) #Steepness function
        
        alfafstarmaxIndx=np.nanargmax(alfafstar) #Index of the maximum value
        fx=fstar[alfafstarmaxIndx]
        fseparation=np.max([0.75*fx,0.9*1.25/Windvel])
        
        if ((fseparation>fmaxswell) or (fseparation==0)): fseparation=fmaxswell #fseperation is about 0.2 in Gulf of Mexico
        
    #Calculating 1D separation frequency of wind sea from swell using Hwang et al. (2012) method
    elif fsepCalcMethod=='hwang':
        
        fup=f[f<=fu] #Frequency up to fu
        fstar=np.zeros(len(fup))
        m1fstar=np.zeros(len(fup))
        mminus1fstar=np.zeros(len(fup))
        for i in range(0,len(fup),1):
            fstar[i]=fup[i]
            m1fstar[i]=np.sum(Syy[i:len(fup)+1]*f[i:len(fup)+1]**1*df)
            mminus1fstar[i]=np.sum(Syy[i:len(fup)+1]*f[i:len(fup)+1]**(-1)*df)

        
        alfafstar=(m1fstar)/np.sqrt(mminus1fstar)
        
        alfafstarmaxIndx=np.nanargmax(alfafstar) #Index of the maximum value
        fm=fstar[alfafstarmaxIndx]
        fseparation=24.2084*fm**3-9.2021*fm**2+1.8906*fm-0.04286
        
        if ((fseparation>fmaxswell) or (fseparation==0)): fseparation=fmaxswell #fseperation is about 0.2 in Gulf of Mexico
        
    
    #Calculating the exact separation frequency of wind sea from swell
    elif fsepCalcMethod=='exact':
    
        #Estimating 1D separation frequency of wind sea from swell using Hwang et al. (2012) method
        fup=f[f<=fu] #Frequency up to fu
        fstar=np.zeros(len(fup))
        m1fstar=np.zeros(len(fup))
        mminus1fstar=np.zeros(len(fup))
        for i in range(0,len(fup),1):
            fstar[i]=fup[i]
            m1fstar[i]=np.sum(Syy[i:len(fup)+1]*f[i:len(fup)+1]**1*df)
            mminus1fstar[i]=np.sum(Syy[i:len(fup)+1]*f[i:len(fup)+1]**(-1)*df)

        
        alfafstar=(m1fstar)/np.sqrt(mminus1fstar)
        
        alfafstarmaxIndx=np.nanargmax(alfafstar) #Index of the maximum value
        fm=fstar[alfafstarmaxIndx]
        fseparation=24.2084*fm**3-9.2021*fm**2+1.8906*fm-0.04286
        
        if ((fseparation>fmaxswell) or (fseparation==0)): fseparation=fmaxswell #fseperation is about 0.2 in Gulf of Mexico
    
        #Calculating the exact separation frequency of wind sea from swell
        fsepIndx=len((np.nonzero(f<=fseparation))[0])-1 #location of fseperation
        fpsea1=(np.sum(Syy[fsepIndx:]**5*f[fsepIndx:]**1*df))/(np.sum(Syy[fsepIndx:]**5*f[fsepIndx:]**0*df)) #Wind sea peak frequency based on fseparation from previous step
        fpswell1=(np.sum(Syy[0:fsepIndx+1]**5*f[0:fsepIndx+1]**1*df))/(np.sum(Syy[0:fsepIndx+1]**5*f[0:fsepIndx+1]**0*df)) #Swell peak frequency based on fseparation from previous step
        Indx1=np.argmin(Syy[(f>fpswell1) & (f<fpsea1)]) #Index of the minimum value
        Indx2=len((np.nonzero(f<=fpswell1))[0])-1 
        fseparation=f[Indx1+Indx2+1]
        if ((fseparation>fmaxswell) or (fseparation==0)): fseparation=fmaxswell #fseperation is about 0.2 in Gulf of Mexico
        #if len(fseparation)==0: fseparation=fmaxswell #fseperation is about 0.2 in Gulf of Mexico
    

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
    #Calculating wind sea and swell wave properties

    fsepIndx=len((np.nonzero(f<=fseparation))[0])-1 #Index of fseperation
    
    m0sea=np.sum(Syy[fsepIndx+1:]*f[fsepIndx+1:]**0*df) #Zero-Moment of wind sea 
    m0swell=np.sum(Syy[0:fsepIndx+1]*f[0:fsepIndx+1]**0*df) #Zero-Moment of swell
    
    Hm0sea=4*np.sqrt(m0sea) #Wind sea zero-Moment wave height
    Hm0swell=4*np.sqrt(m0swell) #Swell zero-Moment wave height
    
    SyymaxseaIndx=np.argmax(Syy[fsepIndx:])
    Tpsea=1/(f[fsepIndx+SyymaxseaIndx]) #Wind sea peak period, fsepIndx+SyymaxseaIndx-1 is the location of sea peak
    
    fpminswellIndx=len((np.nonzero(f<=fpminswell))[0])-1 #Index of fpminswell
    if fpminswellIndx>=fsepIndx: fpminswellIndx=1
    SyymaxswellIndx=np.argmax(Syy[fpminswellIndx:fsepIndx+1])
    Tpswell=1/(f[fpminswellIndx+SyymaxswellIndx]) #Swell peak period, fpminswellIndx+SyymaxswellIndx-1 is the location of swell peak
    
    #Calculating peak frequency from weighted integral (Young, 1995)
    fpsea=(np.sum(Syy[fsepIndx:]**5.*f[fsepIndx:]**1*df))/(np.sum(Syy[fsepIndx:]**5*f[fsepIndx:]**0*df)) #Wind sea peak frequency
    fpswell=(np.sum(Syy[0:fsepIndx+1]**5*f[0:fsepIndx+1]**1*df))/(np.sum(Syy[0:fsepIndx+1]**5*f[0:fsepIndx+1]**0*df))#Swell peak frequency
    
    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':

        val=[fseparation,Hm0,Hm0sea,Hm0swell,Tp,Tpsea,Tpswell,m0,m0sea,m0swell,fp,Tm01,Tm02]
        name=['fseparation','Hm0','Hm0sea','Hm0swell','Tp','Tpsea','Tpswell','m0','m0sea','m0swell','fp','Tm01','Tm02']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            print('{0:10}= {1:0.10f}'.format(name[i],val[i])) 

        #Plotting
        plt.plot(f,Syy,label='Spectrum')
        #plt.hold='True'
        plt.plot([fseparation,fseparation],[np.min(Syy),np.max(Syy)],label='Separation Frequency')
    
        plt.title('Power Spectral Density')
        plt.xlabel('Frequency(Hz)')
        plt.ylabel('Spectral Density(m^2/Hz)')
        plt.legend()

    #--------------------------------------------------------------------------
    #Outputs
    return fseparation, Hm0, Hm0sea, Hm0swell, Tp, Tpsea, Tpswell, m0, m0sea, m0swell, fp, Tm01, Tm02

    #--------------------------------------------------------------------------
