def diagnostictail(fIn, SyyIn, ftail, tailtype='jonswap', tailpower=-5, h=0, transfCalcMethod='approx', kCalcMethod='beji', dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2017-03-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.diagnostictail
    ==========================

    .. code:: python

        fOut, SyyOut, Hm0, fp, Tp, Tm01, Tm02 = scientimate.diagnostictail(fIn, SyyIn, ftail, tailtype='jonswap', tailpower=-5, h=0, transfCalcMethod='approx', kCalcMethod='beji', dispout='no')

    Description
    -----------

    Replace a spectrum tail with JONSWAP (Hasselmann et al., 1973) or TMA Spectrum (Bouws et al., 1985)

    Inputs
    ------

    fIn
        Frequency (Hz), Input
    SyyIn
        Power spectral density (m^2/Hz), Input
    ftail
        Frequency that diagnostic tail apply after that (typically: ftail=2.5fm where fm=1/Tm01)
    tailtype='jonswap'
        | Define type of the diagnostic tail to be applied 
        | 'jonswap': JONSWAP Spectrum tail, 'tma': TMA Spectrum tail
    tailpower=-5
        | Tail power that diagnostic tail apply based on that 
        | tailpower=-3 for shallow water, tailpower=-5 for deep water
    h=0
        Mean water depth in (m)
    transfCalcMethod='approx'
        | Transformation function from JONSWAP into TMA calculation method 
        | 'approx': approximated method, 'tucker': Tucker (1994), 'kitaigordskii': Kitaigordskii et al. (1975) 
    kCalcMethod='beji'
        | Wave number calculation method 
        | 'hunt': Hunt (1979), 'beji': Beji (2013), 'vatankhah': Vatankhah and Aghashariatmadari (2013) 
        | 'goda': Goda (2010), 'exact': calculate exact value 
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    fOut
        Frequency (Hz), Output
    SyyOut
        Power spectral density (m^2/Hz) with diagnostic tail, Output
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
        import numpy as np

        N=2**11 #Total number of points
        fs=8 #Sampling frequency
        df=fs/N #Frequency difference 
        f=np.arange(0,fs/2+df,df) #Frequency vector 
        f[0]=f[1]/2 #Assign temporarily non-zero value to fisrt element of f to prevent division by zero
        Syy=0.016*9.81**2/((2*np.pi)**4*(f**5))*np.exp(-1.25*(0.33/f)**4) #Calculating Spectrum 
        f[0]=0
        Syy[0]=0
        
        fOut,SyyOut,Hm0,fp,Tp,Tm01,Tm02=sm.diagnostictail(f,Syy,0.5,'jonswap',-5,5,'approx','beji','yes')
        
        fOut,SyyOut,Hm0,fp,Tp,Tm01,Tm02=sm.diagnostictail(f,Syy,0.5,'tma',-3,5,'approx','beji','yes')

    References
    ----------

    Beji, S. (2013). 
    Improved explicit approximation of linear dispersion relationship for gravity waves. 
    Coastal Engineering, 73, 11-12.

    Bouws, E.; Günther, H.; Rosenthal, W., and Vincent, C.L., (1985). 
    Similarity of the wind wave spectrum in finite depth water: 1. Spectral form. 
    Journal of Geophysical Research: Oceans, 90(C1), 975-986.

    Goda, Y. (2010). 
    Random seas and design of maritime structures. 
    World scientific.

    Hasselmann, K.; Barnett, T. P.; Bouws, E.; Carlson, H.; Cartwright, D. E.; Enke, K.; Ewing, J. A.; 
    Gienapp, H.; Hasselmann, D. E.; Kruseman, P.; Meerbrug, A.; Muller, P.; Olbers, D. J.; Richter, K.; 
    Sell, W., and Walden, H., (1973). 
    Measurements of wind-wave growth and swell decay during the Joint North Sea Wave Project (JONSWAP). 
    Deutsche Hydrographische Zeitschrift A80(12), 95p.

    Hunt, J. N. (1979). 
    Direct solution of wave dispersion equation. 
    Journal of the Waterway Port Coastal and Ocean Division, 105(4), 457-459.

    Vatankhah, A. R., & Aghashariatmadari, Z. (2013). 
    Improved explicit approximation of linear dispersion relationship for gravity waves: A discussion. 
    Coastal engineering, 78, 21-22.

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
    from scipy import optimize
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
    
    fIn=type2numpy(fIn)
    SyyIn=type2numpy(SyyIn)

    #--------------------------------------------------------------------------

    Syy=SyyIn.copy() #Preserve input Spectrum
    f=fIn.copy() #Preserve input Spectrum
    df=f[1]-f[0] #Frequency difference between consecutive samples, df=fs/N

    #--------------------------------------------------------------------------
    #Calculating PHI for TMA Spectrum

    if tailtype=='tma':
    
        omega=2*np.pi*f*np.sqrt(h/9.81)
    
        #Transformation function from JONSWAP into TMA, approximated method
        if transfCalcMethod=='approx':
            PHI=np.ones(len(omega))
            PHI[omega<=1]=omega[omega<=1]**2/2
            PHI[np.logical_and(omega>1,omega<2)]=1-0.5*(2-omega[np.logical_and(omega>1,omega<2)])**2
            PHI[omega>=2]=1
        
        #Transformation function from JONSWAP into TMA, exact method (Tucker, 1994)
        elif transfCalcMethod=='tucker':
        
            #Calculating wave number (k)
        
            w=2*np.pi*f #Wave angular frequency, w=2*pi/T=2*pi*f
            
            #Deep water
            k0=(w**2)/9.81 #Deep water wave number
            k0h=k0*h
            
            #Estimation of wave number (k) from Hunt (1979)
            if kCalcMethod=='hunt':
                # c2gh=(k0h+(1.0+0.6522*k0h+0.4622*k0h**2+0.0864*k0h**4+0.0675*k0h**5)**(-1))**(-1) #Calculating wave number from Hunt (1979)
                # k=w/(np.sqrt(c2gh*9.81*h))
                d1=0.6666666667; d2=0.3555555556; d3=0.1608465608; d4=0.0632098765; d5=0.0217540484; d6=0.0065407983;
                kh=np.sqrt(k0h**2+(k0h)/(1+d1*k0h+d2*k0h**2+d3*k0h**3+d4*k0h**4+d5*k0h**5+d6*k0h**6)) #Calculating wave number from Hunt (1979)
                k=kh/h #Calculating wave number from Hunt (1979)
            
            
            #Estimation of wave number (k) from Beji (2013)
            elif kCalcMethod=='beji':
                kh=k0h*(1+k0h**1.09*np.exp(-(1.55+1.3*k0h+0.216*k0h**2)))/np.sqrt(np.tanh(k0h)) #Calculating wave number from Beji (2013)
                k=kh/h #Calculating wave number from Beji (2013)
                k[w==0]=0
            
            
            #Estimation of wave number (k) from Vatankhah and Aghashariatmadari (2013)
            elif kCalcMethod=='vatankhah':
                kh=(k0h+k0h**2*np.exp(-(3.2+k0h**1.65)))/np.sqrt(np.tanh(k0h))+k0h*(1-np.exp(-(k0h**0.132)))**(5.0532+2.158*k0h**1.505) #Calculating wave number from Vatankhah and Aghashariatmadari (2013)
                k=kh/h #Calculating wave number from Vatankhah and Aghashariatmadari (2013)
                k[w==0]=0


            #Estimation of wave number (k) from Goad (2010)
            elif kCalcMethod=='goda':
                kh=np.zeros(len(k0h))
                kh[k0h>=1]=k0h[k0h>=1]
                kh[k0h<1]=(k0h[k0h<1])**0.5
                for i in range(0,3,1):
                    kh=kh-((kh-k0h*(np.tanh(kh))**-1)/(1+k0h*((np.tanh(kh))**(-2)-1))) #Calculating wave number from Goda (2010)

                k=kh/h #Calculating wave number from Goda (2010)
                k[w==0]=0
            
            
            #Calculating exact wave number (k)
            elif kCalcMethod=='exact':
            
                #Estimation of wave number (k) from Beji (2013)
                kh=k0h*(1+k0h**1.09*np.exp(-(1.55+1.3*k0h+0.216*k0h**2)))/np.sqrt(np.tanh(k0h)) #Calculating wave number from Beji (2013)
                kini=kh/h #Initial value for k (Wave number from Beji (2013))
                kini[w==0]=0
            
                #Calculating exact wave number (k)
                k=np.zeros(len(f)) #Pre-assigning array to make program run faster
                for i in range(0,len(f),1):
                    if len(h)==1:
                        fun = lambda x : w[i]**2-(9.81*x*np.tanh(x*h)) #w**2=g*k*tanh(kh)
                    else:
                        fun = lambda x : w[i]**2-(9.81*x*np.tanh(x*h[i])) #w**2=g*k*tanh(kh)
            
                    k[i]=sp.optimize.fsolve(fun,kini[i]) #Wave number
        
            PHI=(np.tanh(k*h))**2/(1+(2*k*h)/(np.sinh(2*k*h)))
            PHI[0]=0
        
        #Transformation function from JONSWAP into TMA, exact method (Kitaigordskii et al., 1975)
        elif transfCalcMethod=='kitaigordskii':
            Rwh=np.ones(len(omega))
            for i in range(len(f)-1,0,-1):
                fun=lambda x : (x*np.tanh(omega[i]**2*x)-1) #R(wh).tanh(wh^2*R(wh))=1
                if i==len(f)-1:
                    Rwh[i]=sp.optimize.fsolve(fun,0)
                else:    
                    Rwh[i]=sp.optimize.fsolve(fun,Rwh[i+1])
        
            Rwh[0]=1
            PHI=Rwh**-2*(1+(2*omega**2*Rwh)/(np.sinh(2*omega**2*Rwh)))**-1
            PHI[0]=0

    #--------------------------------------------------------------------------
    #Applying diagnostic tail

    #Index of ftail
    Indxftail=int(np.min((np.nonzero(f>=ftail))[0])) #Indx=Indx.item(0)# Indx=Indx[0]; Indx=Indx.item(0) #convert tulip to numpy array
    
    #Applying JONSWAP Spectrum tail after ftail
    if tailtype=='jonswap':
    
        Syy[f>ftail]=Syy[Indxftail]*(f[f>ftail]/ftail)**tailpower #Adding diagnostic tail
        Syy[Syy<0]=0 #Syy can not be negative
    
    #Applying TMA Spectrum tail after ftail
    elif tailtype=='tma':
    
        Syy[f>ftail]=Syy[Indxftail]*(PHI[f>ftail]/PHI[Indxftail])*(f[f>ftail]/ftail)**tailpower #Adding TMA Spectrum tail
        Syy[Syy<0]=0 #Syy can not be negative
    

    SyyOut=Syy.copy() #Assigning Syy to SyyOut
    fOut=f.copy() #Assigning f to fOut

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
        plt.plot(f,SyyIn,label='Input')
        #plt.hold='True'
        plt.plot(f,SyyOut,label='Output')
    
        plt.title('Power Spectral Density')
        plt.xlabel('Frequency(Hz)')
        plt.ylabel('Spectral Density(m^2/Hz)')
        plt.legend()

    #--------------------------------------------------------------------------
    #Outputs
    return fOut, SyyOut, Hm0, fp, Tp, Tm01, Tm02

    #--------------------------------------------------------------------------