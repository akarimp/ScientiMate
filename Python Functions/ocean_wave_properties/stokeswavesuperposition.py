def stokeswavesuperposition(h, a, T, Phi, fs, duration=10, kCalcMethod='beji', dispout='no'):
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

    scientimate.stokeswavesuperposition
    ===================================

    .. code:: python

        Eta, t, Etaij = scientimate.stokeswavesuperposition(h, a, T, Phi, fs, duration=10, kCalcMethod='beji', dispout='no')

    Description
    -----------

    Superposition second order stokes' waves

    Inputs
    ------

    h
        Mean water depth in (m)
    a
        Wave amplitude in (m)
    T
        Wave mean period in (s)
    Phi
        Phase (radian)
    fs=32
        Sample generation frequency (Hz), number of data points in one second
    duration=10
        Duration time that data will be generated in (s)
    kCalcMethod='beji'
        | Wave number calculation method 
        | 'hunt': Hunt (1979), 'beji': Beji (2013), 'vatankhah': Vatankhah and Aghashariatmadari (2013) 
        | 'goda': Goda (2010), 'exact': calculate exact value 
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    Eta
        Water Surface Level Time Series in (m)
    t
        Time in (s)
    Etaij
        Separated Water Surface Level Time Series in (m)

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np
    
        Eta,t,Etaij=sm.stokeswavesuperposition(5,[0.1,0.2,0.3,0.4],[1,1.5,2,2.5],[np.pi/2,np.pi/4,np.pi/16,np.pi/32],32,10,'beji','yes')
    
        Eta,t,Etaij=sm.stokeswavesuperposition(5,np.array([0.1,0.2,0.3,0.4]),np.array([1,1.5,2,2.5]),np.array([np.pi/2,np.pi/4,np.pi/16,np.pi/32]),32,10,'beji','yes')

    References
    ----------

    Beji, S. (2013). 
    Improved explicit approximation of linear dispersion relationship for gravity waves. 
    Coastal Engineering, 73, 11-12.

    Goda, Y. (2010). 
    Random seas and design of maritime structures. 
    World scientific.

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
    
    h=type2numpy(h)
    a=type2numpy(a)
    T=type2numpy(T)
    Phi=type2numpy(Phi)

    #--------------------------------------------------------------------------

    sample=fs*duration #number of sample in input file
    dt=1/fs
    t=np.linspace(dt,duration,sample) #time
    NoOfWave=len(a) #number of waves to be combined with each other
    H=2*a # Wave height
    
    f=1/T #Wave frequency
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
        
    L=2*np.pi/k
    
    # Generating Waves
    FirstOrder=np.zeros((len(t),NoOfWave)) #Pre-assigning array to make program run faster
    SecondOrder=np.zeros((len(t),NoOfWave)) #Pre-assigning array to make program run faster
    for i in range(0,NoOfWave,1):
        FirstOrder[:,i]=a[i]*np.cos(-w[i]*t+Phi[i])
        SecondOrder[:,i]=(np.pi*H[i]/8)*(H[i]/L[i])*(np.cosh(k[i]*h))*(2+np.cosh(2*k[i]*h))/((np.sinh(k[i]*h))**3)*(np.cos(-2*w[i]*t))

    Etaij=FirstOrder+SecondOrder

    Eta=np.sum(Etaij,axis=1)

#    FirstOrder=np.empty([len(t),0])
#    SecondOrder=np.empty([len(t),0])
#    for i in range(0,NoOfWave,1):
#        FirstOrder=np.column_stack((FirstOrder,(a[i])*np.cos(-w[i]*t+Phi[i])))
#        SecondOrder=np.column_stack((SecondOrder,(np.pi*H[i]/8)*(H[i]/L[i])*(np.cosh(k[i]*h))*(2+np.cosh(2*k[i]*h))/((np.sinh(k[i]*h))**3)*(np.cos(-2*w[i]*t))))
#        # SecondOrder(:,i)=((H(i,1))**2*k(i,1)/16)*(cosh(k(i,1)*h))*(2+cosh(2*k(i,1)*h))/((sinh(k(i,1)*h))**3)*(cos(-2*w(i,1)*t));
#
#    Etaij=FirstOrder+SecondOrder
#
#    Eta=np.empty(len(t))
#    for i in range(0,len(t)):
#        Eta[i]=np.sum(Etaij[i,0:])

    # -------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
        
        plt.subplot(2,1,1)
        for i in range(0,NoOfWave,1):
            plt.plot(t,Etaij[:,i])
            #plt.hold='True'
        plt.title('Generated Waves')
        plt.xlabel('Time (s)')
        plt.ylabel('Water Level (m)')
        plt.xlim(t[0],t[-1])
         
        plt.subplot(2,1,2)
        plt.plot(t,Eta,label='Water Level (m)')
        #plt.hold='True'
        plt.title('Generated Combined Wave')
        plt.xlabel('Time (s)')
        plt.ylabel('Water Level (m)')
        plt.xlim(t[0],t[-1])        

    # -------------------------------------------------------------------------
    #Outputs
    return Eta, t, Etaij

    # -------------------------------------------------------------------------
