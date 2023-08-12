def wavebedstress(h, heightfrombed, d50, Uc=0, H=0, T=0, Rho=1000, kCalcMethod='beji'):
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

    scientimate.wavebedstress
    =========================

    .. code:: python

        Tau, Tauc, Tauw, Ust, z0 = scientimate.wavebedstress(h, heightfrombed, d50, Uc=0, H=0, T=0, Rho=1000, kCalcMethod='beji')

    Description
    -----------

    Calculate the bottom shear velocity and shear stress from current velocity and wave 

    Inputs
    ------

    h
        Water depth in (m) 
    heightfrombed
        Sensor height from bed in (m)
    d50
        Median bed particle diameter in (m)
    Uc=0
        Current velocity at sensor height in (m/s), set equal to 0 if not exist 
    H=0
        | Wave height in (m), H can be approximated and replaced by Hrms for random waves
        | Hrms: root mean square wave height, Hrms=Hm0/sqrt(2)=Hs/sqrt(2) 
        | Hm0: zero moment wave height, Hs: significant wave height
        | set equal to 0 if not exist  
    T=0
        | Wave period in (s), T can be approximated and replaced by mean wave period for random waves 
        | if peak wave frequency (Tp) is used, calculated values represent peak wave 
        | set equal to 0 if not exist  
    Rho=1000
        Water density in (kg/m^3)
    kCalcMethod='beji'
        | Wave number calculation method 
        | 'hunt': Hunt (1979), 'beji': Beji (2013), 'vatankhah': Vatankhah and Aghashariatmadari (2013) 
        | 'goda': Goda (2010), 'exact': calculate exact value 
        | Note: inputs can be as a single value or a 1-D vertical array

    Outputs
    -------

    Tau
        Total bottom shear stress from current velocity and wave
    Tauc
        Bottom shear stress from current velocity
    Tauw
        Bottom shear stress from wave
    Ust
        Bottom shear velocity (m/s)
    z0
        Surface roughness in (N/m^2)

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        Tau,Tauc,Tauw,Ust,z0=sm.wavebedstress(2,1.1,0.0000245,1.5,0.5,3,1000,'beji')

        Tau,Tauc,Tauw,Ust,z0=sm.wavebedstress(2,1.1,0.0000245,[1.5,2],[0.5,0.6],[3,3.1],1000,'exact')

        Tau,Tauc,Tauw,Ust,z0=sm.wavebedstress(2,1.1,0.0000245,np.array([1.5,2]),np.array([0.5,0.6]),np.array([3,3.1]),1000,'exact')

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
    heightfrombed=type2numpy(heightfrombed)
    d50=type2numpy(d50)
    Uc=type2numpy(Uc)
    H=type2numpy(H)
    T=type2numpy(T)
    Rho=type2numpy(Rho)

    #--------------------------------------------------------------------------
    #Calculating bed shear stress due to bed current (Tauc)

    nu=1e-6 #Water viscosity
    ks=2.5*d50 #Nikuradse Equivalent surface roughness, ks=2.5*d50  
    zsensor=heightfrombed #Height that velocity measured at
    
    z0=np.zeros(len(Uc)) 
    Ust=np.zeros(len(Uc)) 
    for i in range(0,len(Uc),1):
        
        z01=0
        z02=ks/30 #Initial guess for z0
        Ust1=0
        Ust2=(0.4*Uc[i])/(np.log(zsensor/z02)) #Initial guess for U*
        
        while np.abs(z02-z01)>0.000001 or np.abs(Ust2-Ust1)>0.000001:  #Trial and error to calculate z0 and U*
            z01=z02
            Ust1=Ust2
            z02=ks/30*(1-np.exp(-Ust2*ks/(27*nu)))+nu/(9*Ust2) #Calculatin z0 from Christoffersen and Jonsson (1985)
            Ust2=(0.4*Uc[i])/(np.log(zsensor/z02)) #calculating friction velocity U*, from logarithmic velocity profile
        
        z0[i]=z02 #Surface roughness
        Ust[i]=Ust2 #Shear velocity
        
    Tauc=Rho*(Ust**2) #Calculating bed shear stress due to current

    #--------------------------------------------------------------------------
    #Calculating wave number (k)

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

    #--------------------------------------------------------------------------
    #Calculating bed shear stress due to wave (Tauw)

    Uwb=(np.pi*H)/(T*np.sinh(k*h)) #Wave orbital veocity at bed using linear wave theory
    A=Uwb*T/(2*np.pi) #Wave orbital excursion at the bed using linear wave theory
    
    fw=np.zeros(len(T)) 
    for i in range(0,len(T),1):
        if A[i]/ks > 1.57:
            fw[i]=0.00251*np.exp(5.21*((A[i]/ks)**(-0.19))) #Calculating wave friction factor fw (Swart, 1974)
        elif A[i]/ks <= 1.57:
            fw[i]=0.3
    
    Tauw=0.5*Rho*(fw*(Uwb**2)) #Bed shear stress due to wave 

    #--------------------------------------------------------------------------
    #Calculating Total bottom shear stress from wave and current

    Tau=Tauc*(1+1.2*((Tauw/(Tauc+Tauw))**3.2))+Tauw #Total bottom shear stress (Soulsby, 1997)

    #--------------------------------------------------------------------------
    #Outputs
    return Tau, Tauc, Tauw, Ust, z0

    #--------------------------------------------------------------------------
