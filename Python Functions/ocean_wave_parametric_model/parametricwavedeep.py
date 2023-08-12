def parametricwavedeep(windvel, Fetch, CalcMethod='jonswap', dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2017-09-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.parametricwavedeep
    ==============================

    .. code:: python

        H, T, Ehat, fphat, m0, Fetchhat = parametricwavedeep(windvel, Fetch, CalcMethod='jonswap', dispout='no')

    Description
    -----------

    Calculate wave properties using parametric wave models in deep water

    Inputs
    ------

    windvel
        | Wind velocity in (m/s)
        | Wind velocity should be measured (or represents velocity) at 10 m above surface
        | Wind velocity should be a 10 minutes averaged wind for all methods, except for for 'cem' and 'spmdeep' methods
        | For 'cem' and 'spmdeep' methods, wind velocity should be converted to duration of sustained wind by using gust factor
    Fetch
        Wind fetch in (m)
    CalcMethod='jonswap'
        | Parametric wave model to be used for wave calculation 
        | 'wilson': Use method by Wislon (1965)
        | 'jonswap': Use method by Hasselmann et al. (1973) known as JONSWAP
        | 'spmdeep': Use method by Shore Protection Manual (SPM),
        |     U.S. Army Corps of Engineers (1984) in deep water
        | 'kahma': Use method by Kahma and Calkoen (1992)
        | 'hwang': Use method by Hwang and Wang (2004)
        | 'cem': Use method by Coastal Engineering Manual (CEM),
        |     U.S. Army Corps of Engineers (2015)
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    H
        | Predicted wave height in (m) 
        | For all methods excepth for 'wilson': H=Hm0 where, Hm0 is a zero-moment wave height
        | For 'wilson' method: H=Hs, where Hs is a significant wave height
    T
        | Predicted wave period in (s) 
        | For all methods excepth for 'wilson': T=Tp where, Tp is a peak wave period
        | For 'wilson' method: T=Ts, where Ts is a significant wave period (Ts=0.95Tp)
    Ehat
        Predicted dimensionless wave energy, Ehat=g^2*m0/U10^4
    fphat
        Predicted dimensionless peak wave frequency, fphat=fp*U10/g
    m0
        Zero-moment of water surface elevation power spectral density in (m^2)
    Fetchhat
        | Dimensionless Fetch: Fetchhat=g*Fetch/U10^2
        | Note, g=9.81: gravitational acceleration
        |     U10: wind velocity
        |     fp: peak wave frequency

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        windvel=10*np.random.rand(100)
        Fetch=10000*np.random.rand(100)
        H,T,Ehat,fphat,m0,Fetchhat=sm.parametricwavedeep(windvel,Fetch,'jonswap','no')

        windvel=10
        Fetch=np.arange(1e3,1e6,1000)
        H,T,Ehat,fphat,m0,Fetchhat=sm.parametricwavedeep(windvel,Fetch,'jonswap','yes')

    References
    ----------

    Department of the Army, Waterways Experiment Station, Corps of Engineers, 
    and Coastal Engineering Research Center (1984), 
    Shore Protection Manual, Washington, 
    D.C., vol. 1, 4th ed., 532 pp.

    Hasselmann, K.; Barnett, T. P.; Bouws, E.; Carlson, H.; Cartwright, D. E.; Enke, K.; Ewing, J. A.; 
    Gienapp, H.; Hasselmann, D. E.; Kruseman, P.; Meerbrug, A.; Muller, P.; Olbers, D. J.; Richter, K.; 
    Sell, W., and Walden, H., (1973). 
    Measurements of wind-wave growth and swell decay during the Joint North Sea Wave Project (JONSWAP). 
    Deutsche Hydrographische Zeitschrift A80(12), 95p.

    Holthuijsen, L. H. (2007). 
    Waves in oceanic and coastal waters. 
    Cambridge university press.

    Hwang, P. A., & Wang, D. W. (2004). 
    Field measurements of duration-limited growth of wind-generated ocean surface waves at young stage of development. 
    Journal of Physical Oceanography, 34(10), 2316-2326.

    Kahma, K. K., & Calkoen, C. J. (1992). 
    Reconciling discrepancies in the observed growth of wind-generated waves. 
    Journal of Physical Oceanography, 22(12), 1389-1405.

    Pierson, W. J., & Moskowitz, L. (1964). 
    A proposed spectral form for fully developed wind seas based on the similarity theory of SA Kitaigorodskii. 
    Journal of geophysical research, 69(24), 5181-5190.

    U.S. Army Corps of Engineers (2015). 
    Coastal Engineering Manual. 
    Engineer Manual 1110-2-1100, Washington, D.C.: U.S. Army Corps of Engineers.

    Wilson, B. W. (1965). 
    Numerical prediction of ocean waves in the North Atlantic for December, 1959. 
    Ocean Dynamics, 18(3), 114-130.

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
    
    windvel=type2numpy(windvel)
    Fetch=type2numpy(Fetch)

    #--------------------------------------------------------------------------
    #Calculating Wave Properties using parametric wave models in deep water

    #Calculating wave properties using Wislon (1965)
    if CalcMethod=='wilson':
    
        #Calculating dimensionless numbers from dimensional number
        Fetchhat=9.81*Fetch/(windvel**2) #Dimensionless Fetch: Fetchhat=g*Fetch/U10^2
    
        #Calculating dimensionless wave properties using Wislon (1965)
        Hshat=0.3*(1-(1+0.004*Fetchhat**(1/2))**(-2)) #Dimensionless significant wave height: Hshat=g*Hs/U10^2
        Tshat=(1.37*2*np.pi)*(1-(1+0.008*Fetchhat**(1/3))**(-5)) #Dimensionless significant wave period: Tshat=g*Ts/U10
    
        #Calculating Fully Developed Condition based on Pierson and Moskowitz (1964)
        #EhatFullyDev=3.64e-3 #Ehat for fully developed condition
        #m0FullyDev=EhatFullyDev*(windvel**4)/(9.81**2) #Zero-moment of power spectral density for fully developed condition
        #Hm0FullyDev=4*np.sqrt(m0FullyDev) #Zero-moment wave height for fully developed condition for fully developed condition
        #Hm0hatFullyDev=9.81*Hm0FullyDev/(windvel**2) #Dimensionless zero-moment wave height for fully developed condition 
        Hm0hatFullyDev=4*np.sqrt(3.64e-3) #Dimensionless zero-moment wave height for fully developed condition 
    
        #fphatFullyDev=0.133 #fphat for fully developed condition
        #fpFullyDev=fphatFullyDev*9.81/windvel #Peak wave frequency for fully developed condition
        #TpFullyDev=1/fpFullyDev #Peak wave period for fully developed condition
        #TsFullyDev=0.95*TpFullyDev #Significant wave period for fully developed condition 
        #TshatFullyDev=9.81*TsFullyDev/windvel #Dimensionless significant wave period for fully developed condition
        TshatFullyDev=0.95*(1/0.133) #Dimensionless significant wave period for fully developed condition
    
        #Checking the predicted values versus fully Deveoped Condition
        Hshat[Hshat>Hm0hatFullyDev]=Hm0hatFullyDev #Check the predicted value be smaller than the ones for fully deveoped condition
        Tshat[Tshat>TshatFullyDev]=TshatFullyDev #Check the predicted value be smaller than the ones for fully deveoped condition
    
        #Calculating dimensional wave propoerties from dimensionless number
        Hs=Hshat*(windvel**2)/9.81 #Significant wave height: Hshat=g*Hs/U10^2
        Ts=Tshat*windvel/9.81 #Significant wave period, Tshat=g*Ts/U10
    
        #Calculating dimensionless wave propoerties
        Hm0=Hs #Zero-moment wave height, Hm0=Hs
        m0=(Hm0/4)**2 #Zero-moment of power spectral density, Hm0=4*(m0)^0.5
        Ehat=9.81**2*m0/(windvel**4) #Dimensionless wave energy: Ehat=g^2*m0/U10^4
    
        Tp=Ts/0.95 #Peak wave period, Ts=0.95*Tp
        fp=1/Tp #Peak wave frequency, fp=1/Tp
        fphat=fp*windvel/9.81 #Dimensionless peak wave frequency: fphat=fp*U10/g
    
        #Assigning wave propoerties
        H=Hs.copy() #Significant wave height
        T=Ts.copy() #Significant wave period
    
    
    #Calculating wave properties using JONSWAP, Hasselmann et al. (1973)
    elif CalcMethod=='jonswap':
    
        #Calculating dimensionless numbers from dimensional number
        Fetchhat=9.81*Fetch/(windvel**2) #Dimensionless Fetch: Fetchhat=g*Fetch/U10^2
    
        #Calculating dimensionless wave properties using JONSWAP
        Ehat=1.6e-7*Fetchhat**(1) #Dimensionless wave energy: Ehat=g^2*m0/U10^4
        fphat=3.5*Fetchhat**(-0.33) #Dimensionless peak wave frequency: fphat=fp*U10/g
    
        #Calculating Fully Developed Condition based on Pierson and Moskowitz (1964)
        EhatFullyDev=3.64e-3 #Ehat for fully developed condition
        fphatFullyDev=0.133 #fphat for fully developed condition
    
        #Checking the predicted values versus fully Deveoped Condition
        Ehat[Ehat>EhatFullyDev]=EhatFullyDev #Check the predicted value be smaller than the ones for fully deveoped condition
        fphat[fphat<fphatFullyDev]=fphatFullyDev #Check the predicted value be larger than the ones for fully deveoped condition
    
        #Calculating dimensional wave propoerties from dimensionless number
        m0=Ehat*(windvel**4)/(9.81**2) #Zero-moment of water surface elevation power spectral density, Ehat=g^2*m0/U10^4
        Hm0=4*np.sqrt(m0) #Zero-moment wave height, Hm0=4*(m0)^0.5
    
        fp=fphat*9.81/(windvel) #Peak wave frequency, fphat=fp*U10/g
        Tp=1/fp #Peak wave period, Tp=1/fp
    
        #Assigning wave propoerties
        H=Hm0.copy() #Peak wave height
        T=Tp.copy() #Peak wave period
    
    
    #Calculating wave properties using Shore Protection Manual (SPM), U.S. Army Corps of Engineers (1984) in deep water
    elif CalcMethod=='spmdeep':
    
        #Calculating adjusted wind velcity
        UA=0.71*windvel**1.23 #Wind stress factor or adjusted wind velcity
        
        #Calculating dimensionless numbers from dimensional number
        Fetchhat=9.81*Fetch/(windvel**2) #Dimensionless Fetch: Fetchhat=g*Fetch/U10^2
        FetchhatUA=9.81*Fetch/(UA**2) #Dimensionless Fetch, FetchhatUA=g.Fetch/UA^2
    
        #Calculating dimensionless wave properties using SPM (1984)
        Hm0hat=1.6e-3*(FetchhatUA)**(1/2) # Dimensionless zero-moment wave height based on UA: Hm0hat=g*Hm0/UA^2
        Tphat=2.857e-1*(FetchhatUA)**(1/3) # Dimensionless peak wave period based on UA: Tphat=g*TpUA
    
        #Calculating Fully Developed Condition based on SPM (1984)
        Hm0hatFullyDev=2.433e-1 #Fully Developed Condition from SPM (1984) based on UA, Holthuijsen(2007), Hm0hat*=(gHm0)/(UA^2)=2.433*10^-1
        #Hm0FullyDev=Hm0hatFullyDev*(UA**2)/9.81 #Significant wave height for fully Developed Condition
        #m0FullyDev=(Hm0FullyDev**2)/16 #Zero-moment of surface elevation power spectral density for fully Developed Condition, Hm0=4*(m0)^0.5
        #EhatFullyDev=9.81**2*m0FullyDev/(windvel**4) #Dimensionless wave energy for fully Developed Condition: Ehat=g^2*m0/U10^4
    
        TphatFullyDev=8.134 #Fully Developed Condition from SPM (1984), Holthuijsen(2007), Tphat*=(gTp)/UA=8.134
        #TpFullyDev=TphatFullyDev*UA/9.81 #Peak wave period for fully Developed Condition
        #fpFullyDev=1/TpFullyDev #Peak wave frequency for fully Developed Condition: fp=1/Tp
        #fphatFullyDev=fpFullyDev*windvel/9.81 #Dimensionless peak wave frequency for fully Developed Condition: fphat=fp*U10/g
    
        #Checking the predicted values versus fully Deveoped Condition
        Hm0hat[Hm0hat>Hm0hatFullyDev]=Hm0hatFullyDev #Check the predicted value be smaller than the ones for fully deveoped condition
        Tphat[Tphat>TphatFullyDev]=TphatFullyDev #Check the predicted value be smaller than the ones for fully deveoped condition
    
        #Calculating dimensional wave propoerties from dimensionless number
        Hm0=Hm0hat*(UA**2)/9.81 #Zero-moment wave height: Hm0hat=g*Hm0/UA^2
        Tp=Tphat*UA/9.81 #Peak wave period: Tphat=g*Tp/UA
    
        #Calculating dimensionless wave propoerties
        m0=(Hm0**2)/16 #Zero-moment of surface elevation power spectral density, Hm0=4*(m0)^0.5
        Ehat=9.81**2*m0/(windvel**4) #Dimensionless wave energy: Ehat=g^2*m0/U10^4
    
        fp=1/Tp #Peak wave frequency: fp=1/Tp
        fphat=fp*windvel/9.81 #Dimensionless peak wave frequency: fphat=fp*U10/g
    
        #Assigning wave propoerties
        H=Hm0.copy() #Peak wave height
        T=Tp.copy() #Peak wave period
    
    
    #Calculating wave properties using Kahma and Calkoen (1992) 
    elif CalcMethod=='kahma':
    
        #Calculating dimensionless numbers from dimensional number
        Fetchhat=9.81*Fetch/(windvel**2) #Dimensionless Fetch: Fetchhat=g*Fetch/U10^2
    
        #Calculating dimensionless wave properties using Kahma and Calkoen (1992)
        Ehat=5.2e-7*Fetchhat**(0.9) #Dimensionless wave energy: Ehat=g^2*m0/U10^4
        fphat=(13.7/(2*np.pi))*Fetchhat**(-0.27) #Dimensionless peak wave frequency: fphat=fp*U10/g
    
        #Calculating Fully Developed Condition based on Pierson and Moskowitz (1964)
        EhatFullyDev=3.64e-3 #Ehat for fully developed condition
        fphatFullyDev=0.133 #fphat for fully developed condition
    
        #Checking the predicted values versus fully Deveoped Condition
        Ehat[Ehat>EhatFullyDev]=EhatFullyDev #Check the predicted value be smaller than the ones for fully deveoped condition
        fphat[fphat<fphatFullyDev]=fphatFullyDev #Check the predicted value be larger than the ones for fully deveoped condition
    
        #Calculating dimensional wave propoerties from dimensionless number
        m0=Ehat*(windvel**4)/(9.81**2) #Zero-moment of water surface elevation power spectral density, Ehat=g^2*m0/U10^4
        Hm0=4*np.sqrt(m0) #Zero-moment wave height, Hm0=4*(m0)^0.5
    
        fp=fphat*9.81/(windvel) #Peak wave frequency, fphat=fp*U10/g
        Tp=1/fp #Peak wave period, Tp=1/fp
    
        #Assigning wave propoerties
        H=Hm0.copy() #Peak wave height
        T=Tp.copy() #Peak wave period
    
    
    #Calculating wave properties using Hwang and Wang (2004)
    elif CalcMethod=='hwang':
    
        #Calculating dimensionless numbers from dimensional number
        Fetchhat=9.81*Fetch/(windvel**2) #Dimensionless Fetch: Fetchhat=g*Fetch/U10^2
    
        #Calculating dimensionless wave properties using Hwang and Wang (2004)
        Ehat=6.191e-7*Fetchhat**(0.8106) #Dimensionless wave energy: Ehat=g^2*m0/U10^4
        fphat=(11.86/(2*np.pi))*Fetchhat**(-0.2368) #Dimensionless peak wave frequency: fphat=fp*U10/g
    
        #Calculating Fully Developed Condition based on Pierson and Moskowitz (1964)
        EhatFullyDev=3.64e-3 #Ehat for fully developed condition
        fphatFullyDev=0.133 #fphat for fully developed condition
    
        #Checking the predicted values versus fully Deveoped Condition
        Ehat[Ehat>EhatFullyDev]=EhatFullyDev #Check the predicted value be smaller than the ones for fully deveoped condition
        fphat[fphat<fphatFullyDev]=fphatFullyDev #Check the predicted value be larger than the ones for fully deveoped condition
    
        #Calculating dimensional wave propoerties from dimensionless number
        m0=Ehat*(windvel**4)/(9.81**2) #Zero-moment of water surface elevation power spectral density, Ehat=g^2*m0/U10^4
        Hm0=4*np.sqrt(m0) #Zero-moment wave height, Hm0=4*(m0)^0.5
    
        fp=fphat*9.81/(windvel) #Peak wave frequency, fphat=fp*U10/g
        Tp=1/fp #Peak wave period, Tp=1/fp
    
        #Assigning wave propoerties
        H=Hm0.copy() #Peak wave height
        T=Tp.copy() #Peak wave period
    
    
    #Calculating wave properties using Coastal Engineering Manual (CEM), U.S. Army Corps of Engineers (2015)
    elif CalcMethod=='cem':
    
        #Calculating shear velcity
        CD=0.001*(1.1+0.035*windvel) #Drag Coefficient
        Ustar=(CD*(windvel**2))**0.5  #Shear Velocity in (m/s)
        
        #Calculating dimensionless numbers from dimensional number
        Fetchhat=9.81*Fetch/(windvel**2) #Dimensionless Fetch: Fetchhat=g*Fetch/U10^2
        FetchhatStar=9.81*Fetch/(Ustar**2) #Dimensionless Fetch, FetchhatStar=g.Fetch/UA^2
    
        #Calculating dimensionless wave properties using CEM (2015)
        Hm0hat=4.13e-2*FetchhatStar**(1/2) #Dimensionless zero-moment wave height based on shear velocity: Hm0hat=g*Hm0/U^*2
        Tphat=0.651*FetchhatStar**(1/3) #Dimensionless peak wave period based on shear velocity: Tphat=g*Tp/U*
    
        #Calculating Fully Developed Condition based on CEM (2015)
        Hm0hatFullyDev=2.115e2 #Fully developed condition for Hm0hat=g*Hm0/U^*2
        #Hm0FullyDev=Hm0hatFullyDev*(Ustar**2)/9.81 #Zero-moment wave height: Hm0hat=g*Hm0/U*^2
        #m0FullyDev=(Hm0FullyDev**2)/16 #Zero-moment of surface elevation power spectral density, Hm0=4*(m0)^0.5
        #EhatFullyDev=9.81**2*m0FullyDev/(windvel**4) #Dimensionless wave energy: Ehat=g^2*m0/U10^4
    
        TphatFullyDev=2.398e2 #Fully developed condition for Tphat=g*Tp/U*
        #TpFullyDev=TphatFullyDev*Ustar/9.81 #Peak wave period: Tphat=g*Tp/U*
        #fpFullyDev=1/TpFullyDev #Peak wave frequency: fp=1/Tp
        #fphatFullyDev=fpFullyDev*windvel/9.81 #Dimensionless peak wave frequency: fphat=fp*U10/g
    
        #Checking the predicted values versus fully Deveoped Condition
        Hm0hat[Hm0hat>Hm0hatFullyDev]=Hm0hatFullyDev #Check the predicted value be smaller than the ones for fully deveoped condition
        Tphat[Tphat>TphatFullyDev]=TphatFullyDev #Check the predicted value be smaller than the ones for fully deveoped condition
    
        #Calculating dimensional wave propoerties from dimensionless number
        Hm0=Hm0hat*(Ustar**2)/9.81 #Zero-moment wave height: Hm0hat=g*Hm0/U^*2
        Tp=Tphat*Ustar/9.81 #Peak wave period: Tphat=g*Tp/U*
    
        #Calculating dimensionless wave propoerties
        m0=(Hm0**2)/16 #Zero-moment of surface elevation power spectral density, Hm0=4*(m0)^0.5
        Ehat=9.81**2*m0/(windvel**4) #Dimensionless wave energy: Ehat=g^2*m0/U10^4
    
        fp=1/Tp #Peak wave frequency: fp=1/Tp
        fphat=fp*windvel/9.81 #Dimensionless peak wave frequency: fphat=fp*U10/g
    
        #Assigning wave propoerties
        H=Hm0.copy() #Peak wave height
        T=Tp.copy() #Peak wave period
    
    
    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':

        #Plotting
        plt.subplot(1,2,1)
        plt.loglog(Fetchhat,Ehat)

        plt.xlabel('Dimensionless Wind Fetch, Fetchhat=g*Fetch/U10^2')
        plt.ylabel('Dimensionless wave energy, Ehat=g^2*m0/U10^4')

        plt.subplot(1,2,2)
        plt.loglog(Fetchhat,fphat)

        plt.xlabel('Dimensionless Wind Fetch, Fetchhat=g*Fetch/U10^2')
        plt.ylabel('Dimensionless peak wave frequency, fphat=fp*U10/g')

    #--------------------------------------------------------------------------
    #Outputs
    return H, T, Ehat, fphat, m0, Fetchhat

    #--------------------------------------------------------------------------
