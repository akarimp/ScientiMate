def parametricwaveshallow(windvel, Fetch, hmean, CalcMethod='karimpour', dispout='no'):
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

    scientimate.parametricwaveshallow
    =================================

    .. code:: python

        H, T, Ehat, fphat, m0, Fetchhat, hhat = scientimate.parametricwaveshallow(windvel, Fetch, hmean, CalcMethod='karimpour', dispout='no')

    Description
    -----------

    Calculate wave properties using parametric wave models in shallow and intermediate water

    Inputs
    ------

    windvel
        | Wind velocity in (m/s)
        | Wind velocity should be measured (or represents velocity) at 10 m above surface
        | Wind velocity should be a 10 minutes averaged wind for all methods, except for for 'spmshallow' methods
        | For 'spmshallow' methods, wind velocity should be converted to duration of sustained wind by using gust factor
    Fetch
        Wind fetch in (m)
    hmean
        Mean water depth along a wind fetch in (m)
    CalcMethod='karimpour'
        | Parametric wave model to be used for wave calculation 
        | 'spmshallow': Use method by Shore Protection Manual (SPM),
        |     U.S. Army Corps of Engineers (1984) in shallow and intermediate water water
        | 'young': Use method by Young and Verhagen (1996) and Young and Babanin (2006)
        | 'karimpour': Use method by Karimpour et al. (2017)
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    H
        | Predicted wave height in (m) 
        | For all methods: H=Hm0 where, Hm0 is a zero-moment wave height
    T
        | Predicted wave period in (s) 
        | For all methods excepth for 'spmshallow': T=Tp where, Tp is a peak wave period
        | For 'spmshallow' method: T=Ts, where Ts is a significant wave period (Ts=0.95Tp)
    Ehat
        Predicted dimensionless wave energy, Ehat=g^2*m0/U10^4
    fphat
        Predicted dimensionless peak wave frequency, fphat=fp*U10/g
    m0
        Zero-moment of water surface elevation power spectral density in (m^2)
    Fetchhat
        Dimensionless wind fetch: Fetchhat=g*Fetch/U10^2
    hhat
        | Dimensionless mean water depth along a wind fetch: hhat=g*hmean/U10^2
        | Note, g=9.81: gravitational acceleration
        |     U10: wind velocity
        |     fp: peak wave frequency
        |     hmean: mean depth along a wind fetch

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        windvel=10*np.random.rand(100)
        Fetch=10000*np.random.rand(100)
        hmean=3*np.random.rand(100)
        H,T,Ehat,fphat,m0,Fetchhat,hhat=sm.parametricwaveshallow(windvel,Fetch,hmean,'karimpour','no')

        windvel=10
        Fetch=np.arange(1e3,1e6,1000)
        hmean=3
        H,T,Ehat,fphat,m0,Fetchhat,hhat=sm.parametricwaveshallow(windvel,Fetch,hmean,'karimpour','yes')

    References
    ----------

    Bretschneider, C. L. (1952). 
    Revised wave forecasting relationships. 
    Coastal Engineering Proceedings, 1(2), 1.

    Bretschneider, C. L. (1958). 
    Revisions in wave forecasting: deep and shallow water. 
    Coastal Engineering Proceedings, 1(6), 3.

    Department of the Army, Waterways Experiment Station, Corps of Engineers, 
    and Coastal Engineering Research Center (1984), 
    Shore Protection Manual, Washington, 
    D.C., vol. 1, 4th ed., 532 pp.

    Karimpour, A., Chen, Q., & Twilley, R. R. (2017). 
    Wind Wave Behavior in Fetch and Depth Limited Estuaries. 
    Scientific reports, 7, 40654.

    Pierson, W. J., & Moskowitz, L. (1964). 
    A proposed spectral form for fully developed wind seas based on the similarity theory of SA Kitaigorodskii. 
    Journal of geophysical research, 69(24), 5181-5190.

    Sverdrup, H. U., & Munk, W. H. (1947). 
    Wind, sea, and swell: theory of relations for forecasting. 
    U.S. Navy Department, Hydrographic Office, Publication No. 601, 44 pp. 

    Young, I. R., & Verhagen, L. A. (1996). 
    The growth of fetch limited waves in water of finite depth. Part 1. Total energy and peak frequency. 
    Coastal Engineering, 29(1-2), 47-78.

    Young, I. R., & Babanin, A. V. (2006). 
    The form of the asymptotic depth‐limited wind wave frequency spectrum. 
    Journal of Geophysical Research: Oceans, 111(C6).

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
    hmean=type2numpy(hmean)

    #--------------------------------------------------------------------------
    #Calculating Wave Properties using parametric wave models in shallow and intermediate water

    #Calculating wave properties using Shore Protection Manual (SPM), U.S. Army Corps of Engineers (1984) in shallow and intermediate water
    if CalcMethod=='spmshallow':
    
        #Calculating adjusted wind velcity
        UA=0.71*windvel**1.23 #Wind stress factor or adjusted wind velcity
        
        #Calculating dimensionless numbers from dimensional number
        Fetchhat=9.81*Fetch/(windvel**2) #Dimensionless Fetch: Fetchhat=g*Fetch/U10^2
        hhat=9.81*hmean/(windvel**2) #Dimensionless depth: hhat=g*hmean/U10^2
        FetchhatUA=9.81*Fetch/(UA**2) #Dimensionless Fetch, FetchhatUA=g.Fetch/UA^2
        hhatUA=9.81*hmean/(UA**2) #Dimensionless Depth, hhatUA=g.hmean/UA^2
    
        #Calculating dimensionless wave properties using SPM (1984)
        Hm0hat=0.283*np.tanh(0.53*(hhatUA)**(3/4))*np.tanh((0.00565*(FetchhatUA)**(1/2))/(np.tanh(0.53*(hhatUA)**(3/4)))) #Dimensionless zero-moment wave height: Hm0hat=Hm0*g/UA^2
        Tshat=7.54*np.tanh(0.833*(hhatUA)**(3/8))*np.tanh((0.0379*(FetchhatUA)**(1/3))/(np.tanh(0.833*(hhatUA)**(3/8)))) #Dimensionless Significant Wave Period, Tshat=Ts*g/UA
    
        #Calculating Dimensionless Asymptotic Limits using SMB method (Munk, 1947; and Bretschneider, 1952, 1958)
        EhatAsympt=1.4e-3*(9.81*hmean/(windvel**2))**1.5 #Asymptotic Limit from SMB
        m0Asympt=EhatAsympt*(windvel**4)/(9.81**2) #Zero-moment of surface elevation power spectral density, Ehat=g^2*m0/U10^4
        Hm0Asympt=4*np.sqrt(m0Asympt) #Zero-moment wave height, Hm0=4*(m0)^0.5
        Hm0hatAsympt=9.81*Hm0Asympt/(UA**2) #Dimensionless zero-moment wave height: Hm0hat=g*Hm0/U10^2
    
        fphatAsympt=0.16*(9.81*hmean/(windvel**2))**-0.375 #Asymptotic Limit from SMB
        fpAsympt=fphatAsympt*9.81/(windvel) #Peak wave frequency, fphat=fp*U10/g
        TpAsympt=1/fpAsympt #Peak Wave Period, Tp=1/fp
        TsAsympt=0.95*TpAsympt #Peak Wave Period (Ts=0.95Tp for Windsea and Ts=Tp for Swell)
        TshatAsympt=TsAsympt*9.81/UA #Dimensionless peak wave frequency: fphat=fp*U10/g
    
        #Check the predicted values versus Asymptotic Limits
        if len(Hm0hatAsympt)==1:
            Hm0hat[Hm0hat>Hm0hatAsympt]=Hm0hatAsympt #Check the predicted value be smaller than the ones for Asymptotic Limits
            Tshat[Tshat>TshatAsympt]=TshatAsympt #Check the predicted value be smaller than the ones for Asymptotic Limits
        else:
            Hm0hat[Hm0hat>Hm0hatAsympt]=Hm0hatAsympt[Hm0hat>Hm0hatAsympt] #Check the predicted value be smaller than the ones for Asymptotic Limits
            Tshat[Tshat>TshatAsympt]=TshatAsympt[Tshat>TshatAsympt] #Check the predicted value be smaller than the ones for Asymptotic Limits
    
        #Calculating Fully Developed Condition based on SPM (1984)
        Hm0hatFullyDev=2.433e-1 #Fully Developed Condition from SPM (1984) based on UA, Holthuijsen(2007), Hm0hat*=(gHm0)/(UA^2)=2.433*10^-1
        #Hm0FullyDev=Hm0hatFullyDev*(UA**2)/9.81 #Significant wave height for fully Developed Condition
        #m0FullyDev=(Hm0FullyDev**2)/16 #Zero-moment of surface elevation power spectral density for fully Developed Condition, Hm0=4*(m0)^0.5
        #EhatFullyDev=9.81**2*m0FullyDev/(windvel**4) #Dimensionless wave energy for fully Developed Condition: Ehat=g^2*m0/U10^4
    
        TphatFullyDev=8.134 #Fully Developed Condition from SPM (1984), Holthuijsen(2007), Tphat*=(gTp)/UA=8.134
        TshatFullyDev=0.95*8.134 #Fully Developed Condition from SPM (1984), Ts=0.95Tp
        #TpFullyDev=TphatFullyDev*UA/9.81 #Peak wave period for fully Developed Condition
        #fpFullyDev=1/TpFullyDev #Peak wave frequency for fully Developed Condition: fp=1/Tp
        #fphatFullyDev=fpFullyDev*windvel/9.81 #Dimensionless peak wave frequency for fully Developed Condition: fphat=fp*U10/g
    
        #Checking the predicted values versus fully Deveoped Condition
        Hm0hat[Hm0hat>Hm0hatFullyDev]=Hm0hatFullyDev #Check the predicted value be smaller than the ones for fully deveoped condition
        Tshat[Tshat>TshatFullyDev]=TshatFullyDev #Check the predicted value be smaller than the ones for fully deveoped condition
    
        #Calculating dimensional wave propoerties from dimensionless number
        Hm0=Hm0hat*(UA**2)/9.81 #Zero-moment wave height: Hm0hat=g*Hm0/UA^2
        Ts=Tshat*UA/9.81 #Significant Wave Period, Tshat=Ts*g/UA
    
        #Calculating dimensionless wave propoerties
        m0=(Hm0**2)/16 #Zero-moment of surface elevation power spectral density, Hm0=4*(m0)^0.5
        Ehat=9.81**2*m0/(windvel**4) #Dimensionless wave energy: Ehat=g^2*m0/U10^4
    
        Tp=Ts/0.95 #Peak Wave Period (Ts=0.95Tp for Windsea and Ts=Tp for Swell)
        fp=1/Tp #Peak Wave Frequency, fp=1/Tp
        fphat=fp*windvel/9.81 #Dimensionless Peak Wave Frequency, fphat=fp*U10/g
    
        #Assigning wave propoerties
        H=Hm0.copy() #Peak wave height
        T=Ts.copy() #Peak wave period
    
    
    #Calculating wave properties using Young and Verhagen (1996) and Young and Babanin (2006)
    elif CalcMethod=='young':
    
        #Calculating dimensionless numbers from dimensional number
        Fetchhat=9.81*Fetch/(windvel**2) #Dimensionless Fetch: Fetchhat=g*Fetch/U10^2
        hhat=9.81*hmean/(windvel**2) #Dimensionless depth: hhat=g*hmean/U10^2
    
        #Calculating dimensionless wave properties using Young and Verhagen (1996)
        Ehat=3.64e-3*(np.tanh(0.493*(hhat)**(0.75))*np.tanh((0.00313*(Fetchhat)**(0.57))/(np.tanh(0.493*(hhat)**(0.75)))))**1.74 #Dimensionless wave energy: Ehat=g^2*m0/U10^4
        fphat=0.133*(np.tanh(0.331*(hhat)**(1.01))*np.tanh((5.215e-4*(Fetchhat)**(0.73))/(np.tanh(0.331*(hhat)**(1.01)))))**-0.37 #Dimensionless peak wave frequency: fphat=fp*U10/g
    
        #Calculating Dimensionless Asymptotic Limits using Young and Verhagen (1996) and Young and Babanin (2006)
        EhatAsympt=1e-3*hhat**1.2 #Asymptotic Limit from Young and Babanin (2006)
        fphatAsympt=0.2*hhat**-0.375 #Asymptotic Limit from Young and Verhagen (1996)
    
        #Check the predicted values versus Asymptotic Limits
        if len(EhatAsympt)==1:
            Ehat[Ehat>EhatAsympt]=EhatAsympt #Check the predicted value be smaller than the ones for Asymptotic Limits
            fphat[fphat<fphatAsympt]=fphatAsympt #Check the predicted value be larger than the ones for Asymptotic Limits
        else:
            Ehat[Ehat>EhatAsympt]=EhatAsympt[Ehat>EhatAsympt] #Check the predicted value be smaller than the ones for Asymptotic Limits
            fphat[fphat<fphatAsympt]=fphatAsympt[fphat<fphatAsympt] #Check the predicted value be larger than the ones for Asymptotic Limits
    
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
    
    
    #Calculating wave properties using Karimpour et al. (2017)
    elif CalcMethod=='karimpour':
    
        #Calculating dimensionless numbers from dimensional number
        Fetchhat=9.81*Fetch/(windvel**2) #Dimensionless Fetch: Fetchhat=g*Fetch/U10^2
        Fetchhat[Fetchhat>20000]=20000 #Predicted values for Fetchhat>20000 should be limited to values for Fetchhat=20000
        hhat=9.81*hmean/(windvel**2) #Dimensionless depth: hhat=g*hmean/U10^2
    
        #Calculating dimensionless wave properties using Karimpour et al. (2017)
        Ehat=3.0e-8*Fetchhat**(2.7*(Fetchhat/hhat)**-0.1) #Dimensionless wave energy: Ehat=g^2*m0/U10^4
        fphat=3.5*Fetchhat**(-0.75*(Fetchhat/hhat)**-0.1) #Dimensionless peak wave frequency: fphat=fp*U10/g
    
        #Calculating Dimensionless Asymptotic Limits using Karimpour et al. (2017)
        EhatAsympt=3e-5*np.tan(1.56255*(np.tanh(3.356*hhat))**0.315) #Asymptotic Limit for Ehat
        fphatAsympt=0.133*(np.tanh(0.832*hhat))**-0.716*(np.tanh(2.623*hhat))**0.461 #Asymptotic Limit for fphat
    
        #Check the predicted values versus Asymptotic Limits
        if len(EhatAsympt)==1:
            Ehat[Ehat>EhatAsympt]=EhatAsympt #Check the predicted value be smaller than the ones for Asymptotic Limits
            fphat[fphat<fphatAsympt]=fphatAsympt #Check the predicted value be larger than the ones for Asymptotic Limits
        else:
            Ehat[Ehat>EhatAsympt]=EhatAsympt[Ehat>EhatAsympt] #Check the predicted value be smaller than the ones for Asymptotic Limits
            fphat[fphat<fphatAsympt]=fphatAsympt[fphat<fphatAsympt] #Check the predicted value be larger than the ones for Asymptotic Limits
    
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
    
        #Recalculating Fetchhat to recover values that set to Fetchhat=20000
        Fetchhat=9.81*Fetch/(windvel**2) #Dimensionless Fetch: Fetchhat=g*Fetch/U10^2
    
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
    return H, T, Ehat, fphat, m0, Fetchhat, hhat

    #--------------------------------------------------------------------------
