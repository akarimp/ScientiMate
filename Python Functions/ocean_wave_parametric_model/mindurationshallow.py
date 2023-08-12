def  mindurationshallow(windvel, Fetch, hmean, CalcMethod='karimpour', dispout='no'):
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

    scientimate.mindurationshallow
    ==============================

    .. code:: python

        tmin, tminhat, Fetchhat, hhat = scientimate.mindurationshallow(windvel, Fetch, hmean, CalcMethod='karimpour', dispout='no')

    Description
    -----------

    Calculate a minimum required wind duration for wave to be fetch limited in shallow and intermediate water

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
        | Parametric wave model to be used for calculation 
        | 'spmshallow': Use method by Shore Protection Manual (SPM),
        |     U.S. Army Corps of Engineers (1984) in shallow and intermediate water water
        | 'karimpour': Use method by Karimpour et al. (2017)
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    tmin
        Minimum required wind duration for wind to be fetch limited in (second)
    tminhat
        Dimensionless minimum required wind duration for wind to be fetch limited: tminhat=g*tmin/U10
    Fetchhat
        Dimensionless wind fetch: Fetchhat=g*Fetch/U10^2
    hhat
        | Dimensionless mean water depth along a wind fetch: hhat=g*hmean/U10^2
        | Note, g=9.81: gravitational acceleration
        |     U10: wind velocity

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        windvel=10*np.random.rand(100)
        Fetch=10000*np.random.rand(100)
        hmean=3*np.random.rand(100)
        tmin,tminhat,Fetchhat,hhat=sm.mindurationshallow(windvel,Fetch,hmean,'karimpour','no')

        windvel=10
        Fetch=np.arange(1e3,1e6,1000)
        hmean=3
        tmin,tminhat,Fetchhat,hhat=sm.mindurationshallow(windvel,Fetch,hmean,'karimpour','yes')

    References
    ----------

    U.S. Army Corps of Engineers (2015). 
    Coastal Engineering Manual. 
    Engineer Manual 1110-2-1100, Washington, D.C.: U.S. Army Corps of Engineers.

    Department of the Army, Waterways Experiment Station, Corps of Engineers, 
    and Coastal Engineering Research Center (1984), 
    Shore Protection Manual, Washington, 
    D.C., vol. 1, 4th ed., 532 pp.

    Karimpour, A., Chen, Q., & Twilley, R. R. (2017). 
    Wind Wave Behavior in Fetch and Depth Limited Estuaries. 
    Scientific reports, 7, 40654.

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
    #Calculating a minimum required wind duration for wave to be fetch limited

    #Calculating a minimum required wind duration for wave to be fetch limited
    #using Shore Protection Manual (SPM), U.S. Army Corps of Engineers (1984) in shallow and intermediate water
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
        Ts=Tshat*UA/9.81 # Significant Wave Period, Tshat=Ts*g/UA
    
        tminstarhat=5.37e2*(9.81*Ts/UA)**(7/3) #Dimensionless minimum required wind duration: tminhat=g*tmin/U10
        tmin=tminstarhat*UA/9.81 #Minimum required wind duration in (s)
        tminhat=9.81*tmin/windvel #Dimensionless minimum required wind duration: tminhat=g*tmin/U10
    
    
    #Calculating a minimum required wind duration for wave to be fetch limited using Karimpour et al. (2017)
    elif CalcMethod=='karimpour':
    
        Fetchhat=9.81*Fetch/(windvel**2) #Dimensionless Fetch: Fetchhat=g*Fetch/U10^2
        Fetchhat[Fetchhat>20000]=20000 #Predicted values for Fetchhat>20000 should be limited to values for Fetchhat=20000
        hhat=9.81*hmean/(windvel**2) #Dimensionless depth: hhat=g*hmean/U10^2
    
        tminhat=2.59*(Fetchhat/hhat)**1*hhat**(2/3) #Dimensionless minimum required wind duration: tminhat=g*tmin/U10
    
        #Calculating a minimum required wind duration for wave to be fetch limited
        #using Coastal Engineering Manual (CEM), U.S. Army Corps of Engineers (2015) in deep water
        tmindeep=77.23*Fetch**0.67/(windvel**0.34*9.81**0.33) #Minimum required wind duration in (s)
        tminhatdeep=9.81*tmindeep/windvel #Dimensionless minimum required wind duration: tminhat=g*tmin/U10
    
        #Checking that tmin in shallow water be smaller than deep water
        if len(tminhatdeep)==1:
            tminhat[tminhat>tminhatdeep]=tminhatdeep
        else:
            tminhat[tminhat>tminhatdeep]=tminhatdeep[tminhat>tminhatdeep]
    
        #Calculating a minimum required wind duration for wave to be fetch limited
        #using Shore Protection Manual (SPM), U.S. Army Corps of Engineers (1984) in deep water
        #for fully developed condition
    
        UA=0.71*windvel**1.23 #Wind stress factor or adjusted wind velcity
        tminhatdeepUAFullyDev=7.15e4 #Dimensionless minimum required wind duration: tminhat=g*tmin/U10
        tmindeepFullyDev=tminhatdeepUAFullyDev*UA/9.81 #Minimum required wind duration in (s)
        tminhatdeepFullyDev=9.81*tmindeepFullyDev/windvel #Dimensionless minimum required wind duration: tminhat=g*tmin/U10
    
        #Checking that tmin in shallow water be smaller than fully developed condition
        if len(tminhatdeepFullyDev)==1:
            tminhat[tminhat>tminhatdeepFullyDev]=tminhatdeepFullyDev
        else:
            tminhat[tminhat>tminhatdeepFullyDev]=tminhatdeepFullyDev[tminhat>tminhatdeepFullyDev]
    
        tmin=tminhat*windvel/9.81 #Minimum required wind duration in (s)
    
        #Recalculating Fetchhat to recover values that set to Fetchhat=20000
        Fetchhat=9.81*Fetch/(windvel**2) #Dimensionless Fetch: Fetchhat=g*Fetch/U10^2
    
    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':

        #Plotting
        plt.loglog(Fetchhat,tminhat)

        plt.xlabel('Dimensionless Wind Fetch, Fetchhat=g*Fetch/U10^2')
        plt.ylabel('Dimensionless minimum duration, tminhat=g*tmin/U10')

    #--------------------------------------------------------------------------
    #Outputs
    return tmin, tminhat, Fetchhat, hhat

    #--------------------------------------------------------------------------
