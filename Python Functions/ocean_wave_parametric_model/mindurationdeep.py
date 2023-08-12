def mindurationdeep(windvel, Fetch, CalcMethod='carter', dispout='no'):
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

    scientimate.mindurationdeep
    ===========================

    .. code:: python

        tmin, tminhat, Fetchhat = scientimate.mindurationdeep(windvel, Fetch, CalcMethod='carter', dispout='no')

    Description
    -----------

    Calculate a minimum required wind duration for wave to be fetch limited in deep water

    Inputs
    ------

    windvel
        | Wind velocity in (m/s)
        | Wind velocity should be measured (or represents velocity) at 10 m above surface
        | Wind velocity should be a 10 minutes averaged wind for all methods, except for for 'cem' and 'spmdeep' methods
        | For 'cem' and 'spmdeep' methods, wind velocity should be converted to duration of sustained wind by using gust factor
    Fetch
        Wind fetch in (m)
    CalcMethod='carter'
        | Parametric wave model to be used for calculation 
        | 'wilson': Use method by Wislon (1965)
        | 'jonswap': Use method by Hasselmann et al. (1973) known as JONSWAP
        | 'carter': Use method by Carter (1982)
        | 'spmdeep': Use method by Shore Protection Manual (SPM),
        |     U.S. Army Corps of Engineers (1984) in deep water
        | 'cem': Use method by Coastal Engineering Manual (CEM),
        |     U.S. Army Corps of Engineers (2015)
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    tmin
        Minimum required wind duration for wind to be fetch limited in (second)
    tminhat
        Dimensionless minimum required wind duration for wind to be fetch limited: tminhat=g*tmin/U10
    Fetchhat
        | Dimensionless Fetch: Fetchhat=g*Fetch/U10^2
        | Note, g=9.81: gravitational acceleration
        |     U10: wind velocity

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        windvel=10*np.random.rand(100)
        Fetch=10000*np.random.rand(100)
        tmin,tminhat,Fetchhat=sm.mindurationdeep(windvel,Fetch,'carter','no')

        windvel=10
        Fetch=np.arange(1e3,1e6,1000)
        tmin,tminhat,Fetchhat=sm.mindurationdeep(windvel,Fetch,'carter','yes')

    References
    ----------

    Carter, D. J. T. (1982). 
    Prediction of wave height and period for a constant wind velocity using the JONSWAP results. 
    Ocean Engineering, 9(1), 17-33.

    Department of the Army, Waterways Experiment Station, Corps of Engineers, 
    and Coastal Engineering Research Center (1984), 
    Shore Protection Manual, Washington, 
    D.C., vol. 1, 4th ed., 532 pp.

    Hasselmann, K.; Barnett, T. P.; Bouws, E.; Carlson, H.; Cartwright, D. E.; Enke, K.; Ewing, J. A.; 
    Gienapp, H.; Hasselmann, D. E.; Kruseman, P.; Meerbrug, A.; Muller, P.; Olbers, D. J.; Richter, K.; 
    Sell, W., and Walden, H., (1973). 
    Measurements of wind-wave growth and swell decay during the Joint North Sea Wave Project (JONSWAP). 
    Deutsche Hydrographische Zeitschrift A80(12), 95p.

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
    #Calculating a minimum required wind duration for wave to be fetch limited

    #Calculating a minimum required wind duration for wave to be fetch limited using Wislon (1965)
    if CalcMethod=='wilson':
    
        Fetchhat=9.81*Fetch/(windvel**2) #Dimensionless Fetch: Fetchhat=g*Fetch/U10^2

        tmin=43*Fetch**0.73/(windvel**0.46*9.81**0.27) #Minimum required wind duration in (s)
        tminhat=9.81*tmin/windvel #Dimensionless minimum required wind duration: tminhat=g*tmin/U10
    
    
    #Calculating a minimum required wind duration for wave to be fetch limited using JONSWAP, Hasselmann et al. (1973)
    elif CalcMethod=='jonswap':
    
        Fetchhat=9.81*Fetch/(windvel**2) #Dimensionless Fetch: Fetchhat=g*Fetch/U10^2

        tmin=65.6*Fetch**0.67/(windvel**0.33*9.81**0.33) #Minimum required wind duration in (s)
        tminhat=9.81*tmin/windvel #Dimensionless minimum required wind duration: tminhat=g*tmin/U10
    
    
    #Calculating a minimum required wind duration for wave to be fetch limited using Carter (1982)
    elif CalcMethod=='carter':
    
        Fetchhat=9.81*Fetch/(windvel**2) #Dimensionless Fetch: Fetchhat=g*Fetch/U10^2

        FetchKm=Fetch/1000 #Wind Fetch in (km)
        tminHr=1.167*FetchKm**0.7/(windvel**0.4) #Minimum required wind duration in (hr)
        tmin=tminHr*3600 #Minimum required wind duration in (second)
        tminhat=9.81*tmin/windvel #Dimensionless minimum required wind duration: tminhat=g*tmin/U10
    
    
    #Calculating a minimum required wind duration for wave to be fetch limited
    #using Shore Protection Manual (SPM), U.S. Army Corps of Engineers (1984) in deep water
    elif CalcMethod=='spmdeep':
    
        Fetchhat=9.81*Fetch/(windvel**2) #Dimensionless Fetch: Fetchhat=g*Fetch/U10^2
        UA=0.71*windvel**1.23 #Wind stress factor or adjusted wind velcity
        FetchhatUA=9.81*Fetch/(UA**2) #Dimensionless Fetch, FetchhatUA=g.Fetch/UA^2
    
        tminstarhat=6.88e1*(FetchhatUA)**(2/3) #Dimensionless minimum required wind duration: tminhat=g*tmin/UA
        tmin=tminstarhat*UA/9.81 #Minimum required wind duration
        tminhat=9.81*tmin/windvel #Dimensionless minimum required wind duration: tminhat=g*tmin/U10
    
    
    #Calculating a minimum required wind duration for wave to be fetch limited
    #using Coastal Engineering Manual (CEM), U.S. Army Corps of Engineers (2015)
    elif CalcMethod=='cem':
    
        Fetchhat=9.81*Fetch/(windvel**2) #Dimensionless Fetch: Fetchhat=g*Fetch/U10^2
        CD=0.001*(1.1+0.035*windvel) #Drag Coefficient
        Ustar=(CD*(windvel**2))**0.5 #Shear Velocity in (m/s)
        FetchhatStar=9.81*Fetch/(Ustar**2) #Dimensionless Fetch, FetchhatStar=g.Fetch/UA^2

        tmin=77.23*Fetch**0.67/(windvel**0.34*9.81**0.33) #Minimum required wind duration in (s)
        tminhat=9.81*tmin/windvel #Dimensionless minimum required wind duration: tminhat=g*tmin/U10
    
    
    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':

        #Plotting
        plt.loglog(Fetchhat,tminhat)

        plt.xlabel('Dimensionless Wind Fetch, Fetchhat=g*Fetch/U10^2')
        plt.ylabel('Dimensionless minimum duration, tminhat=g*tmin/U10')

    #--------------------------------------------------------------------------
    #Outputs
    return tmin, tminhat, Fetchhat

    #--------------------------------------------------------------------------
