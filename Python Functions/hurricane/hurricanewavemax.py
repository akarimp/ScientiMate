def hurricanewavemax(xCenter, yCenter, Pc, Vmax, Rknown, VgRknown, Pn=101325, Rhoa=1.204, VgToVCoeff=0.8, G=0.88, Vt=0, CalcMethod='young13', RVmaxCalcMethod='holland', dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2017-10-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.hurricanewavemax
    ============================

    .. code:: python

        Hsmax, Tmax, RVmax = scientimate.hurricanewavemax(xCenter, yCenter, Pc, Vmax, Rknown, VgRknown, Pn=101325, Rhoa=1.204, VgToVCoeff=0.8, G=0.88, Vt=0, CalcMethod='young13', RVmaxCalcMethod='holland', dispout='no')

    Description
    -----------

    Calculate hurricane maximum wave height and wave period at a location of maximum wind 

    Inputs
    ------

    xCenter
        x (longitude) of hurricane center (track)
    yCenter
        y (latitude) of hurricane center (track)
    Pc
        Hurricane central surface pressure in (Pa)
    Vmax
        | Maximum hurricane 1-min averaged wind velocity at surface level in (m/s)
        | It can be estimated from Vmax=3.44*(1010-Pc*1e-2)^0.644; where Pc is in (Pa)
        | and, Vmax is 1-min averaged wind at a 10-m elevation in (m/s) (Atkinson & Holliday, 1977)
        | Hurricane translation velocity (forward velocity) needs to be removed from Vmax before applying (e.g. Hu et al., 2012)
        | Hurricane translation velocity (forward velocity) needs to be added after rotational velocity is calculated
    Rknown
        | Radius that hurricane wind velocity is known at that radius in (m)
        | Rknown should be larger than radius associated to Vgmax
    VgRknown
        Hurricane wind velocity at the gradient level which is known at radius Rknown in (m/s)
    Pn=101325
        | Ambient surface pressure (external pressure) in (Pa)
        | Standard atmosphere pressure is 101325 (Pa) 
        | Typical value: Pn=101500 (Pa) for the western North Pacific, Pn= 101000 (Pa) for the North Atlantic
        | (Batke et al., 2014)
    Rhoa=1.204
        Air density at the gradient level in (kg/m3)
    VgToVCoeff=0.8
        | Coefficient to convert gradient wind velocity to wind velocity at 10 m above surface as: 
        | V=VgToVCoeff*Vg, if VgToVCoeff=1, then V=Vg
    G=0.88
        | Wind gust factor to convert 1-min averaged wind to 10-min averaged wind
        | e.g. Young (2017); Liu et al. (2017)
        | G=U(600s)/U(60s), therefore U(600s)=G*U(60s)
        | G=1/1.1; based on Powell and Houston (1996)
        | G=1/1.08; based on Harper (2013)
        | G=0.88; based on World Meteorological Organization (2015)
        | G=1/1.08 to 1/1.16; based on Liu et al. (2017)
    Vt=0
        Hurricane central translational velocity in (m/s)
    CalcMethod='young13'
        | Hurricane maximum wave height and period calculation method 
        | 'spm': Use method by Shore Protection Manual (SPM),
        |     U.S. Army Corps of Engineers (1984) in deep water
        | 'young88': Use method by Young (1988)
        | 'young13': Use method by Young & Vinoth (2013)
    RVmaxCalcMethod='holland'
        | Method for calculating a hurricane radius of the maximum wind from known velocity and its radius
        | 'rankine': Velocity is calculated using modiﬁed Rankine vortex model, Depperman (1947)
        | 'holland': Velocity is calculated using using Holland (1980)
        | 'demaria': Velocity is calculated using using DeMaria (1987)
        | 'slosh': Velocity is calculated using using SLOSH model by Jelesnianski et al. (1992)
        | 'emanuel': Velocity is calculated using using Emanuel and Rotunno (2011)
    dispout='no'
        | Define to display outputs or not
        | 'imagesc': 2 dimensional plot using imagesc or imshow
        | 'pcolor': 2 dimensional plot using pcolor
        | 'contour': 2 dimensional contour plot, number of contour=ncolor
        | 'no': not display 
        | Use dispout='no' if calculation mesh is not 2d array
        | if there is more than one time step, only the last one is plotted

    Outputs
    -------

    Hsmax
        Hurricane maximum significant wave height in (m) 
    Tmax
        | Hurricane maximum wave period in (s) 
        | For all methods excepth for 'spm': Tmax=Tpmax where, Tpmax is a maximum peak wave period
        | For 'spm' method: Tmax=Tsmax, where Tsmax is a maximum significant wave period (Ts=0.95Tp)
    RVmax
        Distance (radius) from hurricane center to a location of maximum hurricane wind velocity (m)

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np


        #EXAMPLE 1

        #Longitude of Hurricane Katrine center at max velocity
        longCenter=-88.6

        #Latitude of Hurricane Katrine center at max velocity
        latCenter=26.3

        #Hurricane Katrina centeral pressure (Pa) at max velocity
        Pc=90200

        #Hurricane Katrina translational velocity (m/s) at max velocity
        Vt=5.18467

        #Hurricane Katrina 1-min sustained maximum velocity (m/s) at max velocity
        Vmax=76.5
        Vmax=Vmax-Vt #Removing hurricane translation velocity from Vgmax

        #34 kt (17.49 m/s) wind radii maximum extent in northeastern quadrant in (m) for Hurricane Katrina at max velocity
        Rknown=370400
        VRknown=17.49
        VRknown=VRknown-Vt #Removing hurricane translation velocity from VRknown
        VgRknown=VRknown/0.8 #Converting surface velocity to gradient velocity

        Pn=101325 #Ambient surface pressure (external pressure) in (Pa)
        Rhoa=1.204 #Air density in (kg/m3)

        Hsmax,Tmax,RVmax=sm.hurricanewavemax(longCenter,latCenter,Pc,Vmax,Rknown,VgRknown,Pn,Rhoa,0.8,0.88,Vt,'young13','holland','yes')


        #EXAMPLE 2

        #Longitude of Hurricane Katrine best track
        longtrack=[-75.1,-75.7,-76.2,-76.5,-76.9,-77.7,-78.4,-79.0,-79.6,-80.1,-80.3,-81.3,\
            -82.0,-82.6,-83.3,-84.0,-84.7,-85.3,-85.9,-86.7,-87.7,-88.6,-89.2,-89.6,\
            -89.6,-89.6,-89.6,-89.6,-89.1,-88.6,-88.0,-87.0,-85.3,-82.9]

        #Latitude of Hurricane Katrine best track
        lattrack=[23.1,23.4,23.8,24.5,25.4,26.0,26.1,26.2,26.2,26.0,25.9,25.4,\
            25.1,24.9,24.6,24.4,24.4,24.5,24.8,25.2,25.7,26.3,27.2,28.2,\
            29.3,29.5,30.2,31.1,32.6,34.1,35.6,37.0,38.6,40.1]

        #Hurricane Katrina centeral pressure (Pa)
        Pc=[100800,100700,100700,100600,100300,100000,99700,99400,98800,98400,98300,98700,\
            97900,96800,95900,95000,94200,94800,94100,93000,90900,90200,90500,91300,\
            92000,92300,92800,94800,96100,97800,98500,99000,99400,99600]

        #Hurricane Katrina translational velocity (m/s)
        Vt=np.array([0.00000,3.23091,3.13105,3.86928,4.99513,4.82816,3.27813,2.81998,2.77140,2.53041,\
            1.05928,5.30662,3.60661,2.98269,3.61863,3.43691,3.28168,2.85849,3.20404,4.26279,\
            5.31340,5.18467,5.39195,5.46121,5.66270,1.02958,3.60354,4.63312,8.02540,8.01558,\
            8.12721,8.31580,10.75406,12.28350])
            
        #Hurricane Katrina 1-min sustained maximum velocity (m/s)
        Vmax=np.array([15.3,15.3,15.3,17.850,20.4,22.950,25.5,28.050,30.6,35.7,35.7,33.150,\
            38.250,43.350,45.9,48.450,51.0,51.0,51.0,63.750,73.950,76.5,71.4,63.750,\
            56.1,56.1,53.550,40.8,25.5,20.4,15.3,15.3,15.3,12.750])

        Vmax=Vmax-Vt #Removing hurricane translation velocity from Vmax

        #34 kt (17.49 m/s) wind radii maximum extent in northeastern quadrant in (m) for Hurricane Katrina
        RknownRaw=[0,0,0,111120,111120,111120,111120,111120,129640,np.nan,129640,138900,\
            138900,138900,166680,240760,240760,259280,259280,296320,333360,370400,370400,370400,\
            np.nan,370400,np.nan,185200,138900,138900,0,0,0,0]

        #34 kt (17.49 m/s) wind radii maximum extent in northeastern quadrant in (m) for Hurricane Katrina
        Rknown=[0,0,0,111120,111120,111120,111120,111120,129640,129640,129640,138900,\
            138900,138900,166680,240760,240760,259280,259280,296320,333360,370400,370400,370400,\
            370400,370400,277800,185200,138900,138900,0,0,0,0]
        VRknown=np.ones(34)*17.49
        VRknown=VRknown-Vt #Removing hurricane translation velocity from VRknown
        VgRknown=VRknown/0.8 #Converting surface velocity to gradient velocity

        Pn=101325 #Ambient surface pressure (external pressure) in (Pa)
        Rhoa=1.204 #Air density in (kg/m3)

        Hsmax,Tmax,RVmax=sm.hurricanewavemax(longtrack[3:27],lattrack[3:27],Pc[3:27],Vmax[3:27],Rknown[3:27],VgRknown[3:27],Pn,Rhoa,0.8,0.88,Vt[3:27],'young13','holland','yes')


        #EXAMPLE 3

        longCenter=0 #(Degree)
        latCenter=20 #(Degree)
        Pc=90200 #(Pa)
        Vt=5.18467 #(m/s)
        Vmax=76.5 #(m/s)
        Vmax=Vmax-Vt
        Rknown=370400 #(m)
        VRknown=17.49 #(m/s)
        VRknown=VRknown-Vt 
        VgRknown=VRknown/0.8 #(m/s)
        Pn=101325 #Ambient surface pressure (external pressure) in (Pa)
        Rhoa=1.204 #Air density in (kg/m3)
        Vmax=Vmax-Vt

        Hsmax,Tmax,RVmax=sm.hurricanewavemax(longCenter,latCenter,Pc,Vmax,Rknown,VgRknown,Pn,Rhoa,0.8,0.88,Vt,'young13','holland','yes')

    References
    ----------

    Data

    * www.nhc.noaa.gov/data/
    * www.nhc.noaa.gov/data/hurdat/hurdat2-format-nencpac.pdf
    * coast.noaa.gov/hurricanes
    * www.aoml.noaa.gov/hrd/data_sub/re_anal.html

    Atkinson, G. D., & Holliday, C. R. (1977). 
    Tropical cyclone minimum sea level pressure/maximum sustained wind relationship for the western north Pacific. 
    Monthly Weather Review, 105(4), 421-427.

    Batke, S. P., Jocque, M., & Kelly, D. L. (2014). 
    Modelling hurricane exposure and wind speed on a mesoclimate scale: a case study from Cusuco NP, Honduras. 
    PloS one, 9(3), e91306.

    Department of the Army, Waterways Experiment Station, Corps of Engineers, 
    and Coastal Engineering Research Center (1984), 
    Shore Protection Manual, Washington, 
    D.C., vol. 1, 4th ed., 532 pp.

    Graham and Numm (1959) 
    Meteorological Conditions Pertinent to Standard Project Hurricane, Atlantic and Gulf Coasts of United States.
    National Hurricane Research Project. U.S. Weather Service, Report no. 33.

    Harper, B.A. (2013)
    Best practice in tropical cyclone wind hazard modelling: In search of data and emptying the skeleton cupboard. 
    In Proceedings of the 16th Australasian Wind Engineering Society Workshop, Brisbane, Qld, Australia, 18–19 July 2013

    Holland, G. J. (1980). 
    An analytic model of the wind and pressure profiles in hurricanes. 
    Monthly weather review, 108(8), 1212-1218.

    Hu, K., Chen, Q., & Kimball, S. K. (2012). 
    Consistency in hurricane surface wind forecasting: an improved parametric model. 
    Natural hazards, 61(3), 1029-1050.

    Jelesnianski, C. P., Chen, J., & Shaffer, W. A. (1992). 
    SLOSH: Sea, lake, and overland surges from hurricanes (Vol. 48). 
    US Department of Commerce, National Oceanic and Atmospheric Administration, National Weather Service.

    Liu, Q., Babanin, A., Fan, Y., Zieger, S., Guan, C., & Moon, I. J. (2017). 
    Numerical simulations of ocean surface waves under hurricane conditions: Assessment of existing model performance. 
    Ocean Modelling, 118, 73-93.

    Moon, I. J., Ginis, I., Hara, T., Tolman, H. L., Wright, C. W., & Walsh, E. J. (2003). 
    Numerical simulation of sea surface directional wave spectra under hurricane wind forcing. 
    Journal of physical oceanography, 33(8), 1680-1706.

    Phadke, A. C., Martino, C. D., Cheung, K. F., & Houston, S. H. (2003). 
    Modeling of tropical cyclone winds and waves for emergency management. 
    Ocean Engineering, 30(4), 553-578.

    Powell, M. D., & Houston, S. H. (1996). 
    Hurricane Andrew's landfall in South Florida. Part II: Surface wind fields and potential real-time applications. 
    Weather and Forecasting, 11(3), 329-349.

    Powell, M. D., Vickery, P. J., & Reinhold, T. A. (2003). 
    Reduced drag coefficient for high wind speeds in tropical cyclones. 
    Nature, 422(6929), 279.

    Valamanesh, V., Myers, A. T., Arwade, S. R., Hajjar, J. F., Hines, E., & Pang, W. (2016). 
    Wind-wave prediction equations for probabilistic offshore hurricane hazard analysis. 
    Natural Hazards, 83(1), 541-562.

    Wei, K., Arwade, S. R., Myers, A. T., Valamanesh, V., & Pang, W. (2017). 
    Effect of wind and wave directionality on the structural performance of non‐operational offshore wind turbines supported by jackets during hurricanes. 
    Wind Energy, 20(2), 289-303.

    World Meteorological Organization. Tropical Cyclone Programme, & Holland, G. J. (2015). 
    Global guide to tropical cyclone forecasting. 
    Secretariat of the World Meteorological Organization.

    Young, I. R. (1988). 
    Parametric hurricane wave prediction model. 
    Journal of Waterway, Port, Coastal, and Ocean Engineering, 114(5), 637-652.

    Young, I. R. (2006). 
    Directional spectra of hurricane wind waves. 
    Journal of Geophysical Research: Oceans, 111(C8).

    Young, I. R., & Vinoth, J. (2013). 
    An “extended fetch” model for the spatial distribution of tropical cyclone wind–waves as observed by altimeter. 
    Ocean Engineering, 70, 14-24.

    Young, I.R. (2017)
    A Review of Parametric Descriptions of Tropical Cyclone Wind-Wave Generation.
    Atmosphere 2017, 8, 194.

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
    
    xCenter=type2numpy(xCenter)
    yCenter=type2numpy(yCenter)
    Pc=type2numpy(Pc)
    Vmax=type2numpy(Vmax)
    Rknown=type2numpy(Rknown)
    VgRknown=type2numpy(VgRknown)
    Vt=type2numpy(Vt)

    #--------------------------------------------------------------------------
    #Pre-assigning array

    RVmax=np.zeros(len(xCenter)) #Pre-assigning array

    #--------------------------------------------------------------------------
    #Calculating radius of maximum wind from known velocity and its radius

    #Converting gradient (mean boundary-layer) wind velocity to velocity 10 m above surface
    #e.g. Graham & Numm (1959), Young & Vinoth (2013), Phadke et al. (2003), Powell et al. (2003), Valamanesh et al. (2016), Wei et al. (2017)
    #Shore Protection Manual (SPM), U.S. Army Corps of Engineers (1984)
    #Coastal Engineering Manual (CEM), U.S. Army Corps of Engineers (2015)
    Vgmax=Vmax/VgToVCoeff
    
    deltaP=(Pn-Pc) #Pressure drop from a defined external pressure to the cyclone center
    B=Vgmax**2*Rhoa*np.exp(1)/(deltaP) #Hurricane shape parameter
    B[B>2.5]=2.5 #B is typically between 1 and 2.5
    B[B<0.5]=1 #B is typically between 1 and 2.5
    
    w=7.2921150e-5 #Rotational frequency of the earth (Radian/s)
    f=2*w*np.sin(np.deg2rad(yCenter)) #Coriolis parameter (Coriolis frequency)
    
    for i in range(0,len(xCenter),1):
        
        #Initial guess for RVmax using Rankine vortex model
        #Rankine vortex model for R>=RVmax: V=Vgmax*(RVmax/R)^x, where 0.4<x<0.6 (e,g. Phadke et al., 2003; Holland et al., 2010)
        #Rankine vortex model for R<RVmax: V=Vgmax*(R/RVmax)^x, where 0.4<x<0.6 (e,g. Phadke et al., 2003; Holland et al., 2010)
        RVmaxini=Rknown[i]*(VgRknown[i]/Vgmax[i])**2 #Initial guess
    
        #Calculating radius of maximum wind from known velocity and its radius using modified Rankine vortex model, Depperman (1947) 
        if RVmaxCalcMethod=='rankine':
    
            #Calculating radius of maximum wind from known velocity and its radius
            RVmax[i]=RVmaxini #Radius of maximum wind
    
        #Calculating radius of maximum wind from known velocity and its radius using Holland (1980)
        elif RVmaxCalcMethod=='holland':
    
            #Calculating radius of maximum wind from known velocity and its radius
            fun = lambda x : (VgRknown[i]-(B[i]/Rhoa*(x/Rknown[i])**B[i]*deltaP[i]*np.exp(-(x/Rknown[i])**B[i]))**0.5) #Hurricane radius for Vgmax without coriolis parameter
            RVmax[i]=sp.optimize.fsolve(fun,RVmaxini) #Radius of maximum wind
    
        #Calculating radius of maximum wind from known velocity and its radius using DeMaria (1987), e.g. Holland et al. (2010)
        elif RVmaxCalcMethod=='demaria':
    
            #Calculating radius of maximum wind from known velocity and its radius
            fun = lambda x : (VgRknown[i]-Vgmax[i]*(Rknown[i]/x)*np.exp(1-(Rknown[i]/x))) #Hurricane radius for Vgmax
            RVmax[i]=sp.optimize.fsolve(fun,RVmaxini) #Radius of maximum wind
    
        #Calculating radius of maximum wind from known velocity and its radius using SLOSH model by Jelesnianski et al. (1992)
        elif RVmaxCalcMethod=='slosh':
    
            #Calculating radius of maximum wind from known velocity and its radius
            R1=float(Rknown[i])
            #fun = lambda x : (VgRknown[i]-Vgmax[i]*(2*x*Rknown[i])/((x)**2+(Rknown[i])**2)) #Hurricane radius for Vgmax
            fun = lambda x : (VgRknown[i]-Vgmax[i]*(2*x*R1)/((x)**2+(R1)**2)) #Hurricane radius for Vgmax
            RVmax[i]=sp.optimize.fsolve(fun,RVmaxini) #Radius of maximum wind
    
        #Calculating radius of maximum wind from known velocity and its radius using Emanuel and Rotunno (2011)
        elif RVmaxCalcMethod=='emanuel':
    
            #Calculating radius of maximum wind from known velocity and its radius
            R1=float(Rknown[i])
            #fun = lambda x : (VgRknown[i]-2*Rknown[i]*(x*Vgmax[i]+0.5*f[i]*(x)**2)/((x)**2+(Rknown[i])**2)-0.5*f[i]*Rknown[i]) #Hurricane radius for Vgmax
            fun = lambda x : (VgRknown[i]-2*R1*(x*Vgmax[i]+0.5*f[i]*(x)**2)/((x)**2+(R1)**2)-0.5*f[i]*R1) #Hurricane radius for Vgmax
            RVmax[i]=sp.optimize.fsolve(fun,RVmaxini) #Radius of maximum wind
    
    
    #--------------------------------------------------------------------------
    #Calculating hurricane maximum wave height and wave period at a location of maximum wind 

    #Converting 1-min sustained wind to 10-min averaged wind using gust factor
    #e.g. World Meteorological Organization (2015)
    U10=Vmax*G
    
    
    #Calculating wave properties using Shore Protection Manual (SPM), U.S. Army Corps of Engineers (1984)
    if CalcMethod=='spm':
    
        RVmaxKm=RVmax/1000 #Convert m to km
        deltaPmmHg=deltaP*0.00750061683 #Convert Pa to mmHg
        VtKt=Vt*1.94384 #Convert m/s to knots
        wradhr=2*np.pi/24 #Rotational frequency of the earth (Radian/hr)
        fradhr=2*wradhr*np.sin(np.deg2rad(yCenter)) #Coriolis parameter (Coriolis frequency) in rad/hr
        Umax=0.447*(14.5*deltaPmmHg**0.5-RVmaxKm*0.31*fradhr)
        UR=0.865*Umax+0.5*VtKt
    
        Hs=5.03*np.exp(RVmaxKm*deltaPmmHg/4700)*(1+0.29*1.0*VtKt/(np.sqrt(UR))) #Maximum significant wave height (m)
        Ts=8.6*np.exp(RVmaxKm*deltaPmmHg/9400)*(1+0.145*1.0*VtKt/(np.sqrt(UR))) #Maximum significant wave period (s)
    
        Hsmax=Hs.copy() #Maximum significant wave height (m)
        Tmax=Ts.copy() #Maximum significant wave period (m)
    
    
    #Calculating hurricane wave propoerties based on Young (1988) 
    elif CalcMethod=='young88':
    
        a=-2.175e-3
        b=1.506e-2
        c=-1.223e-1
        d=2.190e-1
        e1=6.737e-1
        f=7.980e-1
    
        Rprime=22.5e3*np.log10(RVmax)-70.8e3
    
        #Calculating an equivalent fetch Feq
        Feq=(a*Vmax**2+b*Vmax*Vt+c*Vt**2+d*Vmax+e1*Vt+f)*Rprime #Equivalent wind fetch
    
        Fetchhat=9.81*Feq/(Vmax**2) #Dimensionless equivalent fetch: Fetchhat=g*Fetch/U10^2
        
        #Calculating hurricane wave height based on JONSWAP (Hasselmann etal.,1973)
        Hshat=0.0016*Fetchhat**0.5 #Dimensionless significant wave height: Hshat=g*Hs/U10^2
        Hs=Hshat*Vmax**2/9.81 #Maximum significant wave height
    
        #Calculating hurricane wave period based on JONSWAP (Hasselmann etal.,1973)
        #fphat=3.5*Fetchhat**(-0.33) #Dimensionless peak wave frequency: fphat=fp*U10/g
        #fp=fphat*9.81/Vmax #Peak wave frequency, fphat=fp*U10/g
        #Tp=1/fp #Maximum peak wave period, Tp=1/fp
    
        #Calculating hurricane wave period based on Young (2006)
        Hm0=Hs.copy() #Zero-moment wave height, Hm0=Hs
        m0=(Hm0/4)**2 #Zero-moment of power spectral density, Hm0=4*(m0)^0.5
        Ehat=9.81**2*m0/(U10**4) #Dimensionless wave energy: Ehat=g^2*m0/U10^4
        fphat=(Ehat/6.365e-6)**(-1/3.3) #Dimensionless peak wave frequency: fphat=fp*U10/g
        fp=fphat*9.81/U10 #Peak wave frequency, fphat=fp*U10/g
        Tp=1/fp #Maximum peak wave period, Tp=1/fp
    
        Hsmax=Hs.copy() #Maximum significant wave height (m)
        Tmax=Tp.copy() #Maximum peak wave period (m)
    
    
    #Calculating hurricane wave propoerties based on Young & Vinoth (2013) 
    elif CalcMethod=='young13':
    
        a=-2.175e-3
        b=1.506e-2
        c=-1.223e-1
        d=8.760e-2
        e1=1.516
        f=1.756
        
        Rprime=22.5e3*np.log10(RVmax)-70.8e3
    
        #Calculating an equivalent fetch Feq
        Feq=(a*Vmax**2+b*Vmax*Vt+c*Vt**2+d*Vmax+e1*Vt+f)*Rprime #Equivalent wind fetch
    
        Fetchhat=9.81*Feq/(Vmax**2) #Dimensionless equivalent fetch: Fetchhat=g*Fetch/U10^2
        
        #Calculating hurricane wave height based on JONSWAP (Hasselmann etal.,1973)
        Hshat=0.0016*Fetchhat**0.5 #Dimensionless significant wave height: Hshat=g*Hs/U10^2
        Hs=Hshat*Vmax**2/9.81 #Maximum significant wave height
    
        #Calculating hurricane wave period based on JONSWAP (Hasselmann etal.,1973)
        #fphat=3.5*Fetchhat**(-0.33) #Dimensionless peak wave frequency: fphat=fp*U10/g
        #fp=fphat*9.81/Vmax #Peak wave frequency, fphat=fp*U10/g
        #Tp=1/fp #Maximum peak wave period, Tp=1/fp
    
        #Calculating hurricane wave period based on Young (2006)
        Hm0=Hs.copy() #Zero-moment wave height, Hm0=Hs
        m0=(Hm0/4)**2 #Zero-moment of power spectral density, Hm0=4*(m0)^0.5
        Ehat=9.81**2*m0/(U10**4) #Dimensionless wave energy: Ehat=g^2*m0/U10^4
        fphat=(Ehat/6.365e-6)**(-1/3.3) #Dimensionless peak wave frequency: fphat=fp*U10/g
        fp=fphat*9.81/U10 #Peak wave frequency, fphat=fp*U10/g
        Tp=1/fp #Maximum peak wave period, Tp=1/fp
    
        Hsmax=Hs.copy() #Maximum significant wave height (m)
        Tmax=Tp.copy() #Maximum peak wave period (m)
    
    
    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
    
        #Plotting
        plt.subplot(1,2,1)
        plt.scatter(RVmax,Hsmax)
    
        plt.xlabel('RVmax (m)')
        plt.ylabel('Hsmax (m)')
    
        plt.subplot(1,2,2)
        plt.scatter(RVmax,Tmax)
    
        plt.xlabel('RVmax (m)')
        plt.ylabel('Tmax (s)')
    

    #--------------------------------------------------------------------------
    #Outputs
    return Hsmax, Tmax, RVmax

    #--------------------------------------------------------------------------
