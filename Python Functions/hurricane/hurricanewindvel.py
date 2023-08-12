def hurricanewindvel(xgrid, ygrid, xCenter, yCenter, Vgmax, Rknown, VgRknown, \
    Rmax=423e3, CalcMethod='slosh', VgToVCoeff=0.8, inflowdirCalcMethod='bretschneider', Vt=0, VtAzmdir=0, backwindCalcMethod='no', distCalcMethod='gc', flattendata='no', savedata='no', dispout='no'):
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

    scientimate.hurricanewindvel
    ============================

    .. code:: python

        Vxgrid, Vygrid, Vgrid, Vgxgrid, Vgygrid, Vggrid, Rgrid, RVmax, thetagrid, thetaVgrid, thetaVgridtan = scientimate.hurricanewindvel(xgrid, ygrid, xCenter, yCenter, Vgmax, Rknown, VgRknown, \
            Rmax=423e3, CalcMethod='slosh', VgToVCoeff=0.8, inflowdirCalcMethod='bretschneider', Vt=0, VtAzmdir=0, backwindCalcMethod='no', distCalcMethod='gc', flattendata='no', savedata='no', dispout='no')


    Description
    -----------

    | Generate hurricane wind velocity data on given (x,y) points
    | For Holland (1980) method use hurricanewindh80 function
    | For Holland (2008) method use hurricanewindh08 function

    Inputs
    ------

    xgrid
        | x (longitude) of points which outputs are calculated at
        | xgrid can be a single point or 1d or 2d array 
    ygrid
        | y (latitude) of points which outputs are calculated at
        | ygrid can be a single point or 1d or 2d array 
    xCenter
        x (longitude) of hurricane center (track)
    yCenter
        y (latitude) of hurricane center (track)
    Vgmax
        | Maximum hurricane 1-min averaged wind velocity at the gradient level in (m/s)
        | It can be estimated from Vmax=3.44*(1010-Pc*1e-2)^0.644; where Pc is in (Pa)
        | and, Vmax is 1-min averaged wind at a 10-m elevation in (m/s) (Atkinson & Holliday, 1977)
        | Hurricane translation velocity (forward velocity) needs to be removed from Vgmax before applying (e.g. Hu et al., 2012)
        | Hurricane translation velocity (forward velocity) needs to be added after rotational velocity is calculated
    Rknown
        | Radius that hurricane wind velocity is known at that radius in (m)
        | Rknown should be larger than radius associated to Vgmax
    VgRknown
        Hurricane wind velocity at the gradient level which is known at radius Rknown in (m/s)
    Rmax=423e3
        | Maximum radius of hurricane from hurricane center in (m)
        | Outputs for points with R>Rmax is set to zero
        | Median values of Rmax is 423 km (e.g. Chavas & Emanuel 2010; Lin & Chavas, 2012)
    CalcMethod='slosh'
        | Hurricane velocity calculation method 
        | 'rankine': Velocity is calculated using modified Rankine vortex model, Depperman (1947)
        | 'demaria': Velocity is calculated using using DeMaria (1987)
        | 'slosh': Velocity is calculated using using SLOSH model by Jelesnianski et al. (1992)
        | 'emanuel': Velocity is calculated using using Emanuel and Rotunno (2011)
        | For Holland (1980) method use hurricanewindh80 function
        | For Holland (2008) method use hurricanewindh08 function
    VgToVCoeff=0.8
        | Coefficient to convert gradient wind velocity to wind velocity at 10 m above surface as: 
        | V=VgToVCoeff*Vg, if VgToVCoeff=1, then V=Vg
    inflowdirCalcMethod='bretschneider'
        | Inflow angle calculation method 
        | 'no': Inflow angle are not calculated, thetaVgrid=thetaVgridtan
        | 'bretschneider': Inflow angle are calculated based on Bretschneider (1972)
        | 'sobey': Inflow angle are calculated based on Sobey et al. (1977)
    Vt=0
        Hurricane central translational velocity in (m/s)
    VtAzmdir=0
        | Hurricane center velocity azimuth (bearing) direction in (Degree)
        | azimuth (bearing) direction which is measured clockwise from the north:
        | 0 (degree): toward North, 90 (degree): toward East, 180 (degree): toward South, 270 (degree): toward West 
    backwindCalcMethod='no'
        | Calculation method for adding background wind velocity due to hurricane motion
        | background wind velocity is added to points with R<=Rmax, e.g. Chavas & Emanuel (2010), Lin & Chavas (2012)
        | 'no': background wind velocity is not added to hurricane velocities
        | 'constant': background wind velocity is added as constant value to hurricane velocities
        | 'slosh': Background wind velocity is calculated and added using SLOSH model by Jelesnianski et al. (1992)
        | 'lin': Background wind velocity is calculated and added based on Lin & Chavas (2012)
    distCalcMethod='gc'
        | Distance calculation method 
        | 'cart': Distances are calculated on cartesian coordinate
        | 'gc': Distances are calculated on Great Circle based on Vincenty formula, Vincenty (1975)
        | Earth radius coonsidered as mean earth radius=6371000 m
    flattendata='no'
        | Define if flat data or not
        | 'no': does not flat the results, outputs are in 3d array
        | 'yes': flat the results, outputs are in 2d array
    savedata='no'
        | Define if save data in a file or not in working folder
        | 'no': does not save, 
        | 'yes': save data as ascii 'dat' file, data are flatten regrdless of flattendata value
    dispout='no'
        | Define to display outputs or not
        | 'imagesc': 2 dimensional plot using imagesc or imshow
        | 'pcolor': 2 dimensional plot using pcolor
        | 'contour': 2 dimensional contour plot, number of contour=ncolor
        | 'quiver': 2 dimensional vector plot 
        | 'no': not display 
        | Use dispout='no' if calculation mesh is not 2d array
        | if there is more than one time step, only the last one is plotted
        | if flattendata='yes'; then dispout is set as dispout='no';

    Outputs
    -------

    Vxgrid
        | Hurricane 1-min averaged wind velocity at 10 m above surface in x (East) direction on defined mesh in (m/s)
        | Gradient wind velocity converted to wind velocity at 10 m above surface by V=VgToVCoeff*Vg
    Vygrid
        | Hurricane 1-min averaged wind velocity at 10 m above surface in y (North) direction on defined mesh in (m/s)
        | Gradient wind velocity converted to wind velocity at 10 m above surface by V=VgToVCoeff*Vg
    Vgrid
        | Resultant hurricane 1-min averaged wind velocity at 10 m above surface (Vx^2+Vy^2)^0.5 on defined mesh in (m/s)
        | Gradient wind velocity converted to wind velocity at 10 m above surface by V=VgToVCoeff*Vg
    Vgxgrid
        Hurricane 1-min averaged gradient wind velocity at the gradient level in x (East) direction on defined mesh in (m/s)
    Vgygrid
        Hurricane 1-min averaged gradient wind velocity at the gradient level in y (North) direction on defined mesh in (m/s)
    Vggrid
        Resultant hurricane 1-min averaged gradient wind velocity at the gradient level on defined mesh in (m/s)
    Rgrid
        Distance (radius) from hurricane center to each point on the grid
    RVmax
        Distance (radius) from hurricane center to a location of maximum hurricane wind velocity (m)
    thetagrid
        Angle from hurricane center to each point on the grid in (Degree)
    thetaVgrid
        | Inflow angle (trigonometric direction) of hurricane velocity at each grid point in (Degree)
        | Inflow angle: angle between the inwardly spiraling surface wind 
        |     and the circular isobars around the hurricane center (Boose et al., 2004)
    thetaVgridtan
        | Angle (trigonometric direction) of hurricane velocity at each grid point in (Degree)
        | thetaVgridtan is tangential angle respect to radius. 
        | Note: Outputs has dimension of [M,N,L] where [M,N] is size of the x-y grid and [L] is number of time steps
        |         If flattendata='yes'; then Outputs has dimension of [M*L,N]
        |     Hurricane translation velocity needs to be added after rotational velocity is calculated 
        |         (e.g. Hu et al., 2012; Lin & Chavas, 2012)
        |     Gradient wind velocity is converted to standard wind height as
        |         wind velocity at 10 m above surface by V=VgToVCoeff*Vg
        |     1-min averaged wind velocity needs to be converted to standard duration such as 
        |         10-min averaged wind by using a gust factor

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np
        import matplotlib.pyplot as plt


        #EXAMPLE 1

        #Creating calculation mesh
        xgrid,ygrid=np.meshgrid(np.linspace(-98,-68,100),np.linspace(16,44,100))

        #Longitude of Hurricane Katrine center at max velocity
        longCenter=-88.6

        #Latitude of Hurricane Katrine center at max velocity
        latCenter=26.3

        #Hurricane Katrina centeral pressure (Pa) at max velocity
        Pc=90200

        #Hurricane Katrina translational velocity (m/s) at max velocity
        Vt=5.18467

        #Hurricane Katrina velocity azimuth (bearing) in (Degree) at max velocity
        VtAzmdir=306.76219

        #Hurricane Katrina 1-min sustained maximum velocity (m/s) at max velocity
        Vmax=76.5
        Vmax=Vmax-Vt #Removing hurricane translation velocity from Vmax
        Vgmax=Vmax/0.8 #Converting surface velocity to gradient velocity

        #34 kt (17.49 m/s) wind radii maximum extent in northeastern quadrant in (m) for Hurricane Katrina at max velocity
        Rknown=370400
        VRknown=17.49
        VRknown=VRknown-Vt #Removing hurricane translation velocity from VRknown
        VgRknown=VRknown/0.8 #Converting surface velocity to gradient velocity

        Pn=101325 #Ambient surface pressure (external pressure) in (Pa)
        Rhoa=1.204 #Air density in (kg/m3)

        Vxgrid,Vygrid,Vgrid,Vgxgrid,Vgygrid,Vggrid,Rgrid,RVmax,thetagrid,thetaVgrid,thetaVgridtan=sm.hurricanewindvel(xgrid,ygrid,longCenter,latCenter,Vgmax,Rknown,VgRknown,\
            423e3,'slosh',0.8,'bretschneider',Vt,VtAzmdir,'slosh','gc','no','no','quiver')

        #Converting 1-min sustained wind to 10-min averaged wind using gust factor
        #e.g. World Meteorological Organization (2015)
        Vxgrid=Vxgrid*0.88
        Vygrid=Vygrid*0.88
        Vgrid=Vgrid*0.88


        #EXAMPLE 2

        #Creating calculation mesh
        xgrid,ygrid=np.meshgrid(np.linspace(-98,-68,100),np.linspace(16,44,100))

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
            
        #Hurricane Katrina velocity azimuth (bearing) in (Degree)
        VtAzmdir=[0.00000,298.67291,311.22135,338.70264,338.13626,309.94476,279.18860,280.65053,270.13245,\
            246.10095,240.96690,241.20181,244.79591,249.93382,244.88325,252.71384,270.14459,280.49918,\
            298.94148,299.05364,299.18896,306.76219,329.36839,340.59069,0.00000,0.00000,0.00000,\
            0.00000,15.67775,15.42254,18.00215,29.63266,39.49673,50.29744]

        #Hurricane Katrina 1-min sustained maximum velocity (m/s)
        Vmax=np.array([15.3,15.3,15.3,17.850,20.4,22.950,25.5,28.050,30.6,35.7,35.7,33.150,\
            38.250,43.350,45.9,48.450,51.0,51.0,51.0,63.750,73.950,76.5,71.4,63.750,\
            56.1,56.1,53.550,40.8,25.5,20.4,15.3,15.3,15.3,12.750])

        Vmax=Vmax-Vt #Removing hurricane translation velocity from Vmax
        Vgmax=Vmax/0.8 #Converting surface velocity to gradient velocity

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

        Vxgrid,Vygrid,Vgrid,Vgxgrid,Vgygrid,Vggrid,Rgrid,RVmax,thetagrid,thetaVgrid,thetaVgridtan=sm.hurricanewindvel(xgrid,ygrid,longtrack[3:27],lattrack[3:27],Vgmax[3:27],Rknown[3:27],VgRknown[3:27],\
            423e3,'slosh',0.8,'bretschneider',Vt,VtAzmdir,'slosh','gc','no','no','quiver')

        #Converting 1-min sustained wind to 10-min averaged wind using gust factor
        #e.g. World Meteorological Organization (2015)
        Vxgrid=Vxgrid*0.88
        Vygrid=Vygrid*0.88
        Vgrid=Vgrid*0.88


        #EXAMPLE 3

        xgrid=np.linspace(0,10,100) #(Degree)
        ygrid=np.ones(100)*20 #(Degree)
        longCenter=0 #(Degree)
        latCenter=20 #(Degree)
        Pc=90200 #(Pa)
        Vt=5.18467 #(m/s)
        VtAzmdir=306.76219 #(Degree) 
        Vmax=76.5 #(m/s)
        Vmax=Vmax-Vt
        Vgmax=Vmax/0.8 #(m/s)
        Rknown=370400 #(m)
        VRknown=17.49 #(m/s)
        VRknown=VRknown-Vt
        VgRknown=VRknown/0.8 #(m/s)
        Pn=101325 #Ambient surface pressure (external pressure) in (Pa)
        Rhoa=1.204 #Air density in (kg/m3)

        Vxgrid,Vygrid,Vgrid,Vgxgrid,Vgygrid,Vggrid,Rgrid,RVmax,thetagrid,thetaVgrid,thetaVgridtan=sm.hurricanewindvel(xgrid,ygrid,longCenter,latCenter,Vgmax,Rknown,VgRknown,\
        423e3,'slosh',0.8,'bretschneider',Vt,VtAzmdir,'slosh','gc','no','no','no')
        plt.plot(Rgrid,Vgrid)

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

    Boose, E. R., Serrano, M. I., & Foster, D. R. (2004). 
    Landscape and regional impacts of hurricanes in Puerto Rico. 
    Ecological Monographs, 74(2), 335-352.

    Bretschneider, C. L. (1972, January). 
    A non-dimensional stationary hurricane wave model. 
    In Offshore Technology Conference. Offshore Technology Conference.

    Chavas, D. R., & Emanuel, K. A. (2010). 
    A QuikSCAT climatology of tropical cyclone size. 
    Geophysical Research Letters, 37(18).

    DeMaria, M. (1987). 
    Tropical cyclone track prediction with a barotropic spectral model. 
    Monthly weather review, 115(10), 2346-2357.

    Department of the Army, Waterways Experiment Station, Corps of Engineers, 
    and Coastal Engineering Research Center (1984), 
    Shore Protection Manual, Washington, 
    D.C., vol. 1, 4th ed., 532 pp.

    Depperman, C. E. (1947). 
    Notes on the origin and structure of Philippine typhoons. 
    Bull. Amer. Meteor. Soc, 28(9), 399-404.

    Emanuel, K., & Rotunno, R. (2011). 
    Self-stratification of tropical cyclone outflow. Part I: Implications for storm structure. 
    Journal of the Atmospheric Sciences, 68(10), 2236-2249.

    Graham and Numm (1959) 
    Meteorological Conditions Pertinent to Standard Project Hurricane, Atlantic and Gulf Coasts of United States.
    National Hurricane Research Project. U.S. Weather Service, Report no. 33.

    Holland, G. J. (1980). 
    An analytic model of the wind and pressure profiles in hurricanes. 
    Monthly weather review, 108(8), 1212-1218.

    Holland, G. (2008). 
    A revised hurricane pressure–wind model. 
    Monthly Weather Review, 136(9), 3432-3445.

    Holland, G. J., Belanger, J. I., & Fritz, A. (2010). 
    A revised model for radial profiles of hurricane winds. 
    Monthly Weather Review, 138(12), 4393-4401.

    Hu, K., Chen, Q., & Kimball, S. K. (2012). 
    Consistency in hurricane surface wind forecasting: an improved parametric model. 
    Natural hazards, 61(3), 1029-1050.

    Jelesnianski, C. P., Chen, J., & Shaffer, W. A. (1992). 
    SLOSH: Sea, lake, and overland surges from hurricanes (Vol. 48). 
    US Department of Commerce, National Oceanic and Atmospheric Administration, National Weather Service.

    Lin, N., & Chavas, D. (2012). 
    On hurricane parametric wind and applications in storm surge modeling. 
    Journal of Geophysical Research: Atmospheres, 117(D9).

    Phadke, A. C., Martino, C. D., Cheung, K. F., & Houston, S. H. (2003). 
    Modeling of tropical cyclone winds and waves for emergency management. 
    Ocean Engineering, 30(4), 553-578.

    Powell, M. D., Vickery, P. J., & Reinhold, T. A. (2003). 
    Reduced drag coefficient for high wind speeds in tropical cyclones. 
    Nature, 422(6929), 279.

    Sobey, R. J., Harper, B. A., & Stark, K. P. (1977). 
    Numerical simulation of tropical cyclone storm surge. 
    James Cook University of North Queensland, Department of Civil & Systems Engineering.

    U.S. Army Corps of Engineers (2015). 
    Coastal Engineering Manual. 
    Engineer Manual 1110-2-1100, Washington, D.C.: U.S. Army Corps of Engineers.

    Valamanesh, V., Myers, A. T., Arwade, S. R., Hajjar, J. F., Hines, E., & Pang, W. (2016). 
    Wind-wave prediction equations for probabilistic offshore hurricane hazard analysis. 
    Natural Hazards, 83(1), 541-562.

    Wei, K., Arwade, S. R., Myers, A. T., Valamanesh, V., & Pang, W. (2017). 
    Effect of wind and wave directionality on the structural performance of non‐operational offshore wind turbines supported by jackets during hurricanes. 
    Wind Energy, 20(2), 289-303.

    World Meteorological Organization. Tropical Cyclone Programme, & Holland, G. J. (2015). 
    Global guide to tropical cyclone forecasting. 
    Secretariat of the World Meteorological Organization.

    Young, I. R., & Vinoth, J. (2013). 
    An 'extended fetch' model for the spatial distribution of tropical cyclone wind–waves as observed by altimeter. 
    Ocean Engineering, 70, 14-24.

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
    if dispout!='no':
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
    Vgmax=type2numpy(Vgmax)
    Rknown=type2numpy(Rknown)
    VgRknown=type2numpy(VgRknown)
    Vt=type2numpy(Vt)
    VtAzmdir=type2numpy(VtAzmdir)

    #--------------------------------------------------------------------------
    #Pre-assigning array

    #xgrid and ygrid are 2d array
    if np.ndim(xgrid)==2: 
        Mgrid,Ngrid=np.shape(xgrid)

    #xgrid and ygrid are 1d array
    else: 
        Mgrid=(np.shape(xgrid))[0]
        Ngrid=1 #Ngrid=1 when xgrid and ygrid are 1d array

    if len(xCenter)>1:
        Rgrid=np.zeros((Mgrid,Ngrid,len(xCenter))) #Pre-assigning array
        thetagrid=np.zeros((Mgrid,Ngrid,len(xCenter))) #Pre-assigning array
        Vgrid=np.zeros((Mgrid,Ngrid,len(xCenter))) #Pre-assigning array
        Vggrid=np.zeros((Mgrid,Ngrid,len(xCenter))) #Pre-assigning array
        thetaVgridtan=np.zeros((Mgrid,Ngrid,len(xCenter))) #Pre-assigning array
        thetaVgrid=np.zeros((Mgrid,Ngrid,len(xCenter))) #Pre-assigning array
        Beta=np.zeros((Mgrid,Ngrid,len(xCenter))) #Pre-assigning array
    else:
        Rgrid=np.zeros((Mgrid,Ngrid,1)) #Pre-assigning array
        thetagrid=np.zeros((Mgrid,Ngrid,1)) #Pre-assigning array
        Vgrid=np.zeros((Mgrid,Ngrid,1)) #Pre-assigning array
        Vggrid=np.zeros((Mgrid,Ngrid,1)) #Pre-assigning array
        thetaVgridtan=np.zeros((Mgrid,Ngrid,1)) #Pre-assigning array
        thetaVgrid=np.zeros((Mgrid,Ngrid,1)) #Pre-assigning array
        Beta=np.zeros((Mgrid,Ngrid,1)) #Pre-assigning array

    
    RVmax=np.zeros(len(xCenter)) #Pre-assigning array
    VtR=np.zeros((Mgrid,Ngrid)) #Pre-assigning array

    #--------------------------------------------------------------------------
    #Calculating distance (radius) from hurricane center to each point

    for i in range(0,len(xCenter),1):
    
        #Calculating distance using cartesian formula
        if distCalcMethod=='cart':
    
            distxy=np.sqrt((xgrid-xCenter[i])**2+(ygrid-yCenter[i])**2) #Calculating distance from (x,y) to (x(1),y(1))
    
            #Calculating angle of the line between start and end points
    
            thetarad=np.arctan2(ygrid-yCenter[i],xgrid-xCenter[i]) #Angle in radian
            theta=np.rad2deg(thetarad) #Angle in degree
    
            #Add 360 to all numbers to have them all positive
            #Use mod(360) to take care of the ones larger than 360
            theta=((theta+360)%360) 
    
            #xgrid and ygrid are 2d array
            if np.ndim(xgrid)==2:
                Rgrid[:,:,i]=distxy.copy() #Distance from hurricane center to each point in (m)
                thetagrid[:,:,i]=theta.copy() #Angle from hurricane center to each point in (degree)

            #xgrid and ygrid are 1d array
            else:
                #Using np.newaxis (or None) to convert shape (Mgrid) into shape (Mgrid,1) when xgrid and ygrid are 1d array
                Rgrid[:,:,i]=distxy[:,np.newaxis].copy() #Distance from hurricane center to each point in (m)
                thetagrid[:,:,i]=theta[:,np.newaxis].copy() #Angle from hurricane center to each point in (degree)
    
        #Calculating distance using Vincenty formula
        elif distCalcMethod=='gc':
    
            #Converting to radian
            lat1rad=np.deg2rad(yCenter[i])
            lon1rad=np.deg2rad(xCenter[i])
            lat2rad=np.deg2rad(ygrid)
            lon2rad=np.deg2rad(xgrid)
    
            deltalatrad21=lat2rad-lat1rad
            deltalonrad21=lon2rad-lon1rad
    
            REarth=6371000 #Earth radius in (m), mean earth radius=6371000 m
            deltasigma=np.arctan2(np.sqrt((np.cos(lat2rad)*np.sin(deltalonrad21))**2+(np.cos(lat1rad)*np.sin(lat2rad)-np.sin(lat1rad)*np.cos(lat2rad)*np.cos(deltalonrad21))**2),np.sin(lat1rad)*np.sin(lat2rad)+np.cos(lat1rad)*np.cos(lat2rad)*np.cos(deltalonrad21)) #Central angle
            arclen=REarth*deltasigma #Total distance of the line 
            distxy=arclen.copy()
    
            #Calculating azimuth (bearing) between start and end of the line
    
            azimuthrad=np.arctan2(np.sin(deltalonrad21)*np.cos(lat2rad),np.cos(lat1rad)*np.sin(lat2rad)-np.sin(lat1rad)*np.cos(lat2rad)*np.cos(deltalonrad21)) #Azimuth (bearing) in radian
            azimuth=np.rad2deg(azimuthrad) #Azimuth (bearing) in degree

            #Add 360 to all numbers to have them all positive
            #Use mod(360) to take care of the ones larger than 360
            azimuth=((azimuth+360)%360) 
    
            #xgrid and ygrid are 2d array
            if np.ndim(xgrid)==2:
                Rgrid[:,:,i]=distxy.copy() #Distance from hurricane center to each point in (m)
                thetagrid[:,:,i]=azimuth.copy() #Angle from hurricane center to each point in (degree)

            #xgrid and ygrid are 1d array
            else:
                #Using np.newaxis (or None) to convert shape (Mgrid) into shape (Mgrid,1) when xgrid and ygrid are 1d array
                Rgrid[:,:,i]=distxy[:,np.newaxis].copy() #Distance from hurricane center to each point in (m)
                thetagrid[:,:,i]=azimuth[:,np.newaxis].copy() #Angle from hurricane center to each point in (degree)
    

    #--------------------------------------------------------------------------
    #Calculating hurricane wind velocity

    w=7.2921150e-5 #Rotational frequency of the earth (Radian/s)
    f=2*w*np.sin(np.deg2rad(yCenter)) #Coriolis parameter (Coriolis frequency)
    
    for i in range(0,len(xCenter),1):
        
        #Initial guess for RVmax using Rankine vortex model
        #Rankine vortex model for R>=RVmax: V=Vgmax*(RVmax/R)^x, where 0.4<x<0.6 (e,g. Phadke et al., 2003 Holland et al., 2010)
        #Rankine vortex model for R<RVmax: V=Vgmax*(R/RVmax)^x, where 0.4<x<0.6 (e,g. Phadke et al., 2003 Holland et al., 2010)
        RVmaxini=Rknown[i]*(VgRknown[i]/Vgmax[i])**2 #Initial guess
    
        #Calculating hurricane wind velocity using modified Rankine vortex model, Depperman (1947) 
        if CalcMethod=='rankine':
    
            #Calculating radius of maximum wind from known velocity and its radius
            RVmax[i]=RVmaxini #Radius of maximum wind
    
            R1=Rgrid[:,:,i].copy()
            VgR=np.zeros((Mgrid,Ngrid)) #Pre-assigning array
    
            #Rankine vortex model for R<RVmax: V=Vgmax*(R/RVmax)^x, where 0.4<x<0.6 (e,g. Phadke et al., 2003; Holland et al., 2010)
            VgR[R1<RVmax[i]]=Vgmax[i]*(R1[R1<RVmax[i]]/RVmax[i])**1 #Gradient wind velocity at radius R
    
            #Rankine vortex model for R>=RVmax: V=Vgmax*(RVmax/R)^x, where 0.4<x<0.6 (e,g. Phadke et al., 2003; Holland et al., 2010)
            VgR[R1>=RVmax[i]]=Vgmax[i]*(RVmax[i]/R1[R1>=RVmax[i]])**0.5 #Gradient wind velocity at radius R
        
            VgR[Rgrid[:,:,i]>=Rmax]=0 #e.g. Chavas & Emanuel (2010), Lin & Chavas (2012)
            VgR[VgR<0]=0
    
        #Calculating hurricane wind velocity using DeMaria (1987), e.g. Holland et al. (2010)
        elif CalcMethod=='demaria':
    
            #Calculating radius of maximum wind from known velocity and its radius
            fun = lambda x : (VgRknown[i]-Vgmax[i]*(Rknown[i]/x)*np.exp(1-(Rknown[i]/x))) #Hurricane radius for Vgmax
            RVmax[i]=sp.optimize.fsolve(fun,RVmaxini) #Radius of maximum wind
    
            VgR=Vgmax[i]*(Rgrid[:,:,i]/RVmax[i])*np.exp(1-(Rgrid[:,:,i]/RVmax[i])) #Gradient wind velocity at radius R
            VgR[Rgrid[:,:,i]>=Rmax]=0 #e.g. Chavas & Emanuel (2010), Lin & Chavas (2012)
            VgR[VgR<0]=0
    
         
        #Calculating hurricane wind velocity using SLOSH model by Jelesnianski et al. (1992)
        elif CalcMethod=='slosh':
    
            #Calculating radius of maximum wind from known velocity and its radius
            R1=float(Rknown[i])
            #fun = lambda x : (VgRknown[i]-Vgmax[i]*(2*x*Rknown[i])/((x)**2+(Rknown[i])**2)) #Hurricane radius for Vgmax
            fun = lambda x : (VgRknown[i]-Vgmax[i]*(2*x*R1)/((x)**2+(R1)**2)) #Hurricane radius for Vgmax
            RVmax[i]=sp.optimize.fsolve(fun,RVmaxini) #Radius of maximum wind
    
            VgR=Vgmax[i]*(2*RVmax[i]*Rgrid[:,:,i])/((RVmax[i])**2+(Rgrid[:,:,i])**2) #Gradient wind velocity at radius R
            VgR[Rgrid[:,:,i]>=Rmax]=0 #e.g. Chavas & Emanuel (2010), Lin & Chavas (2012)
            VgR[VgR<0]=0
    
    
        #Calculating hurricane wind velocity using Emanuel and Rotunno (2011)
        elif CalcMethod=='emanuel':
    
            #Calculating radius of maximum wind from known velocity and its radius
            R1=float(Rknown[i])
            #fun = lambda x : (VgRknown[i]-2*Rknown[i]*(x*Vgmax[i]+0.5*f[i]*(x)**2)/((x)**2+(Rknown[i])**2)-0.5*f[i]*Rknown[i]) #Hurricane radius for Vgmax
            fun = lambda x : (VgRknown[i]-2*R1*(x*Vgmax[i]+0.5*f[i]*(x)**2)/((x)**2+(R1)**2)-0.5*f[i]*R1) #Hurricane radius for Vgmax
            RVmax[i]=sp.optimize.fsolve(fun,RVmaxini) #Radius of maximum wind
    
            VgR=2*Rgrid[:,:,i]*(RVmax[i]*Vgmax[i]+0.5*f[i]*(RVmax[i])**2)/((RVmax[i])**2+(Rgrid[:,:,i])**2)-0.5*f[i]*Rgrid[:,:,i] #Gradient wind velocity at radius R
            VgR[Rgrid[:,:,i]>=Rmax]=0 #e.g. Chavas & Emanuel (2010), Lin & Chavas (2012)
            VgR[VgR<0]=0
    
         
        #Converting gradient (mean boundary-layer) wind velocity to velocity 10 m above surface
        #e.g. Graham & Numm (1959), Young & Vinoth (2013), Phadke et al. (2003), Powell et al. (2003), Valamanesh et al. (2016), Wei et al. (2017)
        #Shore Protection Manual (SPM), U.S. Army Corps of Engineers (1984)
        #Coastal Engineering Manual (CEM), U.S. Army Corps of Engineers (2015)
        VR=VgToVCoeff*VgR
    
        Vggrid[:,:,i]=VgR.copy() #Assigning hurricane gradient velocity at each track point
        Vgrid[:,:,i]=VR.copy() #Assigning hurricane surface velocity at each track point
    

    #--------------------------------------------------------------------------
    #Calculating hurricane velocity tangential and inflow angle and inflow velocity in x (East) and y (North) directions
    #Tangential velocity has a right angle respect to radius
    #thetaVgridtan is tangential angle respect to radius 
    #thetaVgrid is inflow angle calculated based on Bretschneider (1972)
    #Inflow angle: angle between the inwardly spiraling surface wind 
    #and the circular isobars around the hurricane center (Boose et al., 2004)

    for i in range(0,len(xCenter),1):
    
        #Calculating inflow angle based on Bretschneider (1972), e.g Phadke et al.,(2003), 
        if inflowdirCalcMethod=='bretschneider':
            #Jeong, C. K., Panchang, V., & Demirbilek, Z. (2012). 
            #Parametric adjustments to the rankine vortex wind model for Gulf of Mexico hurricanes. 
            #Journal of Offshore Mechanics and Arctic Engineering, 134(4), 041102.
            
            #Beta is an inflow angle inward from the tangential flow (Martino et al., 2001)
            Beta[Rgrid<RVmax[i]]=10+(1+Rgrid[Rgrid<RVmax[i]]/RVmax[i])
            Beta[((Rgrid>=RVmax[i]) & (Rgrid<1.2*RVmax[i]))]=20+25*(Rgrid[((Rgrid>=RVmax[i]) & (Rgrid<1.2*RVmax[i]))]/RVmax[i]-1)
            Beta[Rgrid>=1.2*RVmax[i]]=25
    
        #Calculating inflow angle based on Sobey et al. (1977), e.g. Martino et al. (2001)
        elif inflowdirCalcMethod=='sobey':
            #Sobey, R. J., Harper, B. A., & Stark, K. P. (1977). 
            #Numerical simulation of tropical cyclone storm surge. 
            #James Cook University of North Queensland, Department of Civil & Systems Engineering.
    
            #Martino, C. D., Cheung, K. F., Phadke, A. C., & Houston, S. H. (2001, January). 
            #Modeling of hurricane waves in Hawaiian waters. 
            #In The Eleventh International Offshore and Polar Engineering Conference. International Society of Offshore and Polar Engineers.    
            
            #Beta is an inflow angle inward from the tangential flow (Martino et al., 2001)
            Beta[Rgrid<RVmax[i]]=10*(Rgrid[Rgrid<RVmax[i]]/RVmax[i])
            Beta[((Rgrid>=RVmax[i]) & (Rgrid<1.2*RVmax[i]))]=10+75*(Rgrid[((Rgrid>=RVmax[i]) & (Rgrid<1.2*RVmax[i]))]/RVmax[i]-1)
            Beta[Rgrid>=1.2*RVmax[i]]=25
    
        else:
            Beta=0
        
    
    #Calculating velocity vector direction using cartesian formula
    if distCalcMethod=='cart':
    
        #Hurricane is in northern hemisphere
        if np.mean(yCenter)>=0:
            #In northern hemisphere, velocity vector rotate counter clockwise with right angle respect to radius
            thetaVgridtan=thetagrid+90 #Angle of velocity vector in degree
            thetaVgrid=thetaVgridtan+Beta #Inflow angle of velocity vector in degree
    
        #Hurricane is in southern hemisphere
        elif np.mean(yCenter)<0:
            #In southern hemisphere, velocity vector rotate clockwise with right angle respect to radius
            thetaVgridtan=thetagrid-90 #Angle of velocity vector in degree
            thetaVgrid=thetaVgridtan-Beta #Inflow angle of velocity vector in degree
    
    
        #Add 360 to all numbers to have them all positive
        #Use mod(360) to take care of the ones larger than 360
        thetaVgridtan=((thetaVgridtan+360)%360) 
        thetaVgrid=((thetaVgrid+360)%360) 
    
    #Calculating velocity vector direction using great circle formula
    elif distCalcMethod=='gc':
    
        #Converting azimuth (bearing) to trigonometric direction
        thetagrid=-thetagrid+90
    
        #Hurricane is in northern hemisphere
        if np.mean(yCenter)>=0:
            #In northern hemisphere, velocity vector rotate counter clockwise with right angle respect to radius
            thetaVgridtan=thetagrid+90 #Angle of velocity vector in degree
            thetaVgrid=thetaVgridtan+Beta #Inflow angle of velocity vector in degree
    
        #Hurricane is in southern hemisphere
        elif mean(yCenter)<0:
            #In southern hemisphere, velocity vector rotate clockwise with right angle respect to radius
            thetaVgridtan=thetagrid-90 #Angle of velocity vector in degree
            thetaVgrid=thetaVgridtan-Beta #Inflow angle of velocity vector in degree
    
    
        #Add 360 to all numbers to have them all positive
        #Use mod(360) to take care of the ones larger than 360
        thetaVgridtan=((thetaVgridtan+360)%360) 
        thetaVgrid=((thetaVgrid+360)%360) 
    
    
    Vgxgrid=Vggrid*np.cos(np.deg2rad(thetaVgrid)) #Hurricane gradient velocity in x (East) direction
    Vgygrid=Vggrid*np.sin(np.deg2rad(thetaVgrid)) #Hurricane gradient velocity in y (North) direction
    
    Vxgrid=Vgrid*np.cos(np.deg2rad(thetaVgrid)) #Hurricane velocity in x (East) direction
    Vygrid=Vgrid*np.sin(np.deg2rad(thetaVgrid)) #Hurricane velocity in y (North) direction

    #--------------------------------------------------------------------------
    #Adding background wind velocity due to hurricane motion to hurricane wind velocity

    #Converting azimuth (bearing) to trigonometric direction
    VtTrigdir=-VtAzmdir+90
    
    #Add 360 to all numbers to have them all positive
    #Use mod(360) to take care of the ones larger than 360
    VtTrigdir=((VtTrigdir+360)%360) 
    
    #Adding background wind velocity as constant value to all hurricane velocities
    if backwindCalcMethod=='constant':
    
        for i in range(0,len(xCenter),1):
            VtR[:]=Vt[i] #Translational wind velocity at radius R
            VtR[Rgrid[:,:,i]>=Rmax]=0 #e.g. Chavas & Emanuel (2010), Lin & Chavas (2012)
    
            Vgxgrid[:,:,i]=Vgxgrid[:,:,i]+VtR*np.cos(np.deg2rad(VtTrigdir[i])) #Hurricane gradient velocity in x (East) direction
            Vgygrid[:,:,i]=Vgygrid[:,:,i]+VtR*np.sin(np.deg2rad(VtTrigdir[i])) #Hurricane gradient velocity in y (North) direction
    
            Vxgrid[:,:,i]=Vxgrid[:,:,i]+VtR*np.cos(np.deg2rad(VtTrigdir[i])) #Hurricane velocity in x (East) direction
            Vygrid[:,:,i]=Vygrid[:,:,i]+VtR*np.sin(np.deg2rad(VtTrigdir[i])) #Hurricane velocity in y (North) direction
    

        Vggrid=np.sqrt(Vgxgrid**2+Vgygrid**2) #Assigning hurricane gradient velocity at each track point
        Vgrid=np.sqrt(Vxgrid**2+Vygrid**2) #Assigning hurricane velocity at each track point
    
    #Adding background wind velocity using SLOSH model by Jelesnianski et al. (1992)
    elif backwindCalcMethod=='slosh':
    
        for i in range(0,len(xCenter),1):
            VtR=Vt[i]*(RVmax[i]*Rgrid[:,:,i])/((RVmax[i])**2+(Rgrid[:,:,i])**2) #Translational wind velocity at radius R
            VtR[VtR<0]=0
            VtR[Rgrid[:,:,i]>=Rmax]=0 #e.g. Chavas & Emanuel (2010), Lin & Chavas (2012)
    
            Vgxgrid[:,:,i]=Vgxgrid[:,:,i]+VtR*np.cos(np.deg2rad(VtTrigdir[i])) #Hurricane gradient velocity in x (East) direction
            Vgygrid[:,:,i]=Vgygrid[:,:,i]+VtR*np.sin(np.deg2rad(VtTrigdir[i])) #Hurricane gradient velocity in y (North) direction
    
            Vxgrid[:,:,i]=Vxgrid[:,:,i]+VtR*np.cos(np.deg2rad(VtTrigdir[i])) #Hurricane velocity in x (East) direction
            Vygrid[:,:,i]=Vygrid[:,:,i]+VtR*np.sin(np.deg2rad(VtTrigdir[i])) #Hurricane velocity in y (North) direction

    
        Vggrid=np.sqrt(Vgxgrid**2+Vgygrid**2) #Assigning hurricane gradient velocity at each track point
        Vgrid=np.sqrt(Vxgrid**2+Vygrid**2) #Assigning hurricane velocity at each track point
    
    #Adding background wind velocity using Lin & Chavas (2012)
    elif backwindCalcMethod=='lin':
    
        for i in range(0,len(xCenter),1):
            #Hurricane is in northern hemisphere
            if np.mean(yCenter)>=0:
                #In northern hemisphere, velocity vector rotate counter clockwise with right angle respect to radius
    
                VtR[:]=0.55*Vt[i] #Translational wind velocity at radius R
                VtR[Rgrid[:,:,i]>=Rmax]=0 #e.g. Chavas & Emanuel (2010), Lin & Chavas (2012)
    
                Vgxgrid[:,:,i]=Vgxgrid[:,:,i]+VtR*np.cos(np.deg2rad(VtTrigdir[i]+20)) #Hurricane gradient velocity in x (East) direction
                Vgygrid[:,:,i]=Vgygrid[:,:,i]+VtR*np.sin(np.deg2rad(VtTrigdir[i]+20)) #Hurricane gradient velocity in y (North) direction
    
                Vxgrid[:,:,i]=Vxgrid[:,:,i]+VtR*np.cos(np.deg2rad(VtTrigdir[i]+20)) #Hurricane velocity in x (East) direction
                Vygrid[:,:,i]=Vygrid[:,:,i]+VtR*np.sin(np.deg2rad(VtTrigdir[i]+20)) #Hurricane velocity in y (North) direction
    
            #Hurricane is in southern hemisphere
            elif np.mean(yCenter)<0:
                #In southern hemisphere, velocity vector rotate clockwise with right angle respect to radius
    
                VtR[:]=0.55*Vt[i] #Translational wind velocity at radius R
                VtR[Rgrid[:,:,i]>=Rmax]=0 #e.g. Chavas & Emanuel (2010), Lin & Chavas (2012)
    
                Vgxgrid[:,:,i]=Vgxgrid[:,:,i]+VtR*np.cos(np.deg2rad(VtTrigdir[i]-20)) #Hurricane gradient velocity in x (East) direction
                Vgygrid[:,:,i]=Vgygrid[:,:,i]+VtR*np.sin(np.deg2rad(VtTrigdir[i]-20)) #Hurricane gradient velocity in y (North) direction
    
                Vxgrid[:,:,i]=Vxgrid[:,:,i]+VtR*np.cos(np.deg2rad(VtTrigdir[i]-20)) #Hurricane velocity in x (East) direction
                Vygrid[:,:,i]=Vygrid[:,:,i]+VtR*np.sin(np.deg2rad(VtTrigdir[i]-20)) #Hurricane velocity in y (North) direction
    

        Vggrid=np.sqrt(Vgxgrid**2+Vgygrid**2) #Assigning hurricane gradient velocity at each track point
        Vgrid=np.sqrt(Vxgrid**2+Vygrid**2) #Assigning hurricane velocity at each track point
    

    #--------------------------------------------------------------------------
    #Flating results

    if flattendata=='yes':
    
        #Vgrid is 3d array
        if np.ndim(Vgrid)==3: 
            M,N,L=np.shape(Vgrid)
    
        #Vgrid is 2d array
        elif np.ndim(Vgrid)==2: 
            M,N=np.shape(Vgrid)
            L=1 #L=1 when Vgrid is 2d array
    
        #Vgrid is 1d array
        else: 
            M=(np.shape(Vgrid))[0]
            N=1 #N=1 when Vgrid is 1d array
            L=1 #L=1 when Vgrid is 1d array

        Vxgrid=np.reshape(Vxgrid,(M*L,N))
        Vygrid=np.reshape(Vygrid,(M*L,N))
        Vgrid=np.reshape(Vgrid,(M*L,N))
    
        dispout='no' #Results of flatten array cannot be plotted
    

    #--------------------------------------------------------------------------
    #Saving data

    if savedata=='yes':
    
        if flattendata=='yes':
    
            np.savetxt('Vx.dat',Vxgrid,delimiter=' ')
            np.savetxt('Vy.dat',Vygrid,delimiter=' ')
            np.savetxt('V.dat',Vgrid,delimiter=' ')
    
        #Flating results if they were not flatten before
        elif flattendata=='no':
    
            #Vgrid is 3d array
            if np.ndim(Vgrid)==3: 
                M,N,L=np.shape(Vgrid)
        
            #Vgrid is 2d array
            elif np.ndim(Vgrid)==2: 
                M,N=np.shape(Vgrid)
                L=1 #L=1 when Vgrid is 2d array
        
            #Vgrid is 1d array
            else: 
                M=(np.shape(Vgrid))[0]
                N=1 #N=1 when Vgrid is 1d array
                L=1 #L=1 when Vgrid is 1d array

            Vxgridrshp=np.reshape(Vxgrid,(M*L,N))
            Vygridrshp=np.reshape(Vygrid,(M*L,N))
            Vgridrshp=np.reshape(Vgrid,(M*L,N))
     
            np.savetxt('Vx.dat',Vxgridrshp,delimiter=' ')
            np.savetxt('Vy.dat',Vygridrshp,delimiter=' ')
            np.savetxt('V.dat',Vgridrshp,delimiter=' ')
    

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='imagesc': #2 dimensional plot using imagesc or imshow
        
        plt.imshow(Vgrid[:,:,-1],extent=[np.nanmin(xgrid),np.nanmax(xgrid),np.nanmin(ygrid),np.nanmax(ygrid)],interpolation=None,cmap=plt.get_cmap(),aspect='auto',origin='lower')
    
    
    elif dispout=='pcolor': #2 dimensional plot using pcolor
    
        plt.pcolormesh(xgrid,ygrid,Vgrid[:,:,-1],cmap=plt.get_cmap())    
    

    elif dispout=='contour': #2 dimensional contour plot
    
        plt.contourf(xgrid,ygrid,Vgrid[:,:,-1],cmap=plt.get_cmap())
    
    
    elif dispout=='quiver': #Surface plot
        
        plt.contour(xgrid,ygrid,Vgrid[:,:,-1],cmap=plt.get_cmap())
        plt.colorbar()
        plt.quiver(xgrid,ygrid,Vxgrid[:,:,-1],Vygrid[:,:,-1])
    
    
    #Setting plot properties
    if ((dispout=='imagesc') or (dispout=='pcolor') or (dispout=='contour') or (dispout=='quiver')):
    
        #Setting axis limits
        #plt.xlim(np.nanmin(xgrid),np.nanmax(xgrid))
        #plt.ylim(np.nanmin(ygrid),np.nanmax(ygrid))
    
        #Plotting colorbar
        if dispout!='quiver':
            plt.colorbar()
    
        #Setting label
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Velocity (m/s)')
        

    #--------------------------------------------------------------------------
    #Converting output array with size of (Mgrid,Ngrid,1) to array with size of (Mgrid,Ngrid)

    if ((len(xCenter)==1) and (np.ndim(xgrid)>1)):
        #xgrid and ygrid are 2d array
        #Converting array with size of (Mgrid,Ngrid,1) to array with size of (Mgrid,Ngrid)
        Vxgrid=np.reshape(Vxgrid,(Mgrid,Ngrid)) #Reshaping array
        Vygrid=np.reshape(Vygrid,(Mgrid,Ngrid)) #Reshaping array
        Vgrid=np.reshape(Vgrid,(Mgrid,Ngrid)) #Reshaping array
        Vgxgrid=np.reshape(Vgxgrid,(Mgrid,Ngrid)) #Reshaping array
        Vgygrid=np.reshape(Vgygrid,(Mgrid,Ngrid)) #Reshaping array
        Vggrid=np.reshape(Vggrid,(Mgrid,Ngrid)) #Reshaping array
        Rgrid=np.reshape(Rgrid,(Mgrid,Ngrid)) #Reshaping array
        thetagrid=np.reshape(thetagrid,(Mgrid,Ngrid)) #Reshaping array
        thetaVgrid=np.reshape(thetaVgrid,(Mgrid,Ngrid)) #Reshaping array
        thetaVgridtan=np.reshape(thetaVgridtan,(Mgrid,Ngrid)) #Reshaping array

    elif ((len(xCenter)==1) and (np.ndim(xgrid)==1)):
        #xgrid and ygrid are 1d array
        #Converting output array with size of (Mgrid,1,1) to array with size of (Mgrid,)
        Vxgrid=np.reshape(Vxgrid,(Mgrid)) #Reshaping array
        Vygrid=np.reshape(Vygrid,(Mgrid)) #Reshaping array
        Vgrid=np.reshape(Vgrid,(Mgrid)) #Reshaping array
        Vgxgrid=np.reshape(Vgxgrid,(Mgrid)) #Reshaping array
        Vgygrid=np.reshape(Vgygrid,(Mgrid)) #Reshaping array
        Vggrid=np.reshape(Vggrid,(Mgrid)) #Reshaping array
        Rgrid=np.reshape(Rgrid,(Mgrid)) #Reshaping array
        thetagrid=np.reshape(thetagrid,(Mgrid)) #Reshaping array
        thetaVgrid=np.reshape(thetaVgrid,(Mgrid)) #Reshaping array
        thetaVgridtan=np.reshape(thetaVgridtan,(Mgrid)) #Reshaping array
    
    #--------------------------------------------------------------------------
    #Output
    return Vxgrid, Vygrid, Vgrid, Vgxgrid, Vgygrid, Vggrid, Rgrid, RVmax, thetagrid, thetaVgrid, thetaVgridtan

    #--------------------------------------------------------------------------
