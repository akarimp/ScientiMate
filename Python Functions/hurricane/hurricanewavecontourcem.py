def hurricanewavecontourcem(xgrid, ygrid, xCenter, yCenter, Hmax, RVmax, VtAzmdir=0, CalcMethod='cem', distCalcMethod='gc', dispout='no'):
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

    scientimate.hurricanewavecontourcem
    ===================================

    .. code:: python

        Hgrid = scientimate.hurricanewavecontourcem(xgrid, ygrid, xCenter, yCenter, Hmax, RVmax, VtAzmdir=0, CalcMethod='cem', distCalcMethod='gc', dispout='no')

    Description
    -----------

    | Calculate hurricane wave height field (contours) on given mesh using 
    | method from Shore Protection Manual (SPM), U.S. Army Corps of Engineers (1984), and Young (1988)
    | For method from Young (1988) use hurricanewavecontoury88
    | For method from Hwang (2016) and Hwang & Walsh (2016) use hurricanewavecontourh16

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
    Hmax
        Hurricane maximum wave height
    RVmax
        Distance (radius) from hurricane center to a location of maximum hurricane wind velocity (m)
    VtAzmdir=0
        | Hurricane center velocity azimuth (bearing) direction in (Degree)
        | azimuth (bearing) direction which is measured clockwise from the north:
        | 0 (degree): toward North, 90 (degree): toward East, 180 (degree): toward South, 270 (degree): toward West 
    CalcMethod='cem'
        | Hurricane wave contour calculation method 
        | 'spm': Use method by Shore Protection Manual (SPM),
        |     U.S. Army Corps of Engineers (1984) in deep water
        | 'cem': Use method by Coastal Engineering Manual (CEM),
        |     U.S. Army Corps of Engineers (2015)
        |     For hurricane with maximum velocity of Vmax>40 (m/s) and translational velocity of Vt<12 (m/s)
        | 'cemnext': Use method by Coastal Engineering Manual (CEM),
        |     U.S. Army Corps of Engineers (Next Release)
        |     For hurricane with maximum velocity of Vmax=40 (m/s) and translational velocity of Vt=7.5 (m/s)
        | 'youngvt2': Use method by Young (1988)
        |     For hurricane with maximum velocity of Vmax=40 (m/s) and translational velocity of Vt=2.5 (m/s)
        | 'youngvt5': Use method by Young (1988)
        |     For hurricane with maximum velocity of Vmax=40 (m/s) and translational velocity of Vt=5 (m/s)
        | 'youngvt7': Use method by Young (1988)
        |     For hurricane with maximum velocity of Vmax=40 (m/s) and translational velocity of Vt=7.5 (m/s)
        | 'youngvt10': Use method by Young (1988)
        |     For hurricane with maximum velocity of Vmax=40 (m/s) and translational velocity of Vt=10 (m/s)
    distCalcMethod='gc'
        | Distance calculation method 
        | 'cart': Distances are calculated on cartesian coordinate
        | 'gc': Distances are calculated on Great Circle based on Vincenty formula, Vincenty (1975)
        | Earth radius coonsidered as mean earth radius=6371000 m
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

    Hgrid
        Hurricane wave height on grid mesh

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np


        #EXAMPLE 1

        #Creating calculation mesh
        xgrid,ygrid=np.meshgrid(np.linspace(-98,-68,100),np.linspace(16,44,100))

        #Longitude of Hurricane Katrine center at max velocity
        longCenter=-88.6

        #Latitude of Hurricane Katrine center at max velocity
        latCenter=26.3

        #Hurricane Katrina maximum significant wave height (m) at max velocity
        Hmax=24.9821

        #Hurricane Katrina radius from hurricane center to a location of maximum hurricane wind velocity (m) at max velocity
        RVmax=6.2750e+004

        #Hurricane Katrina velocity azimuth (bearing) in (Degree) at max velocity
        VtAzmdir=306.76219

        Hgrid=sm.hurricanewavecontourcem(xgrid,ygrid,longCenter,latCenter,Hmax,RVmax,VtAzmdir,'cem','gc','contour')


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

        #Hurricane Katrina maximum significant wave height
        Hmax=[0,0,0,4.3788,4.9295,5.5527,6.2110,6.8516,7.5428,9.1513,8.5021,8.6332,10.1511,11.3434,\
            12.3171,13.5606,14.1226,14.4931,14.1972,19.9683,24.0121,24.9821,23.0419,19.9342,16.5366,\
            14.5246,14.8050,0,0,0,0,0,0,0]

        #Hurricane Katrina radius from hurricane center to a location of maximum hurricane wind velocity (m)
        RVmax=[0,0,0,8.0290e+004,5.6029e+004,4.2063e+004,3.6769e+004,3.3849e+004,3.1352e+004,3.3405e+004,3.3773e+004,\
            3.2657e+004,3.1122e+004,2.7037e+004,2.6512e+004,3.3476e+004,3.0881e+004,4.0266e+004,3.2433e+004,\
            5.1747e+004,5.7297e+004,6.2750e+004,5.3376e+004,4.3074e+004,3.1790e+004,4.3114e+004,2.7800e+004,\
            0,0,0,0,0,0,0]

        #Hurricane Katrina velocity azimuth (bearing) in (Degree)
        VtAzmdir=[0.00000,298.67291,311.22135,338.70264,338.13626,309.94476,279.18860,280.65053,270.13245,\
        246.10095,240.96690,241.20181,244.79591,249.93382,244.88325,252.71384,270.14459,280.49918,\
        298.94148,299.05364,299.18896,306.76219,329.36839,340.59069,0.00000,0.00000,0.00,\
            0.00000,15.67775,15.42254,18.00215,29.63266,39.49673,50.29744]

        Hgrid=sm.hurricanewavecontourcem(xgrid,ygrid,longtrack[3:27],lattrack[3:27],Hmax[3:27],RVmax[3:27],VtAzmdir[3:27],'cem','gc','contour')


    References
    ----------

    Data
    
    * www.nhc.noaa.gov/data/
    * www.nhc.noaa.gov/data/hurdat/hurdat2-format-nencpac.pdf
    * coast.noaa.gov/hurricanes
    * www.aoml.noaa.gov/hrd/data_sub/re_anal.html

    Department of the Army, Waterways Experiment Station, Corps of Engineers, 
    and Coastal Engineering Research Center (1984), 
    Shore Protection Manual, Washington, 
    D.C., vol. 1, 4th ed., 532 pp.

    U.S. Army Corps of Engineers (2015). 
    Coastal Engineering Manual. 
    Engineer Manual 1110-2-1100, Washington, D.C.: U.S. Army Corps of Engineers.

    Young, I. R. (1988). 
    Parametric hurricane wave prediction model. 
    Journal of Waterway, Port, Coastal, and Ocean Engineering, 114(5), 637-652.

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
    from scipy import interpolate
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
    Hmax=type2numpy(Hmax)
    RVmax=type2numpy(RVmax)
    VtAzmdir=type2numpy(VtAzmdir)

    #--------------------------------------------------------------------------
    #Defining required functions

    def cart2pol(x,y):
        rho=np.sqrt(x**2+y**2)
        theta=np.arctan2(y,x)
        return theta,rho
    
    def pol2cart(theta,rho):
        x=rho*np.cos(theta)
        y=rho*np.sin(theta)
        return x,y

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
        xgrid_R_RVmax=np.zeros((Mgrid,Ngrid,len(xCenter))) #Pre-assigning array
        ygrid_R_RVmax=np.zeros((Mgrid,Ngrid,len(xCenter))) #Pre-assigning array
        xgrid_RVmax=np.zeros((Mgrid,Ngrid,len(xCenter))) #Pre-assigning array
        ygrid_RVmax=np.zeros((Mgrid,Ngrid,len(xCenter))) #Pre-assigning array
        Hgrid=np.zeros((Mgrid,Ngrid,len(xCenter))) #Pre-assigning array
    else:
        Rgrid=np.zeros((Mgrid,Ngrid,1)) #Pre-assigning array
        thetagrid=np.zeros((Mgrid,Ngrid,1)) #Pre-assigning array
        xgrid_R_RVmax=np.zeros((Mgrid,Ngrid,1)) #Pre-assigning array
        ygrid_R_RVmax=np.zeros((Mgrid,Ngrid,1)) #Pre-assigning array
        xgrid_RVmax=np.zeros((Mgrid,Ngrid,1)) #Pre-assigning array
        ygrid_RVmax=np.zeros((Mgrid,Ngrid,1)) #Pre-assigning array
        Hgrid=np.zeros((Mgrid,Ngrid,1)) #Pre-assigning array

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
    #Contour (isobar) of hurricane waves based on Shore Protection Manual (SPM), U.S. Army Corps of Engineers (1984)
    #Contour are for hurricane moving toward North

    #Contour (isobar) of hurricane waves based on Shore Protection Manual (SPM), U.S. Army Corps of Engineers (1984)
    if CalcMethod=='spm':
    
        #x of contour in form of R/RVmax ratio
        xcontour_R_RVmax=[-2.1397,-2.1078,-2.0751,-2.0069,-1.8662,-1.7955,-1.5478,-1.0859,-0.5517,\
            -0.0170,0.4113,0.9826,1.3758,1.8768,2.1277,2.3427,2.6656,3.0248,\
            3.5284,3.9606,4.6095,5.0072,5.4420,5.6236,5.8418,5.9881,6.1715,\
            6.2475,6.2145,6.0751,5.8290,5.6878,5.3338,5.0500,4.7304,4.0536,\
            3.4472,2.9829,2.5177,2.2312,1.8367,1.2983,0.5791,0.0033,-0.6456,\
            -1.0064,-1.5124,-1.8026,-2.0218,-2.0967,-2.1368,-1.6379,-1.5707,-1.4661,\
            -1.2545,-0.9354,-0.4731,0.0969,0.5610,0.9900,1.5269,1.8853,2.4237,\
            2.9629,3.5751,3.9720,4.3336,4.6235,4.8784,5.0257,5.0644,4.9967,\
            4.8568,4.5743,4.2907,3.9711,3.6506,3.3655,2.9372,2.5799,2.0072,\
            1.8283,1.5415,1.1109,0.6081,-0.0036,-0.5447,-0.9785,-1.3411,-1.4874,\
            -1.5989,-0.5627,-0.5279,-0.4211,-0.2783,-0.0266,0.0459,0.0825,0.0121,\
            -0.1304,-0.2377,-0.4174,-0.5258,-1.1361,-1.1032,-0.9991,-0.7166,-0.3261,\
            0.0302,0.4943,1.0667,1.3893,1.7483,1.9999,2.4316,2.7921,3.0088,\
            3.4428,3.5878,3.8058,3.9892,3.8501,3.6379,3.3185,2.9986,2.8207,\
            2.3574,1.9648,1.5713,1.4282,1.2131,1.1056,0.7467,0.2436,-0.0086,\
            -0.3692,-0.7310,-0.9861,-1.0966,-0.8494,-0.8149,-0.6731,-0.4233,-0.2799,\
            0.0084,0.0812,0.1183,0.0838,-0.0224,-0.2005,-0.3795,-0.5951,-0.7756,\
            0.2259,0.1884,0.1154,-0.0656,-0.2820,-0.3906,-0.4272,-0.3566,-0.1431,\
            0.1420,0.4633,0.7850,0.9998,1.1790,1.4302,1.6816,1.9694,2.2216,\
            2.5104,2.7276,2.8728,3.0191,3.0573,3.0239,2.9187,2.7419,2.5645,\
            2.3158,1.9950,1.5663,1.2084,0.9933,0.6702,0.2747,-0.0495,-0.3384,\
            -0.4108,-0.4475,-0.4122,-0.1987,-0.0208,0.1563,0.2975,0.2965,0.2232,\
            0.0058,-0.1755,-0.1766,0.0010,0.3937,0.7513,0.8588,1.0739,1.3970,\
            1.7567,2.0450,2.5155,2.6254,2.6631,2.6291,2.5239,2.2760,1.9558,\
            1.6348,1.2774,1.0985,0.8834,0.7401,0.4887,0.1289,-0.0875,-0.1966,\
            -0.1971,-0.0909,0.1223,0.2636,0.4051,0.3674,0.2214,0.0037,0.0029,\
            0.2166,0.5025,0.6815,0.8248,1.0043,1.2916,1.5435,1.7957,1.9404,\
            2.1945,2.2688,2.2351,2.0585,1.8100,1.5249,1.2034,0.9888,0.8813,\
            0.6304,0.3793,0.1634,-0.0171,-0.0179,0.1595,0.3360,0.4767,0.4749,\
            0.4372,0.3276,0.2546,0.3252,0.5397,0.6831,1.0064,1.2223,1.4747,\
            1.8002,1.9093,1.9462,1.8764,1.6993,1.3791,1.0223,0.7005,0.5213,\
            0.3418,0.2335,0.2686,0.4093,0.4781,0.5843,0.5798,0.5062,0.6130,\
            0.7563,1.0080,1.2966,1.5504,1.5878,1.5183,1.2703,0.9851,0.6992,\
            0.5917,0.5195,0.5187,0.5890,0.7635,0.7598,0.7946,0.8301,0.9379,\
            1.0826,1.2284,1.2653,1.2310,1.1251,0.9830,0.9119,0.8405,0.7685,\
            0.7677,1.0000]
    
        #y of contour in form of R/RVmax ratio
        ycontour_R_RVmax=[-0.0153,0.5384,0.9813,1.4612,1.8308,1.9786,2.4223,2.9775,3.4592,3.8671,4.1270,\
            4.4243,4.5733,4.6859,4.6868,4.6876,4.6519,4.5425,4.2861,3.9925,3.4783,3.0000,\
            2.3374,2.0060,1.5640,1.1586,0.5689,-0.0212,-0.4273,-0.9813,-1.6464,-1.9421,-2.5707,\
            -2.9777,-3.3848,-3.9777,-4.3858,-4.6089,-4.7213,-4.7593,-4.7238,-4.6151,-4.2857,-3.9557,\
            -3.4415,-3.1107,-2.5222,-2.0436,-1.4540,-1.0114,-0.4212,-0.0134,0.6141,1.0204,1.5009,\
            1.9818,2.5001,2.9819,3.2419,3.3911,3.5038,3.5051,3.3964,3.1770,2.7733,2.4058,\
            1.9643,1.5226,0.9700,0.4171,0.0113,-0.5425,-1.0227,-1.6141,-2.0580,-2.4651,-2.7614,\
            -2.9839,-3.2438,-3.3927,-3.5055,-3.5431,-3.5441,-3.4719,-3.3262,-2.9963,-2.5186,-2.0036,\
            -1.4145,-1.0092,-0.4561,-0.0095,0.1383,0.2494,0.3237,0.2139,0.1035,-0.0071,-0.1918,\
            -0.3031,-0.3404,-0.2672,-0.1569,-0.0116,0.3945,0.8745,1.4660,1.9840,2.2806,2.5406,\
            2.6903,2.6915,2.6190,2.5092,2.2894,1.9955,1.7749,1.2230,1.0022,0.5971,0.0073,\
            -0.5836,-0.9903,-1.4342,-1.8044,-1.9896,-2.3603,-2.5832,-2.6953,-2.7327,-2.7335,-2.7339,\
            -2.6615,-2.4788,-2.2952,-2.0014,-1.5230,-0.9335,-0.5280,-0.0105,0.1741,0.3960,0.5446,\
            0.5451,0.3248,0.1774,-0.0069,-0.1916,-0.3765,-0.5247,-0.5623,-0.4893,-0.3054,-0.0065,\
            0.2147,0.3990,0.6566,0.8403,0.9875,1.0981,1.2459,1.4681,1.6906,1.8762,1.9881,\
            2.0258,2.0265,1.9905,1.9176,1.7711,1.5875,1.2934,0.9990,0.7412,0.3359,0.0039,\
            -0.3283,-0.6608,-0.9936,-1.2525,-1.5487,-1.8081,-1.9942,-2.0693,-2.0701,-1.9975,-1.8145,\
            -1.5943,-1.3001,-1.1897,-1.0791,-1.0052,-0.7830,-0.5979,-0.3020,-0.0063,0.1413,0.3625,\
            0.6938,0.9883,1.1359,1.3579,1.5808,1.6928,1.6932,1.6940,1.6214,1.4382,1.2179,\
            0.5923,0.2606,0.0024,-0.2560,-0.5885,-0.9953,-1.3286,-1.5512,-1.7001,-1.7376,-1.7384,\
            -1.7390,-1.6661,-1.4829,-1.2992,-1.0782,-1.0044,-0.8195,-0.5604,-0.2647,-0.0059,0.2523,\
            0.6207,0.9889,1.0996,1.2849,1.3967,1.4342,1.4348,1.3985,1.3258,1.1791,0.9956,\
            0.8116,0.3697,0.0010,-0.2943,-0.6640,-0.9970,-1.2195,-1.3683,-1.4428,-1.4432,-1.4442,\
            -1.4082,-1.2983,-1.1145,-1.0038,-0.7448,-0.3751,-0.0056,0.2527,0.5108,0.8056,0.9899,\
            1.1377,1.2123,1.2129,1.1033,0.9934,0.7730,0.3683,0.1473,-0.0002,-0.2588,-0.5546,\
            -0.8879,-1.1106,-1.2225,-1.2232,-1.1869,-1.0766,-0.9658,-0.5963,-0.1901,-0.0052,0.6221,\
            0.8801,0.9912,0.9917,0.8820,0.6247,0.2197,-0.0015,-0.2970,-0.7038,-0.9263,-1.0380,\
            -1.0384,-0.9649,-0.8542,-0.6694,-0.0046,0.5120,0.6598,0.6968,0.6603,0.4763,0.1448,\
            -0.0027,-0.2242,-0.4460,-0.6311,-0.7051,-0.7423,-0.70570,-0.5950,0.0000]
    
        #Contour values in form of H/Hmax ratio
        zcontour_H_Hmax=[0.40,0.40,0.40,0.40,0.40,0.40,0.40,0.40,0.40,0.40,0.40,\
            0.40,0.40,0.40,0.40,0.40,0.40,0.40,0.40,0.40,0.40,0.40,\
            0.40,0.40,0.40,0.40,0.40,0.40,0.40,0.40,0.40,0.40,0.40,\
            0.40,0.40,0.40,0.40,0.40,0.40,0.40,0.40,0.40,0.40,0.40,\
            0.40,0.40,0.40,0.40,0.40,0.40,0.40,0.50,0.50,0.50,0.50,\
            0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,\
            0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,\
            0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,\
            0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,\
            0.50,0.50,0.50,0.50,0.60,0.60,0.60,0.60,0.60,0.60,0.60,\
            0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,\
            0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,\
            0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,\
            0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.70,\
            0.70,0.70,0.70,0.70,0.70,0.70,0.70,0.70,0.70,0.70,0.70,\
            0.70,0.70,0.70,0.70,0.70,0.70,0.70,0.70,0.70,0.70,0.70,\
            0.70,0.70,0.70,0.70,0.70,0.70,0.70,0.70,0.70,0.70,0.70,\
            0.70,0.70,0.70,0.70,0.70,0.70,0.70,0.70,0.75,0.75,0.75,\
            0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,\
            0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,\
            0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.80,0.80,\
            0.80,0.80,0.80,0.80,0.80,0.80,0.80,0.80,0.80,0.80,0.80,\
            0.80,0.80,0.80,0.80,0.80,0.80,0.80,0.80,0.80,0.80,0.80,\
            0.80,0.80,0.80,0.80,0.80,0.80,0.85,0.85,0.85,0.85,0.85,\
            0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,\
            0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.90,0.90,\
            0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90,\
            0.90,0.90,0.90,0.90,0.95,0.95,0.95,0.95,0.95,0.95,0.95,\
            0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,1.00]
    
        xcontour_R_RVmax=np.array(xcontour_R_RVmax)
        ycontour_R_RVmax=np.array(ycontour_R_RVmax)
        zcontour_H_Hmax=np.array(zcontour_H_Hmax)

    #--------------------------------------------------------------------------
    #Contour (isobar) of hurricane waves based on Coastal Engineering Manual (CEM), U.S. Army Corps of Engineers (2015)
    #Contour are for hurricane moving toward North

    #Contour (isobar) of hurricane waves based on Coastal Engineering Manual (CEM), U.S. Army Corps of Engineers (2015)
    #For hurricane with maximum velocity of Vmax>40 (m/s) and translational velocity of Vt<12 (m/s)
    if CalcMethod=='cem':
    
        #x of contour in form of R/RVmax ratio
        xcontour_R_RVmax=[-0.2907,-0.2481,-0.2481,-0.1628,-0.0775,0.0078,0.0078,-0.0349,-0.1202,\
            -0.9729,-0.8450,-0.6744,-0.5465,-0.3760,-0.1628,0.0504,0.1783,0.2209,\
            0.1783,0.0930,0.0504,-0.2054,-0.3760,-0.5465,-0.7597,-0.9302,-0.9729,\
            -4.4264,-4.2558,-4.1705,-4.0426,-3.7442,-3.4884,-3.2752,-2.9767,-2.5504,\
            -2.2093,-2.0388,-2.0388,-2.1240,-2.2946,-2.4225,-4.4690,-4.2984,-4.1705,\
            -4.0000,-3.7868,5.6357,5.8915,6.1473,6.3605,6.7016,6.8721,6.9574,\
            7.1705,7.3411,7.4690,-4.2984,-4.2984,-4.2984,-4.2558,-4.2132,-4.0853,\
            -4.0000,-3.9574,-3.8721,-3.8721,-4.0426,-4.1279,-4.2132,-4.4690,-4.3837,\
            -4.2984,-4.1705,-4.0426,-4.0000,-3.8721,-3.7016,-3.4884,-3.3178,-3.1899,\
            -3.0194,-2.7636,-2.5930,-2.5078,-2.4225,-2.4225,-2.4651,-2.4651,-2.2946,\
            -1.9961,-1.6977,-1.3992,-0.7597,-0.4612,0.0504,0.5620,1.0310,1.2868,\
            1.6705,1.9690,2.3101,2.6512,2.9070,3.0775,3.1628,3.1628,4.0155,\
            4.3992,4.6977,4.9535,5.2519,5.3372,5.4651,5.8062,5.8488,6.1047,\
            6.4031,6.5310,6.7868,7.0000,7.0853,7.2132,7.4690,7.4690,7.2984,\
            7.2132,7.0426,6.8721,6.6163,6.2326,5.9341,5.6783,0.0504,-0.1202,\
            -0.2481,-0.2907,-0.3333,-0.6318,-0.7171,-0.8450,-0.8450,-0.7597,-0.5891,\
            -0.2907,0.0078,0.2209,0.4341,0.3062,-0.0349,-0.5465,-0.6744,-0.8876,\
            -1.1008,-1.4845,-1.6124,-1.8682,-2.0814,-2.3372,-2.9341,-3.1473,-3.3178,\
            -3.5736,-3.6163,-3.6163,-3.7016,-3.8721,-4.4690,-2.5930,-2.5504,-2.5504,\
            -2.5930,-2.5504,-2.5078,-2.5078,-2.6357,-3.0194,-3.0194,-2.8915,-2.5930,\
            -2.3798,-2.1240,-1.8256,-1.4419,-1.1860,-1.1008,-0.8876,-0.7171,-0.5465,\
            -0.4186,-0.1202,0.0504,0.3062,0.5620,1.2016,1.7558,1.8411,2.2674,\
            2.5233,2.7791,2.9070,3.0775,3.2481,3.4612,3.5465,3.8023,4.1434,\
            4.4419,4.6977,4.9109,5.0388,5.3372,5.4651,5.5930,5.7209,5.7636,\
            5.8488,6.0194,6.1473,6.2326,6.4031,6.4884,6.4884,6.4457,6.4031,\
            6.3178,6.2752,6.2326,6.1899,6.0194,5.8062,5.4651,5.1240,4.8256,\
            4.4845,4.1860,3.8023,3.3760,2.7364,2.4380,2.0543,1.5853,1.2868,\
            0.9884,0.6899,0.4767,0.3062,0.2636,0.2636,0.2636,0.3488,0.3915,\
            0.4767,0.5620,0.5620,0.5194,0.4341,0.3062,0.1357,-0.1202,-0.2907,\
            -0.5465,-0.8023,-1.1434,-1.3140,-1.4845,-1.6124,-1.8682,-1.9961,-2.1240,\
            -2.3372,-2.5078,-2.5930,-1.4845,-1.4845,-1.2713,-1.3992,-1.2713,-1.0581,\
            -0.8876,-0.6744,-0.4612,-0.2907,-0.0349,0.4341,0.9884,0.9884,0.9031,\
            0.8605,0.9884,1.2442,1.7132,1.9264,2.0969,2.2674,2.4380,2.9070,\
            3.1202,3.3760,3.6744,3.8450,3.9302,3.9729,4.2713,4.3566,4.4845,\
            4.6550,4.8256,4.9961,5.1240,5.2093,5.2519,5.2519,5.2093,5.1240,\
            4.9961,4.8682,4.7403,4.5698,4.3566,4.1008,3.7171,3.3333,2.9922,\
            2.6085,2.2248,2.0116,1.7132,1.4574,1.2442,1.2016,1.1163,0.9031,\
            0.7752,0.6899,0.6899,0.6899,0.6899,0.6473,0.6047,0.4341,0.0930,\
            -0.0775,-0.2481,-0.5039,-0.8876,-1.1434,-1.2713,-1.3992,-0.2481,-0.3760,\
            -0.3760,-0.3333,-0.2481,-0.0775,0.3488,0.8605,1.1163,1.3295,1.5853,\
            1.8411,2.0116,2.0969,2.1822,2.4806,2.8643,3.2054,3.4612,3.5039,\
            3.5039,3.4186,3.3760,3.3760,3.5039,3.6744,3.7171,3.7597,3.7597,\
            3.7171,3.6318,3.4186,3.2481,3.0775,2.8217,2.6085,2.3527,2.0969,\
            1.8837,1.6705,1.5000,1.2868,1.1589,1.0310,0.9031,0.8605,0.7752,\
            0.7326,0.6473,0.4341,0.1783,-0.0775,0.2636,0.2636,0.3488,0.4341,\
            0.5194,0.6473,0.8178,0.9457,1.5000,1.6705,1.8411,2.0969,2.2674,\
            2.3953,2.4380,2.4380,2.3527,2.1822,2.0116,1.8837,1.7558,1.6705,\
            1.5426,1.3295,1.0310,0.9031,0.8605,0.8178,0.6473,0.5194]
    
        #y of contour in form of R/RVmax ratio
        ycontour_R_RVmax=[-0.5907,-0.5049,-0.4620,-0.4623,-0.4626,-0.5059,-0.6348,-0.7206,-0.7203,\
            -0.8888,-0.5456,-0.2884,-0.2030,-0.2037,-0.2475,-0.3342,-0.4207,-0.5497,\
            -0.8074,-1.0649,-1.1936,-1.3645,-1.4068,-1.4061,-1.3193,-1.1468,-1.0177,\
            -2.8949,-3.0674,-3.3256,-3.4550,-3.4561,-3.5431,-3.7587,-4.0177,-4.2342,\
            -4.3645,-4.4940,-4.6229,-4.7515,-4.8798,-4.9222,5.1834,5.3546,5.5690,\
            5.7402,5.9112,5.9174,5.7445,5.5716,5.3989,5.1398,4.8813,4.7091,\
            4.7083,4.8365,4.9220,1.6164,1.7453,1.8742,2.0029,2.0457,2.0452,\
            2.0019,1.8728,1.6147,1.4428,1.1427,1.1860,1.3582,3.5076,3.7651,\
            4.0656,4.3229,4.4943,4.5371,4.4506,4.2781,4.1483,4.1047,4.1042,\
            4.1465,4.2744,4.4027,4.5312,4.9176,5.0895,5.2615,5.4334,5.6046,\
            5.7324,5.8601,5.9019,5.8994,5.8553,5.7244,5.5935,5.5487,5.5047,\
            5.5032,5.5020,5.3288,5.2416,5.2406,5.2829,5.4974,5.9271,5.9237,\
            5.5785,5.1906,4.8888,4.5869,4.4576,4.2423,3.9831,3.8111,3.3804,\
            2.7777,2.5194,2.4754,2.3886,2.4742,3.0323,3.3751,-3.1562,-3.4993,\
            -3.6708,-3.8420,-4.0562,-4.2700,-4.4834,-4.7400,-4.9968,-4.9749,-4.7593,\
            -4.5010,-4.2860,-4.1140,-3.7261,-3.5969,-3.1667,-2.6081,-2.3506,-2.0075,\
            -1.7509,-1.5372,-1.0654,-0.4646,-0.2493,-0.0331,0.0978,0.0553,-0.2876,\
            -0.4157,-0.4142,-0.4566,-0.7564,-0.9704,-1.0554,-1.2679,-1.2670,-1.2234,\
            -1.0505,-0.9644,-0.5777,-0.3625,-0.2330,0.0272,0.5784,0.8361,0.9650,\
            1.3519,1.6095,1.7383,1.8672,2.0395,2.4278,2.5137,2.6421,2.8558,\
            2.9409,2.9829,3.0676,3.2810,3.4089,3.5375,3.9663,4.3094,4.5236,\
            4.6090,4.6938,4.6931,4.6491,4.6052,4.5597,4.4716,4.4713,4.1688,\
            4.0389,3.9520,3.9515,4.0368,4.2080,4.3360,4.3357,4.2058,4.1185,\
            4.0314,3.9015,3.6858,3.4705,3.0397,2.8243,2.6519,2.2647,2.0497,\
            1.8775,1.6620,1.4037,1.0596,0.7152,0.3711,0.1133,-0.1873,-0.6169,\
            -1.0032,-1.2609,-1.7334,-2.0770,-2.6349,-2.9778,-3.3202,-3.6197,-3.7904,\
            -4.0039,-4.1746,-4.2590,-4.3433,-4.3408,-4.2967,-4.1663,-3.9496,-3.7336,\
            -3.4316,-3.1297,-2.9140,-2.6555,-2.4405,-2.1397,-1.9679,-1.5815,-1.3238,\
            -1.0234,-0.6370,-0.2933,-0.0353,0.0940,0.2663,0.3959,0.5258,0.4835,\
            0.3556,0.1847,0.0572,0.1008,0.1874,0.1879,0.0170,-0.1544,-0.1539,\
            0.0618,0.2344,0.4066,0.9179,1.0897,1.5616,1.9488,2.1631,2.3771,\
            2.5054,2.5905,2.7186,2.7609,2.8028,2.9729,2.9707,3.0996,3.3148,\
            3.4439,3.5293,3.5283,3.4835,3.3967,3.2672,3.0087,2.8791,2.8343,\
            2.7475,2.5317,2.3157,2.1431,1.9280,1.7989,1.5829,1.4536,1.3242,\
            1.4095,1.5377,1.5371,1.3217,0.8917,0.4189,0.2040,-0.2684,-0.7408,\
            -1.1270,-1.5562,-1.9424,-2.2425,-2.5425,-2.7993,-3.0126,-3.1400,-3.1817,\
            -3.1802,-3.0927,-3.0060,-2.8329,-2.5311,-2.1436,-2.0145,-1.8853,-1.6696,\
            -1.4113,-1.0672,-0.8523,-0.6375,-0.4656,-0.1647,0.0074,0.2658,0.6109,\
            0.7405,0.6982,0.6133,0.4429,0.4439,0.5733,0.7027,0.9560,1.0854,\
            1.2143,1.3431,1.5576,1.7717,1.8990,1.9829,2.0679,2.1530,2.1520,\
            2.0650,1.8925,1.7632,1.6770,1.5469,1.3735,1.2433,1.1564,1.1132,\
            1.0273,0.8557,0.6840,0.4692,0.2109,-0.0046,-0.2196,-0.3487,-0.4776,\
            -0.7353,-1.0787,-1.4646,-1.7217,-1.8929,-2.0208,-2.0200,-1.9760,-1.9750,\
            -1.9312,-1.8445,-1.7149,-1.4992,-1.3268,-1.1974,-0.8962,-0.6811,-0.2082,\
            -0.0361,0.1791,0.4377,0.7395,0.9124,0.9540,1.0829,1.2115,1.2541,\
            1.2538,1.2103,1.0807,1.0373,1.0351,0.9485,0.7760,0.6031,0.4306,\
            0.2152,-0.0857,-0.3006,-0.6440,-0.9871,-1.2442,-1.2437,-1.1573,-0.9851,\
            -0.9846,-0.8548,-0.5959,-0.3805,-0.1225,0.0925,0.4798,0.6952]
    
        #Contour values in form of H/Hmax ratio
        zcontour_H_Hmax=[0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.4,0.4,\
            0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,\
            0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,\
            0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,\
            0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,\
            0.4,0.4,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.7,0.7,\
            0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,\
            0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,\
            0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,\
            0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,\
            0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,\
            0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,\
            0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.8,0.8,0.8,\
            0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,\
            0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,\
            0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,\
            0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,\
            0.8,0.8,0.8,0.8,0.8,0.9,0.9,0.9,0.9,0.9,0.9,\
            0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,\
            0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,\
            0.9,0.9]
    
        xcontour_R_RVmax=np.array(xcontour_R_RVmax)
        ycontour_R_RVmax=np.array(ycontour_R_RVmax)
        zcontour_H_Hmax=np.array(zcontour_H_Hmax)

    #--------------------------------------------------------------------------
    #Contour (isobar) of hurricane waves based on Coastal Engineering Manual (CEM), U.S. Army Corps of Engineers
    #Contour extracted from a figure prepared by Ian Young for next release of Coastal Engineering Manual (CEM)
    #Contour are for hurricane moving toward North

    #Contour (isobar) of hurricane waves based on Coastal Engineering Manual (CEM), U.S. Army Corps of Engineers
    #Contour extracted from a figure prepared by Ian Young for next release of Coastal Engineering Manual (CEM)
    #For hurricane with maximum velocity of Vmax=40 (m/s) and translational velocity of Vt=7.5 (m/s)
    if CalcMethod=='cemnext':

        #x of contour in form of x/RVmax ratio
        xcontour_x_RVmax=[-4.5052,-4.3674,-4.2526,-4.1378,-3.4029,-3.0125,-2.6681,-2.3236,-2.0021,\
            -1.5887,-1.1524,-0.8309,-0.7161,4.0605,4.1294,4.2213,4.4739,4.7265,\
            5.5073,5.6910,5.9207,6.1503,6.3570,6.6096,6.8392,7.0000,7.2067,\
            7.4823,-4.5052,-4.2985,-3.8852,-3.7933,-3.6096,-3.3800,-3.0125,-2.7140,\
            -1.5887,-1.1983,-0.8998,-0.6013,-0.2797,0.0877,0.3633,0.5010,1.0981,\
            1.3507,1.4885,1.5804,1.6952,1.9478,2.2004,2.3841,2.6827,3.0731,\
            3.4635,3.7161,3.9687,4.1983,4.5428,4.8184,5.1858,5.9436,6.1503,\
            6.3800,7.0000,7.2985,7.4134,7.4823,-0.9916,-1.2213,-1.4050,-1.4969,\
            -1.4050,-1.4050,-2.1858,-2.3236,-2.5992,-2.8977,-2.5992,-2.4384,-2.2547,\
            -2.0710,-1.8873,-1.7954,-1.6806,-1.7954,-2.1628,-2.3006,-2.4154,-2.4843,\
            -2.6451,-2.8058,-2.8518,-2.9207,-2.9896,-3.1733,-3.4259,-3.4259,-3.5177,\
            -3.6555,-3.8622,-4.0230,-4.2296,-4.3904,-4.5052,-4.5052,-4.4134,-4.3445,\
            -4.2067,-3.7474,-3.5637,-3.4029,-3.2422,-3.0125,-2.5992,-2.0021,-1.7724,\
            -1.5887,-1.4739,-1.2672,-0.9916,-0.8079,-0.2109,0.2255,0.7766,1.0522,\
            1.3507,1.4885,1.5804,1.7411,1.9478,2.1086,2.4071,2.7286,2.9353,\
            3.1649,3.3486,3.6701,3.9457,4.2213,4.5887,5.3236,5.3925,5.8058,\
            6.1273,6.4718,6.5866,7.0919,7.2985,7.3904,7.4823,7.4823,7.2067,\
            7.0230,6.9081,0.8914,0.7077,0.6159,0.5929,0.5929,0.5240,0.5240,\
            0.5010,0.3862,0.3862,0.4092,0.4092,0.3633,0.2255,0.1106,-0.0042,\
            -0.0960,-0.3257,-0.4864,-0.7850,-0.9228,-1.0146,-1.2213,-1.3361,-1.4969,\
            -1.8643,-2.0021,-2.2088,-3.0125,-3.1962,-3.5637,-3.9081,-4.0000,-4.0919,\
            -4.2296,-4.3445,-4.4363,-4.4593,-4.4823,-4.5052,-4.1608,-4.1148,-4.0230,\
            -3.7933,-3.6555,-3.5177,-3.3800,-3.2881,-3.1733,-3.1503,-2.9896,-2.8747,\
            -2.6910,-2.5992,-2.4154,-2.2088,-2.0480,-1.9102,-1.6806,-1.1524,-0.3257,\
            0.1795,0.4781,0.7077,0.9374,1.2129,1.3737,1.3967,1.4885,1.5804,\
            1.8559,2.2234,2.5449,2.9123,3.0271,3.0960,3.2338,3.3946,4.1524,\
            4.4739,4.7035,4.9791,5.2088,5.3925,5.6681,5.9207,6.0814,6.2881,\
            6.4029,6.5407,6.6096,6.6555,6.6785,6.6785,6.6555,6.3800,6.1733,\
            5.9436,5.6221,5.5762,5.2088,4.8643,4.5658,4.1065,3.9687,3.7850,\
            3.6013,3.1190,2.7516,2.6827,2.5678,1.8330,1.6722,1.6033,1.4885,\
            1.3967,1.1900,1.1670,1.0981,1.0752,1.0292,1.0292,0.8914,0.8685,\
            0.7996,0.8225,0.8225,0.7996,0.6618,0.3403,0.1566,0.0877,-0.0042,\
            0.0877,0.2025,0.2944,0.3403,0.3633,0.3173,0.2025,0.1106,-0.0042,\
            -0.0960,-0.1649,-0.5324,-0.6013,-0.6242,-0.7390,-0.9228,-1.1065,-1.2213,\
            -1.2672,-1.4280,-1.6117,-1.8184,-2.0021,-2.2088,-2.5532,-2.7829,-2.8977,\
            -3.0355,-3.1503,-3.1962,-3.3800,-3.5177,-3.7015,-3.9081,-4.0000,-4.1148,\
            -4.1378,-1.7265,-1.8413,-1.8873,-1.8873,-1.8643,-1.7954,-1.4509,-1.3361,\
            -1.0835,-0.7161,-0.6472,-0.4864,-0.2568,-0.0960,0.1106,0.3862,0.5010,\
            0.6848,0.8914,1.0752,1.8100,2.0856,2.3152,2.5449,2.7975,3.0960,\
            3.4175,3.7850,4.0835,4.3820,4.6806,4.8643,5.2088,5.3925,5.4843,\
            5.5762,5.8288,5.8288,5.6910,5.5992,5.2088,5.0939,5.0251,4.7035,\
            4.3361,4.0835,3.8768,3.7161,3.5094,3.2568,2.9582,2.6827,1.9248,\
            1.6033,1.3967,1.2589,1.1441,1.0063,0.7077,0.6159,0.5699,0.5470,\
            0.5470,0.4781,0.2025,-0.0042,-0.1879,-0.3946,-0.6013,-0.8079,-0.8998,\
            -0.9687,-0.6701,-0.6472,-0.3716,-0.0731,0.7077,1.0752,1.3278,1.4885,\
            1.6033,1.9937,2.2234,2.4990,2.6827,2.7745,3.0042,3.2797,3.4635,\
            3.6472,3.7850,3.8998,3.9687,4.0146,4.0835,4.2443,4.3591,4.3361,\
            4.2672,4.1983,4.0835,3.8998,3.6931,3.4635,3.3027,3.0042,2.8434,\
            2.7516,1.9478,1.6952,1.5344,1.4196,1.1211,0.8685,0.7307,0.6848,\
            0.6618,0.5929,0.4092,0.2944,0.2025,-0.0042,-0.2338,-0.4175,-0.5783,\
            0.2714,0.4092,0.5240,0.5929,0.7766,0.9833,1.0981,1.3048,1.3967,\
            1.5344,1.6722,1.8330,2.0167,2.3841,2.4990,2.8894,3.0042,3.0960,\
            3.1879,3.3716,3.6013,3.5783,3.4864,3.2797,3.0042,2.8205,2.6597,\
            2.3152,1.8789,1.6493,1.6033,1.3507,1.0522,0.7077,0.5929,0.5240,\
            0.5010,0.4092,0.9603,0.9833,1.0752,1.2129,1.4426,1.6952,1.8559,\
            1.9248,1.9248,1.8789,1.7871,1.4885,1.2129,1.0981,1.0063,0.9603,\
            0.8685]

        #y of contour in form of y/RVmax ratio
        ycontour_y_RVmax=[4.5532,4.7599,4.9896,5.3340,5.5637,5.6555,5.7244,5.7933,5.8392,\
            5.8852,5.9311,5.9541,5.9770,5.9770,5.9770,5.9541,5.9081,5.8622,\
            5.6555,5.5866,5.4718,5.3570,5.1962,5.0355,4.8977,4.8058,4.6221,\
            4.3695,2.4175,2.9687,3.7495,3.7495,3.7265,3.7495,3.8184,3.8873,\
            3.8873,3.9332,4.0021,4.2088,4.3466,4.4614,4.5303,4.5992,5.2192,\
            5.3570,5.4259,5.4259,5.3800,5.2422,5.1044,4.9666,4.8518,4.6910,\
            4.5303,4.4384,4.3466,4.1858,3.9562,3.7954,3.5887,3.2672,3.1754,\
            3.0146,2.8079,2.5324,2.4175,2.3486,-5.0000,-4.8163,-4.7015,-4.7015,\
            -4.6096,-4.4718,-4.4259,-4.4718,-4.4948,-4.5177,-4.5407,-4.5637,-4.5637,\
            -4.5637,-4.6096,-4.6785,-4.7015,-4.7244,-4.7244,-4.6785,-4.6785,-4.7015,\
            -4.7244,-4.7704,-4.8163,-4.8392,-4.9081,-5.0000,-5.0000,-4.6096,-4.5177,\
            -4.4718,-4.4489,-4.4489,-4.4948,-4.5407,-4.5407,-0.0626,0.2589,0.4885,\
            0.7641,1.5219,1.8205,2.0501,2.3486,2.6242,2.8309,3.0146,3.0146,\
            3.0376,3.0376,3.0146,3.0146,2.9916,2.9916,3.0605,3.6576,4.0710,\
            4.3925,4.4614,4.4614,4.3925,4.2317,4.1169,3.8873,3.7724,3.7035,\
            3.5658,3.4280,3.2443,3.0605,2.8768,2.6931,2.3027,2.2338,2.0731,\
            1.8894,1.6367,1.5678,1.1086,0.8100,0.6952,0.5115,-4.2881,-4.6326,\
            -4.9081,-5.0230,-5.0230,-4.7474,-4.5177,-4.4259,-4.3111,-3.3925,-3.0021,\
            -2.8184,-2.6117,-2.3591,-2.2213,-2.0835,-1.9457,-1.6931,-1.5094,-1.4405,\
            -1.4405,-1.5553,-1.6242,-2.1065,-2.2902,-2.4050,-2.4739,-2.4969,-2.5658,\
            -2.3591,-2.3132,-2.3361,-2.5428,-2.6117,-2.7724,-2.9102,-2.9102,-2.8873,\
            -2.7495,-2.5428,-2.3820,-2.2213,-1.9228,-1.6701,-1.0501,-0.8894,-0.6597,\
            -0.5219,-0.2923,-0.0626,0.3048,0.5574,0.7871,0.8559,1.0856,1.2923,\
            1.4990,1.5908,1.8664,2.0960,2.2109,2.2568,2.3027,2.3946,2.4175,\
            2.3946,2.3716,2.7390,3.0835,3.4280,3.5887,3.6347,3.7035,3.6806,\
            3.5198,3.3361,3.1754,2.9687,2.8998,2.8309,2.7390,2.6242,2.1879,\
            1.9812,1.7975,1.6367,1.5219,1.4071,1.1775,0.9937,0.7641,0.4885,\
            0.2818,-0.0626,-0.3382,-0.6597,-0.7516,-1.4864,-1.5783,-2.4280,-2.7724,\
            -3.0939,-3.5073,-3.5992,-3.8977,-4.0814,-4.2192,-4.3570,-4.3800,-4.3800,\
            -4.3800,-4.2422,-4.1044,-4.0585,-4.0355,-4.3800,-4.5177,-4.6785,-4.8163,\
            -4.7474,-4.4259,-4.3111,-4.0125,-3.6910,-3.3925,-2.9102,-2.5428,-2.4050,\
            -1.8998,-1.8309,-1.6931,-1.5094,-1.3946,-1.1420,-0.9812,-0.8434,-0.7056,\
            -0.6597,-0.5678,-0.4530,-0.3152,-0.2004,-0.0856,0.0292,0.0981,0.0981,\
            0.0292,-0.0626,-0.6597,-0.7745,-0.8894,-1.0960,-1.3027,-1.3946,-1.4175,\
            -1.3946,-1.3946,-1.3257,-1.2109,-1.0960,-1.0501,-1.0501,-1.0960,-1.1190,\
            -1.2109,-1.3257,-1.4175,-1.5553,-1.6931,-1.9687,-2.1983,-2.2443,-2.0835,\
            -1.9687,-0.0397,0.0981,0.4426,0.7182,0.9019,1.0856,1.4760,1.5678,\
            1.7745,1.9123,1.9353,1.9353,1.9582,1.9582,1.9353,1.8664,1.8664,\
            2.0960,2.4175,2.6472,2.8768,2.7850,2.7390,2.6472,2.4635,2.2568,\
            2.0042,1.7745,1.5678,1.0167,0.7641,0.6952,0.4885,0.4426,0.3967,\
            0.2818,-0.5678,-0.7056,-1.0271,-1.1420,-1.5094,-1.7620,-1.8768,-2.2443,\
            -2.5198,-2.7035,-2.7724,-2.8184,-2.8184,-2.7724,-2.6806,-2.5658,-2.3132,\
            -2.1754,-1.9916,-1.6701,-1.4635,-1.3027,-1.0960,-0.9812,-0.7975,-0.6827,\
            -0.0856,0.0752,0.3967,0.5115,0.4656,0.3278,0.2129,0.0063,-0.1315,\
            -0.2004,0.9937,1.0856,1.2923,1.4530,1.6367,2.0271,2.1879,2.2568,\
            2.2568,2.1420,2.1420,2.1190,1.9812,1.9582,1.8205,1.5908,1.4530,\
            1.2463,1.1315,0.9708,0.7871,0.6033,0.4196,0.1441,-0.7516,-1.0271,\
            -1.3257,-1.5094,-1.6242,-1.7390,-1.8539,-1.9228,-1.9228,-1.8539,-1.8309,\
            -1.8079,-1.5553,-1.5553,-1.4864,-1.4175,-1.0731,-0.7516,-0.5219,-0.3612,\
            -0.0626,0.1441,0.4196,0.4885,0.6033,0.7641,0.8100,0.8330,0.8789,\
            0.9937,1.0167,1.1086,1.1315,1.2923,1.4530,1.4990,1.4990,1.5219,\
            1.5219,1.4760,1.3382,1.2004,1.0626,1.0397,1.2234,1.2463,1.2004,\
            1.0856,0.6952,-0.2004,-0.3152,-0.4760,-0.6597,-0.7975,-0.9353,-1.0042,\
            -1.0731,-1.0731,-1.0042,-1.0042,-0.8664,-0.6827,0.1670,0.3967,0.6033,\
            0.8100,0.9019,0.9478,0.9937,0.9937,0.8789,0.6722,0.3967,0.2129,\
            0.0981,0.0292,-0.1086,-0.2234,-0.4760,-0.3841,-0.3612,-0.2923,-0.2004,\
            0.1211]

        #Contour values in form of H/Hmax ratio
        zcontour_H_Hmax=[0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,\
            0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,\
            0.2,0.2,0.2,0.2,0.2,0.2,0.3,0.3,0.3,0.3,0.3,\
            0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,\
            0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,\
            0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,\
            0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,\
            0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,\
            0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,\
            0.3,0.3,0.3,0.3,0.3,0.3,0.4,0.4,0.4,0.4,0.4,\
            0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,\
            0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,\
            0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,\
            0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,\
            0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,\
            0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,\
            0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,\
            0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,\
            0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,\
            0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,\
            0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,\
            0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.8,\
            0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,\
            0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,\
            0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,\
            0.8,0.8,0.8,0.8,0.9,0.9,0.9,0.9,0.9,0.9,0.9,\
            0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9]

        xcontour_x_RVmax=np.array(xcontour_x_RVmax)
        ycontour_y_RVmax=np.array(ycontour_y_RVmax)
        zcontour_H_Hmax=np.array(zcontour_H_Hmax)

    #--------------------------------------------------------------------------
    #Contour (isobar) of hurricane waves based on Young (1988)
    #Contour are for hurricane moving toward North
    #Data extracted from Young (2017)

    #Contour (isobar) of hurricane waves based on Young (1988)
    #For hurricane with maximum velocity of Vmax=40 (m/s) and translational velocity of Vt=2.5 (m/s)
    if CalcMethod=='youngvt2':
    
        #x of contour in form of x/RVmax ratio
        xcontour_x_RVmax=[-2.1529,-2.5350,-2.9809,-4.0000,-4.1274,-4.3822,5.6815,5.8089,5.9363,6.1274,\
            6.9554,7.4650,7.5287,7.3376,7.0191,6.2548,5.7452,-4.4459,-4.1274,-3.9363,\
            -3.9363,-3.4904,-3.3631,-3.3631,-3.0446,-2.5350,-2.0255,-1.8344,-1.6433,-1.4522,\
            -1.1338,-0.8153,-0.3057,0.1401,0.5223,0.6497,0.8408,0.9682,1.4140,1.5414,\
            1.6051,1.7325,1.9873,2.2420,2.6879,3.0064,3.3885,3.7070,4.0892,4.2803,\
            4.5350,4.7898,4.9172,4.9809,5.2357,5.5541,5.6178,5.8089,5.9363,6.0637,\
            6.1911,6.2548,6.1911,6.1274,5.7452,5.5541,5.4904,5.2994,4.9809,4.6624,\
            4.5350,4.3439,4.0892,3.7070,3.1338,2.6242,1.7325,1.0318,0.7134,0.4586,\
            -0.0510,-0.4331,-0.8153,-2.1529,-2.4713,-2.8535,-3.1083,-3.2357,-3.3631,-3.6178,\
            -3.8089,-3.8089,-4.0000,-4.3822,-2.0892,-2.5350,-2.5987,-2.5987,-2.4076,-2.2166,\
            -1.9618,-1.5796,-1.0701,-0.6242,-0.2420,0.3949,1.0318,1.9873,2.4968,2.9427,\
            3.4522,3.8981,4.1529,4.1529,4.5350,4.5987,4.6624,4.8535,4.8535,4.6624,\
            3.9618,3.7707,3.3885,2.8153,2.1783,1.6051,0.9682,0.3312,0.0127,-0.1783,\
            0.0127,0.2675,0.5223,0.6497,0.7134,0.7134,0.3312,0.0127,-0.3057,-0.6879,\
            -1.0064,-1.1975,-1.5159,-1.2611,-1.0064,-0.6242,-0.2420,0.0764,0.3312,0.5860,\
            1.0318,1.6051,1.8599,3.0701,3.2611,3.4522,3.5159,3.5159,3.2611,2.7516,\
            2.1783,1.7325,1.4777,1.0955,0.7134,0.6497,0.7134,0.8408,0.8408,0.6497,\
            0.3312,0.0764,-0.3057,-0.9427,-0.1146,-0.1146,0.0127,0.3949,0.6497,0.8408,\
            1.1592,1.4777,1.7325,2.1783,1.8599,1.5414,1.3503,1.0955,1.0318,0.9682,\
            0.7771,0.6497,0.3312]
    
        #y of contour in form of y/RVmax ratio
        ycontour_y_RVmax=[-4.9343,-4.8175,-4.5839,-3.4745,-2.9489,-2.6569,5.9270,5.3431,4.8175,4.2336,\
            3.0657,2.7153,-2.2482,-2.5401,-2.9489,-4.3504,-4.9927,1.0219,1.1971,1.3723,\
            1.4891,1.8978,2.0730,2.1898,3.5328,3.7080,3.8832,4.1752,4.7591,5.0511,\
            5.1679,5.1679,4.9927,4.8175,4.7007,4.8175,5.2847,5.5766,5.9854,5.9854,\
            5.8686,5.6934,5.4599,5.2847,5.0511,4.8175,4.5255,4.3504,4.1752,4.1168,\
            4.1168,3.7080,3.4161,3.0073,2.8905,2.7153,2.4234,2.0146,1.5474,0.9635,\
            0.3212,-0.3796,-0.7883,-0.9635,-2.3066,-2.6569,-2.9489,-3.1241,-3.2993,-3.4745,\
            -3.6496,-4.0000,-4.2920,-4.5255,-4.7007,-4.7591,-4.8759,-4.8759,-4.7007,-4.5839,\
            -4.4672,-4.3504,-3.7664,-2.8321,-2.3650,-2.1314,-1.8978,-1.5474,-1.2555,-0.7883,\
            -0.4380,-0.1460,0.0292,-0.0876,-0.0876,0.3796,0.8467,1.0803,1.4891,2.0146,\
            2.4234,2.8321,3.4161,3.7080,3.7664,3.6496,4.0584,3.6496,3.4745,3.5328,\
            3.1825,2.7153,2.3650,2.1314,1.6642,1.4307,1.0219,-0.4380,-0.6715,-1.1971,\
            -2.1314,-2.4234,-2.8321,-3.1241,-3.3577,-3.5328,-3.1241,-2.5985,-2.0146,-1.4891,\
            -1.3723,-1.0219,-0.7299,-0.4964,-0.2628,0.0292,0.4380,0.6715,0.6131,0.3796,\
            0.0876,0.0292,-0.0876,1.5474,1.9562,2.4234,2.6569,2.7737,2.6569,2.6569,\
            3.0073,2.8321,2.5401,1.7810,1.2555,0.8467,0.3212,-0.3796,-0.9635,-1.4891,\
            -1.8394,-1.9562,-1.9562,-1.7810,-1.4891,-1.1387,-0.7883,-0.4380,-0.0292,0.3212,\
            0.6131,0.7883,0.8467,1.0219,0.9635,1.1387,1.3139,1.4891,1.6642,1.7810,\
            1.7810,1.6058,1.3139,-0.0292,-0.3212,-0.5547,-0.4964,-0.4964,-0.3796,0.0292,\
            0.3796,0.6715,0.8467]
    
        #Contour values in form of H/Hmax ratio
        zcontour_H_Hmax=[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.7,0.7,0.7,0.7,0.7,\
            0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,\
            0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,\
            0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,\
            0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,\
            0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,\
            0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,\
            0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.9,0.9,\
            0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,\
            0.9,0.9,0.9,0.9,0.9,0.9]
    
        xcontour_x_RVmax=np.array(xcontour_x_RVmax)
        ycontour_y_RVmax=np.array(ycontour_y_RVmax)
        zcontour_H_Hmax=np.array(zcontour_H_Hmax)

    #--------------------------------------------------------------------------
    #Contour (isobar) of hurricane waves based on Young (1988)
    #Contour are for hurricane moving toward North
    #Data extracted from Young (2017)

    #Contour (isobar) of hurricane waves based on Young (1988)
    #For hurricane with maximum velocity of Vmax=40 (m/s) and translational velocity of Vt=5 (m/s)
    if CalcMethod=='youngvt5':
    
        #x of contour in form of x/RVmax ratio
        xcontour_x_RVmax=[-1.0270,-1.0257,-0.9710,-0.9693,-1.1828,-1.1811,-1.1788,-1.5011,-1.7697,-3.2254,\
            -3.8711,-4.0336,-4.3602,-4.4670,4.5066,4.7756,5.2604,5.4766,5.6921,5.9611,\
            7.1996,7.5226,2.0000,1.5686,1.0290,0.7057,0.5446,-0.0413,-0.0400,-0.1998,\
            -0.1988,0.0745,0.1863,0.3512,0.5157,0.5174,0.2488,-0.1279,-0.5056,-0.7772,\
            -0.9404,-1.5360,-1.8613,-1.9705,-2.2958,-2.5123,-2.7282,-3.0512,-3.1580,-3.2097,\
            -3.4219,-3.6351,-3.9557,-4.1695,-4.2219,-4.2192,-4.1628,-4.0520,-3.9429,-3.7808,\
            -3.5652,-3.4031,-3.2382,-3.1828,-3.1795,-3.1761,-3.2281,-3.2812,-3.2802,-2.7923,\
            -2.2518,-1.8741,-1.3896,-0.8508,-0.3663,0.0107,0.3885,0.5506,0.7692,0.9334,\
            1.0435,1.3151,1.7499,2.9391,3.0448,3.3661,4.0635,4.6033,5.0344,5.5236,\
            5.7922,6.1131,6.2743,6.5405,6.7534,6.8051,6.9105,7.0687,7.5034,7.5014,\
            7.3379,7.3362,7.3890,7.4413,7.4940,7.4830,7.3191,7.2090,7.1523,7.0959,\
            6.9864,6.0057,5.6257,5.2456,4.7027,-1.4786,-2.0171,-2.2307,-2.3905,-2.3895,\
            -2.3881,-2.2783,-2.0524,-2.1031,-2.1004,-1.9365,-1.6659,-1.3420,-0.8571,-0.2649,\
            0.4898,0.7071,1.0351,1.0898,1.3614,2.0124,2.2286,2.6060,2.9266,3.5189,\
            4.0561,4.3240,4.6453,5.0210,5.3470,5.5656,5.8348,6.0507,6.2095,6.2609,\
            6.2548,6.3065,6.2511,6.2484,6.1924,6.0269,5.2614,5.0428,4.7699,4.2813,\
            4.1189,3.8472,3.5212,3.2503,2.8175,2.3314,1.8999,1.5229,1.1479,0.8800,\
            0.6661,0.6691,0.4563,0.4573,0.5137,0.6235,0.6268,0.5738,0.4133,0.1440,\
            -0.0715,-0.5043,-0.8840,-0.3408,-0.5029,-0.9344,-1.2036,-1.3101,-1.3080,-1.1972,\
            -1.0865,-0.9750,-0.7034,-0.4321,-0.2699,0.1078,0.4311,0.5933,0.8108,1.0287,\
            1.4615,1.7868,2.0581,2.4895,2.7571,2.9176,3.1328,3.4568,3.6723,3.9946,\
            4.1551,4.3676,4.4724,4.6302,4.6826,4.2951,4.0772,3.7512,3.0999,2.4502,\
            2.0178,1.5316,1.3708,1.2100,0.9951,0.8894,0.8931,0.8407,0.7880,0.7349,\
            0.6822,0.5755,0.3606,0.1994,0.0383,-0.1239,-0.3918,-0.3908,-0.2263,-0.0624,\
            0.2085,0.2085,0.5328,0.7510,1.0223,1.2916,1.5605,1.8321,2.0494,2.2642,\
            2.5318,3.3346,3.4407,3.4366,3.2721,3.1079,2.5654,2.2951,1.9171,1.5394,\
            1.4323,1.0559,1.0029,0.9498,0.8968,0.8440,0.7910,0.6305,0.4160,0.1474,\
            -0.0141,-0.2300,0.4190,0.6906,0.9622,1.1784,1.5014,1.7173,1.9335,2.2552,\
            2.4150,2.5170,1.3816,1.0586,0.9522,0.9532,0.7927,0.6332,0.5802]
    
        #y of contour in form of y/RVmax ratio
        ycontour_y_RVmax=[-4.9938,-4.7950,-4.6957,-4.4472,-4.0497,-3.8012,-3.4534,-3.1553,-2.9068,-2.3602,\
            -1.9130,-1.9627,-2.3106,-2.1118,5.9876,5.7888,5.5404,5.5404,5.4410,5.2422,\
            4.5466,4.3478,-4.9938,-4.8447,-4.6957,-4.5466,-4.3975,-3.1056,-2.9068,-2.5590,\
            -2.4099,-1.9627,-1.4161,-1.0186,-0.6708,-0.4224,-0.1739,0.0745,0.1739,-0.0248,\
            -0.1739,-0.3230,-0.4720,-0.6211,-0.7702,-0.8199,-0.7702,-0.5714,-0.3727,-0.0248,\
            0.5714,1.0186,1.5652,1.9130,2.1615,2.5590,2.9068,3.3043,3.4534,3.4534,\
            3.3540,3.3540,3.7516,3.9503,4.4472,4.9441,5.2422,5.3913,5.5404,5.7391,\
            5.7391,5.6398,5.3416,5.0932,4.7950,4.5963,4.4969,4.4969,4.8447,5.1429,\
            5.4410,5.6398,5.9876,5.9876,5.6398,5.1925,4.3975,4.2981,4.0994,4.4969,\
            4.2484,3.7516,3.6025,3.0062,2.5093,2.1615,1.7640,1.1677,1.5155,1.2174,\
            1.0186,0.7702,0.5714,0.3230,0.1242,-1.5155,-1.7640,-2.0621,-2.4596,-2.8075,\
            -3.0062,-4.1491,-4.3975,-4.6460,-4.9938,0.1739,0.4720,0.8696,1.2174,1.3665,\
            1.5652,1.8137,3.2547,3.7516,4.1491,4.3975,4.4472,4.3975,4.1491,3.8012,\
            3.5031,3.6522,4.1988,4.2981,4.4969,4.8447,4.8447,4.6957,4.1491,3.8012,\
            3.3043,2.9565,2.5093,2.1118,2.3602,2.7081,2.5590,2.5093,2.0124,1.6149,\
            0.7205,0.3727,0.1739,-0.2236,-0.5217,-1.0186,-2.3106,-2.6584,-3.0559,-3.3540,\
            -3.4037,-3.6025,-3.8509,-3.9503,-4.0000,-3.9503,-3.8012,-3.6025,-3.1056,-2.7578,\
            -2.4099,-1.9627,-1.4658,-1.3168,-0.9689,-0.7205,-0.2236,-0.0745,0.1739,0.3230,\
            0.4224,0.3727,0.1739,0.5714,0.5714,0.7205,0.8696,1.1180,1.4161,1.8137,\
            2.2112,2.7081,2.9068,3.0559,3.0559,2.9565,2.8075,2.8075,3.0062,3.2547,\
            3.3043,3.4534,3.6025,3.4534,3.0559,2.8075,2.6584,2.6087,2.5093,2.2112,\
            1.9627,1.4161,0.9193,0.2733,0.0248,-1.3168,-1.5652,-1.8137,-2.2112,-2.3602,\
            -2.3602,-2.3106,-2.1118,-1.9130,-1.7143,-1.3665,-0.8199,-0.5714,-0.3727,-0.2236,\
            -0.0248,0.1739,0.3727,0.5217,0.6708,0.6708,1.0186,1.1677,1.5155,1.7640,\
            1.8634,1.8634,1.8634,2.1615,2.3106,2.1615,1.9627,2.1615,2.3106,2.1118,\
            1.7143,0.5217,0.2236,-0.3727,-0.7205,-1.0186,-1.3168,-1.3168,-1.2671,-1.1677,\
            -1.0186,-0.7205,-0.5714,-0.4224,-0.2733,-0.0745,0.0745,0.3230,0.5714,0.8199,\
            0.9193,0.9689,1.0186,1.2174,1.4161,1.4161,1.2174,1.1677,1.1677,0.7702,\
            0.4224,-0.4720,-0.5217,-0.3230,-0.0745,0.0745,0.3230,0.7205,0.8696]
    
        #Contour values in form of H/Hmax ratio
        zcontour_H_Hmax=[0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,\
            0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.7,0.7,0.7,0.7,\
            0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,\
            0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,\
            0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,\
            0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,\
            0.7,0.7,0.7,0.7,0.7,0.8,0.8,0.8,0.8,0.8,0.8,\
            0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,\
            0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,\
            0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.9,0.9,0.9,\
            0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,\
            0.9,0.9,0.9]
    
        xcontour_x_RVmax=np.array(xcontour_x_RVmax)
        ycontour_y_RVmax=np.array(ycontour_y_RVmax)
        zcontour_H_Hmax=np.array(zcontour_H_Hmax)

    #--------------------------------------------------------------------------
    #Contour (isobar) of hurricane waves based on Young (1988)
    #Contour are for hurricane moving toward North
    #Data extracted from Young (2017)

    #Contour (isobar) of hurricane waves based on Young (1988)
    #For hurricane with maximum velocity of Vmax=40 (m/s) and translational velocity of Vt=7.5 (m/s)
    if CalcMethod=='youngvt7':
    
        #x of contour in form of x/RVmax ratio
        xcontour_x_RVmax=[-4.4762,-4.2646,-4.1058,-4.0000,-3.8942,-3.7354,-3.6296,-3.5238,-3.3122,-3.1005,\
            -2.6772,-2.2540,-1.3545,-1.0370,-0.7196,0.6561,1.0265,1.3439,1.5556,1.7672,\
            2.0847,2.4550,2.9841,3.5661,3.9365,4.1481,5.4180,5.8942,6.1058,6.4233,\
            6.9524,7.2169,7.4286,0.8677,0.6561,0.6032,0.6032,0.5503,0.5503,0.4974,\
            0.3915,0.3915,0.2328,0.0212,-0.0847,-0.5079,-0.8254,-1.0370,-1.4074,-2.9418,\
            -3.4180,-3.8413,-4.0529,-4.2646,-4.4762,-4.4762,-4.2646,-4.0000,-3.6296,-3.3122,\
            -3.0476,-2.7831,-2.4656,-0.9312,-0.1376,0.2328,0.4974,0.6561,0.7090,0.9735,\
            1.1852,1.3439,1.4497,1.6085,1.8730,2.1376,2.4021,2.7725,3.0899,3.3545,\
            3.7778,4.2011,4.5714,4.8360,6.1587,6.5820,6.9524,7.4286,7.4286,6.9524,\
            1.5026,1.2910,1.1852,1.1323,1.0794,1.0265,1.0265,0.8677,0.8148,0.8148,\
            0.6561,0.1799,0.0212,0.2857,0.3386,0.3386,0.1270,-0.0317,-0.2434,-0.2963,\
            -0.2434,-0.5079,-0.6138,-0.6138,-0.8254,-0.9841,-1.1429,-1.4074,-1.6720,-1.9894,\
            -2.2011,-2.5714,-2.8889,-3.0476,-3.1534,-3.4180,-3.6296,-3.8413,-4.0000,-4.1058,\
            -4.1587,-4.1587,-4.0000,-3.7884,-3.5238,-3.2593,-3.1005,-2.7831,-2.4656,-2.2011,\
            -1.9894,-1.6720,-1.3016,0.5503,0.7619,1.1323,1.5026,1.8730,2.3492,2.8783,\
            3.3545,3.9365,4.5185,4.8889,5.3122,5.8942,6.2116,6.3704,6.6349,6.6349,\
            6.5820,6.5820,6.4233,6.1587,5.8942,5.6296,5.5767,5.3122,4.8889,4.4127,\
            4.0423,3.6720,3.4074,2.9312,2.6667,2.5079,2.3492,1.9788,1.6614,-1.0370,\
            -1.2487,-1.5132,-1.7249,-1.8307,-1.8836,-1.8836,-1.8307,-1.4603,-0.9841,0.5503,\
            0.8148,1.0794,1.2910,1.5026,1.6614,2.0847,2.2963,2.4550,2.7196,3.1429,\
            3.4074,3.5661,4.0423,4.2540,4.4127,4.6772,4.9418,5.2593,5.3651,5.5238,\
            5.6825,5.7884,5.8413,5.8413,5.7884,5.6825,5.3122,4.4127,4.0423,3.7249,\
            3.4603,3.0370,2.3492,2.0317,1.6085,1.3968,1.2381,1.0794,0.9206,0.6561,\
            0.5503,0.5503,0.4974,0.3915,0.1799,0.0212,-0.0847,-0.2434,-0.4021,-0.5608,\
            -0.7725,-0.8783,-0.9312,-0.6667,-0.6667,-0.4021,0.6032,0.9206,1.1323,1.5026,\
            1.6085,1.9259,2.5079,2.7196,3.0370,3.3545,3.5132,3.6720,3.8307,3.9365,\
            4.0423,4.2011,4.3069,4.2011,3.9894,3.5661,3.3016,2.8783,2.5079,1.9788,\
            1.7143,1.4497,1.2910,1.0265,0.7090,0.6561,0.6032,0.4444,0.2857,0.1799,\
            0.0212,-0.2963,-0.4550,-0.6138,0.3386,0.4444,0.6561,1.0265,1.3439,1.6085,\
            1.8201,2.0317,2.4550,2.8783,3.0370,3.1958,3.4074,3.4603,3.5661,3.5661,\
            3.3545,2.6667,1.8730,1.3968,1.0265,0.9206,0.8148,0.7090,0.5503,0.4974,\
            0.9735,0.7619,0.7619,0.8677,0.9735,1.0794,1.2381,1.3968,1.8730]
    
        #y of contour in form of y/RVmax ratio
        ycontour_y_RVmax=[2.4970,3.0303,3.5636,3.7091,3.7576,3.7576,3.7091,3.7091,3.7576,3.8061,\
            3.9030,3.9030,3.9030,3.9515,4.1455,4.7758,5.1636,5.3576,5.4061,5.3576,\
            5.1636,4.9212,4.7273,4.4848,4.3394,4.1939,3.4667,3.2727,3.1758,2.9818,\
            2.8364,2.5939,2.4000,-4.9697,-4.6303,-4.4848,-4.2424,-3.8061,-3.3212,-2.7879,\
            -2.6424,-2.0121,-1.6727,-1.4303,-1.4303,-1.6242,-2.1576,-2.4000,-2.4970,-2.4970,\
            -2.6909,-2.8848,-2.8848,-2.6424,-2.3030,0.1212,0.6545,1.0909,1.7212,2.2061,\
            2.5939,2.7394,2.8848,2.9818,2.9818,3.0788,3.1758,3.4182,3.5636,3.9515,\
            4.2424,4.3879,4.4364,4.4364,4.2909,4.0970,3.9030,3.7576,3.6121,3.4182,\
            3.1758,2.8848,2.6909,2.5455,1.8667,1.5758,1.2364,0.6061,-4.3394,-4.9697,\
            -4.8242,-4.5818,-4.4364,-4.2424,-4.0485,-3.4667,-2.8364,-2.4970,-1.9636,-1.5758,\
            -1.3818,-0.9939,-0.7030,-0.4606,-0.3152,-0.1212,0.0727,0.0727,-0.1697,-0.4121,\
            -0.6061,-0.6061,-0.8000,-0.8970,-1.1879,-1.3333,-1.3818,-1.3818,-1.2848,-1.0909,\
            -1.0424,-1.0424,-1.0909,-1.1879,-1.3333,-1.5758,-1.8667,-2.1091,-2.2061,-2.0606,\
            -1.9636,-1.0424,-0.6545,-0.5091,-0.1212,0.6061,0.9455,1.3818,1.7697,2.1091,\
            2.2545,2.3030,2.3515,2.4485,2.8364,3.3212,3.7091,3.5152,3.2727,2.9818,\
            2.6424,2.3515,1.9152,1.7212,1.4303,0.9939,0.6061,0.3152,-0.9939,-1.5273,\
            -1.9152,-2.0606,-2.3515,-2.7394,-3.1273,-3.4667,-3.5636,-3.8061,-4.0485,-4.2424,\
            -4.3394,-4.3394,-4.2909,-4.1455,-4.0485,-4.0000,-4.0485,-4.2909,-4.4848,-0.2182,\
            -0.2182,-0.1697,-0.0242,0.1212,0.3636,0.8000,1.0424,1.4788,1.8182,1.9152,\
            2.3030,2.6424,2.8364,2.9333,2.9333,2.7879,2.7394,2.6909,2.4970,2.2061,\
            2.0121,1.9152,1.6242,1.2848,0.9939,0.7515,0.6545,0.4606,0.4606,0.3636,\
            0.1697,-0.1697,-0.4121,-0.5576,-0.8000,-1.0424,-1.3818,-2.4485,-2.6909,-2.7879,\
            -2.7879,-2.6909,-2.4485,-2.3515,-2.1576,-1.9636,-1.6242,-1.3818,-1.2364,-1.0424,\
            -0.8000,-0.1212,0.0242,0.1697,0.4121,0.5091,0.5091,0.4606,0.3152,0.2182,\
            0.0727,-0.0727,-0.1697,0.9939,1.0909,1.2848,1.5273,1.8667,2.0606,2.2545,\
            2.2545,2.1576,2.1091,1.9636,1.8182,1.5273,1.3818,1.1879,1.0909,0.8485,\
            0.6061,0.2182,-1.1879,-1.4788,-1.6727,-1.9152,-1.9152,-1.8182,-1.7212,-1.5758,\
            -1.5758,-1.4303,-1.2364,-0.9455,-0.4121,-0.0727,0.1697,0.3636,0.5091,0.6545,\
            0.7515,0.8485,0.8485,0.8970,0.9939,1.0424,1.1879,1.4788,1.5273,1.5273,\
            1.3333,1.1879,1.0424,1.2364,1.2364,1.0424,0.5576,0.3152,0.1212,-0.3152,\
            -0.6061,-1.0424,-1.0909,-0.8970,-0.6545,-0.5576,-0.1697,0.1697,0.5091,0.8000,\
            -0.2182,0.3636,0.5576,0.8485,0.9939,0.9939,0.8485,0.7030,-0.0727]
    
        #Contour values in form of H/Hmax ratio
        zcontour_H_Hmax=[0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,\
            0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,\
            0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,\
            0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,\
            0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,\
            0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,\
            0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,\
            0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,\
            0.4,0.4,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,\
            0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,\
            0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,\
            0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.8,0.8,\
            0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,\
            0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,\
            0.8,0.8,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9]
    
        xcontour_x_RVmax=np.array(xcontour_x_RVmax)
        ycontour_y_RVmax=np.array(ycontour_y_RVmax)
        zcontour_H_Hmax=np.array(zcontour_H_Hmax)

    #--------------------------------------------------------------------------
    #Contour (isobar) of hurricane waves based on Young (1988)
    #Contour are for hurricane moving toward North
    #Data extracted from Young (2017)

    #Contour (isobar) of hurricane waves based on Young (1988)
    #For hurricane with maximum velocity of Vmax=40 (m/s) and translational velocity of Vt=10 (m/s)
    if CalcMethod=='youngvt10':
    
        #x of contour in form of x/RVmax ratio
        xcontour_x_RVmax=[-4.4737,-4.2632,-4.0526,-3.8421,-3.6842,-3.3684,-3.3684,-3.4737,-3.4737,-3.2632,-3.2105,\
            -3.1579,-3.0526,-2.9474,-2.1053,-1.6842,-1.2105,-0.8947,-0.4737,0.0000,0.6316,1.1579,\
            1.7895,2.0526,2.3684,2.9474,3.6316,5.0526,5.4737,5.8947,6.3684,6.7895,7.0526,\
            7.2105,7.4211,1.1053,0.7368,0.5263,0.5263,0.6316,0.6842,0.5789,0.5263,0.3684,\
            0.2632,0.1053,0.1053,0.4737,0.4737,0.2632,0.0000,-0.9474,-1.2105,-1.5263,-1.8421,\
            -2.5789,-3.3684,-3.6842,-3.7895,-4.0000,-4.1053,-4.2105,-4.2632,-4.2632,-4.2105,-4.1579,\
            -4.0526,-3.8947,-3.6316,-3.3684,-3.1579,-2.7895,-2.5263,-2.3684,-2.1579,-2.1053,-1.8947,\
            -1.7368,-1.4211,-0.9474,-0.7895,-0.5263,-0.1053,0.4737,0.8421,1.2105,1.8421,2.4737,\
            3.0526,3.5789,4.1053,5.2632,5.6316,6.1053,6.4737,6.8421,6.9474,7.1053,7.2105,\
            7.4211,2.5263,2.2632,2.1053,1.0000,0.8421,0.7368,0.6316,0.6316,0.5789,0.5789,\
            0.4737,0.3158,0.0000,-0.2105,-0.3684,-0.5789,-0.8947,-1.1053,-1.4737,-1.8421,-2.1579,\
            -2.4211,-2.7895,-3.0000,-3.1053,-3.2632,-3.3158,-3.3158,-3.2632,-3.1053,-2.9474,-2.6316,\
            -2.3684,-2.0526,-1.6842,-1.0000,-0.6842,-0.4737,-0.2105,0.0000,0.3158,0.8947,1.1053,\
            1.3684,1.6316,2.0526,2.2632,2.6316,2.9474,3.3158,3.6842,3.8421,4.0000,4.2632,\
            4.6316,4.9474,5.8947,6.1579,6.3158,6.7368,6.8947,7.0000,7.0526,7.1053,7.1053,\
            7.0000,6.7895,6.6842,6.5789,6.3158,6.0526,0.0000,-0.0526,-0.0526,-0.3158,-0.3158,\
            -0.1053,0.1053,0.1579,0.1579,-0.0526,-0.3158,-0.3158,-0.1579,0.1053,0.4737,0.7368,\
            0.9474,1.5263,1.8421,2.2632,2.6316,3.0526,3.4737,3.7368,4.1053,4.4211,4.8421,\
            5.0526,5.2105,5.4211,6.0526,6.0526,6.0000,5.8421,5.6842,5.4737,5.2105,4.8947,\
            4.4737,4.1053,3.6842,3.4211,2.9474,2.6316,2.4211,2.0526,1.7368,1.5263,1.3158,\
            1.0000,0.8947,0.7895,0.7368,0.6842,0.6316,0.4737,0.3158,1.0000,0.9474,0.9474,\
            1.0000,1.0000,0.8947,0.6842,0.6316,0.5789,0.5789,0.6842,0.6842,0.4737,0.4211,\
            0.4211,0.4737,0.6842,1.2105,1.5263,1.7368,2.1053,2.4737,2.8421,3.2632,3.6842,\
            4.1053,4.4211,4.8421,5.0000,5.0000,4.5789,3.7368,3.0526,2.7368,2.4211,1.7368,\
            1.4211,1.9474,1.6316,1.2632,1.0000,0.8421,0.7368,0.7368,0.7895,0.7895,0.7368,\
            0.7368,0.8421,0.9474,1.3158,1.6842,2.1053,2.6316,3.3158,3.5789,4.0000,4.1579,\
            4.1579,4.1579,4.1053,3.8947,3.6842,3.3158,2.0000,1.6842,1.5263,1.1579,0.8947,\
            0.8947,0.8947,0.9474,1.0000,1.1053,1.4737,1.7895,2.1053,2.7368,2.9474,2.8947]
    
        #y of contour in form of y/RVmax ratio
        ycontour_y_RVmax=[-1.5460,-1.2025,-0.8589,-0.6626,-0.5153,-0.1227,0.1227,0.4172,0.5644,1.0061,\
            1.3497,1.6442,1.9877,2.1840,3.2638,3.5583,3.9018,4.1472,4.3926,4.5890,\
            4.7362,4.8834,4.9816,4.9816,4.9325,4.7853,4.6380,4.0000,3.7055,3.3129,\
            2.9202,2.4785,2.1840,1.7914,1.3988,-4.9325,-4.6380,-4.5399,-4.3436,-4.1963,\
            -3.8528,-3.5092,-2.8221,-2.9693,-3.4601,-3.8037,-4.0491,-4.3926,-4.5399,-4.7853,\
            -4.8834,-4.9325,-4.5890,-4.2454,-4.2454,-4.4908,-4.4908,-4.3436,-4.1963,-4.0491,\
            -4.0491,-3.8037,-3.6074,-3.2638,-3.0675,-2.7239,-2.4294,-2.2331,-1.8896,-1.6442,\
            -1.3497,-0.9080,-0.4663,-0.2209,0.2699,0.5153,0.9080,1.3497,1.9387,2.5276,\
            2.6258,2.8712,3.0184,3.2638,3.4110,3.5092,3.5092,3.4110,3.3129,3.0675,\
            2.7730,1.8405,1.4479,0.8589,0.4663,-0.0736,-0.3681,-0.6626,-1.0061,-1.3988,\
            -4.9816,-4.6380,-4.4417,-3.4110,-3.0675,-2.7730,-2.3804,-2.1350,-1.8896,-1.6933,\
            -1.6442,-1.7423,-1.9877,-2.2331,-2.4294,-2.9693,-3.1166,-3.1656,-3.2147,-3.3620,\
            -3.5092,-3.6074,-3.6074,-3.5583,-3.5092,-3.1656,-3.0675,-2.9202,-2.7239,-2.3804,\
            -2.2331,-1.7914,-1.3988,-0.9571,-0.4172,0.8098,1.3006,1.5460,1.7914,2.0368,\
            2.1840,2.4294,2.4785,2.4785,2.4294,2.3804,2.3313,2.1840,2.0368,1.8405,\
            1.5460,1.3988,1.3006,1.0552,0.7607,0.5153,-0.5644,-0.9080,-1.2025,-2.0368,\
            -2.4785,-2.8712,-3.1166,-3.4110,-3.6074,-3.8528,-4.2454,-4.3926,-4.5399,-4.6871,\
            -4.8834,-0.9080,-0.8098,-0.4663,-0.0736,0.0736,0.2209,0.3681,0.4663,0.6135,\
            0.8098,0.9571,1.0552,1.2025,1.3497,1.4479,1.4969,1.5951,1.6442,1.5460,\
            1.4479,1.2515,1.0061,0.7117,0.4663,0.2699,-0.0736,-0.4663,-0.6626,-0.9080,\
            -1.1534,-2.3804,-2.9202,-3.0675,-3.3129,-3.5583,-3.7055,-3.7546,-3.8037,-3.9018,\
            -4.0000,-4.0982,-4.0982,-4.0000,-3.8037,-3.5583,-3.2147,-3.0184,-2.7730,-2.7730,\
            -2.8221,-2.5767,-2.2331,-1.4969,-1.3006,-1.1043,-0.9571,-0.9080,-2.1350,-2.0368,\
            -1.8896,-1.5951,-1.4969,-1.2025,-0.9080,-0.7117,-0.6135,-0.4172,-0.1227,0.0736,\
            0.4663,0.7117,0.9080,1.0552,1.1043,1.1043,1.1043,1.0061,0.9080,0.7607,\
            0.5153,0.2209,-0.2209,-0.6626,-1.0061,-1.3988,-1.7914,-2.0368,-2.7730,-2.9202,\
            -2.9202,-2.9202,-2.7730,-2.2331,-2.0368,-1.7914,-1.5951,-1.3497,-1.0552,-0.8589,\
            -0.6135,-0.3681,-0.1227,0.1227,0.3190,0.5644,0.7117,0.8098,0.8098,0.7117,\
            0.5153,0.2209,-0.5153,-0.9080,-1.3006,-1.4479,-1.6933,-1.9877,-2.1350,-2.2331,\
            -2.3313,-2.3313,-1.1534,-1.0552,-1.0552,-0.8589,-0.6135,-0.4172,-0.1227,0.1227,\
            0.4663,0.5644,0.4663,0.3190,0.1227,-1.1534,-1.3988,-1.5460]
    
        #Contour values in form of H/Hmax ratio
        zcontour_H_Hmax=[0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,\
            0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,\
            0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,\
            0.3,0.3,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,\
            0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,\
            0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,\
            0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,\
            0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,\
            0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,\
            0.4,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,\
            0.5,0.5,0.5,0.5,0.5,0.5,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,\
            0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.7,0.7,0.7,\
            0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,\
            0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,\
            0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,\
            0.7,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,\
            0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,\
            0.8,0.8,0.8,0.8,0.8,0.8,0.9,0.9,0.9,0.9,0.9,\
            0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9]
    
        xcontour_x_RVmax=np.array(xcontour_x_RVmax)
        ycontour_y_RVmax=np.array(ycontour_y_RVmax)
        zcontour_H_Hmax=np.array(zcontour_H_Hmax)

    #--------------------------------------------------------------------------
    #Converting azimuth (bearing) to trigonometric direction

    #Converting azimuth (bearing) to trigonometric direction
    VtTrigdir=-VtAzmdir+90
    
    #Add 360 to all numbers to have them all positive
    #Use mod(360) to take care of the ones larger than 360
    VtTrigdir=((VtTrigdir+360)%360) 
    
    #Converting azimuth (bearing) to trigonometric direction
    if distCalcMethod=='gc':
    
        #Converting azimuth (bearing) to trigonometric direction
        thetagrid=-thetagrid+90
    
        #Add 360 to all numbers to have them all positive
        #Use mod(360) to take care of the ones larger than 360
        thetagrid=((thetagrid+360)%360) 
    

    #--------------------------------------------------------------------------
    #Calculating new grid location

    #Calculating new grid location as a function of R/RVmax based on 
    #Shore Protection Manual (SPM), U.S. Army Corps of Engineers (1984)
    #Coastal Engineering Manual (CEM), U.S. Army Corps of Engineers (2015)
    if ((CalcMethod=='spm') or (CalcMethod=='cem')):
    
        for i in range(0,len(xCenter),1):
    
            #Calculating R/RVmax for each grid location
            Rgrid_RVmax=Rgrid[:,:,i]/RVmax[i]
    
            #Axis rotation respect to North (90 degree)
            Beta=VtTrigdir[i]-90
    
            #Adjusting direction of each grid point based on moving direction of hurricane
            #SPM contours are for hurricane moving toward North
            thetagridtoVt=thetagrid[:,:,i]-Beta
    
            #Converting polar to Cartesian system with axis rotation equal to Beta
            xgrid_R_RVmax[:,:,i],ygrid_R_RVmax[:,:,i]=pol2cart(np.deg2rad(thetagridtoVt),Rgrid_RVmax)
    
    #Calculating new grid location as a function of x/RVmax based on Young (1988)
    elif ((CalcMethod=='cemnext') or (CalcMethod=='youngvt2') or (CalcMethod=='youngvt5') or (CalcMethod=='youngvt7') or (CalcMethod=='youngvt10')):
    
        for i in range(0,len(xCenter),1):
    
            #Axis rotation respect to North (90 degree)
            Beta=VtTrigdir[i]-90
    
            #Adjusting direction of each grid point based on moving direction of hurricane
            #SPM contours are for hurricane moving toward North
            thetagridtoVt=thetagrid[:,:,i]-Beta
    
            #Converting polar to Cartesian system with axis rotation equal to Beta
            xgridprime,ygridprime=pol2cart(np.deg2rad(thetagridtoVt),Rgrid[:,:,i])
    
            #Calculating x/RVmax and y/RVmax for each grid location
            xgrid_RVmax[:,:,i]=xgridprime/RVmax[i]
            ygrid_RVmax[:,:,i]=ygridprime/RVmax[i]
    

    #--------------------------------------------------------------------------
    #Calculating hurricane wave height contour

    #Calculating hurricane wave height contour based on 
    #Shore Protection Manual (SPM), U.S. Army Corps of Engineers (1984)
    #Coastal Engineering Manual (CEM), U.S. Army Corps of Engineers (2015)
    if ((CalcMethod=='spm') or (CalcMethod=='cem')):
    
        for i in range(0,len(xCenter),1):
            
            #Hurricane is in northern hemisphere
            if np.mean(yCenter)>=0:
                #In northern hemisphere, velocity vector rotate counter clockwise with right angle respect to radius
    
                #Calculating H/Hmax ratio for each grid point
                H_Hmax_grid=sp.interpolate.griddata((xcontour_R_RVmax,ycontour_R_RVmax),zcontour_H_Hmax,(xgrid_R_RVmax[:,:,i],ygrid_R_RVmax[:,:,i]))
    
            #Hurricane is in southern hemisphere
            elif np.mean(yCenter)<0:
                #In southern hemisphere, velocity vector rotate clockwise with right angle respect to radius
    
                #Calculating H/Hmax ratio for each grid point
                H_Hmax_grid=sp.interpolate.griddata((-1.*xcontour_R_RVmax,ycontour_R_RVmax),zcontour_H_Hmax,(xgrid_R_RVmax[:,:,i],ygrid_R_RVmax[:,:,i]))
    
    
            H_Hmax_grid[np.isnan(H_Hmax_grid)==1]=0
            
            #Calcuating wave height at each grid point
            Hgrid[:,:,i]=Hmax[i]*H_Hmax_grid
    
    
    #Calculating new grid location as a function of x/RVmax based on Young (1988)
    elif ((CalcMethod=='cemnext') or (CalcMethod=='youngvt2') or (CalcMethod=='youngvt5') or (CalcMethod=='youngvt7') or (CalcMethod=='youngvt10')):
    
        for i in range(0,len(xCenter),1):
            
            #Hurricane is in northern hemisphere
            if np.mean(yCenter)>=0:
                #In northern hemisphere, velocity vector rotate counter clockwise with right angle respect to radius
    
                #Calculating H/Hmax ratio for each grid point
                H_Hmax_grid=sp.interpolate.griddata((xcontour_x_RVmax,ycontour_y_RVmax),zcontour_H_Hmax,(xgrid_RVmax[:,:,i],ygrid_RVmax[:,:,i]))
    
            #Hurricane is in southern hemisphere
            elif np.mean(yCenter)<0:
                #In southern hemisphere, velocity vector rotate clockwise with right angle respect to radius
    
                #Calculating H/Hmax ratio for each grid point
                H_Hmax_grid=sp.interpolate.griddata((-1*xcontour_x_RVmax,ycontour_y_RVmax),zcontour_H_Hmax,(xgrid_RVmax[:,:,i],ygrid_RVmax[:,:,i]))
    
    
            H_Hmax_grid[np.isnan(H_Hmax_grid)==1]=0
            
            #Calcuating wave height at each grid point
            Hgrid[:,:,i]=Hmax[i]*H_Hmax_grid
    
    Hgrid[Hgrid<0]=0
    Hgrid[Hgrid>1]=1
    Hgrid[np.isnan(Hgrid)==True]=0

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='imagesc': #2 dimensional plot using imagesc or imshow
        
        plt.imshow(Hgrid[:,:,-1],extent=[np.nanmin(xgrid),np.nanmax(xgrid),np.nanmin(ygrid),np.nanmax(ygrid)],interpolation=None,cmap=plt.get_cmap(),aspect='auto',origin='lower')
        plt.plot([xCenter[-1],xCenter[-1]],[np.min(ygrid),np.max(ygrid)],color=[1,1,1])
        plt.plot([np.min(xgrid),np.max(xgrid)],[yCenter[-1],yCenter[-1]],color=[1,1,1])
    
    elif dispout=='pcolor': #2 dimensional plot using pcolor
    
        plt.pcolormesh(xgrid,ygrid,Hgrid[:,:,-1],cmap=plt.get_cmap())    
        plt.plot([xCenter[-1],xCenter[-1]],[np.min(ygrid),np.max(ygrid)],color=[1,1,1])
        plt.plot([np.min(xgrid),np.max(xgrid)],[yCenter[-1],yCenter[-1]],color=[1,1,1])
    
    elif dispout=='contour': #2 dimensional contour plot
    
        plt.contourf(xgrid,ygrid,Hgrid[:,:,-1],cmap=plt.get_cmap())
        plt.plot([xCenter[-1],xCenter[-1]],[np.min(ygrid),np.max(ygrid)],color=[1,1,1])
        plt.plot([np.min(xgrid),np.max(xgrid)],[yCenter[-1],yCenter[-1]],color=[1,1,1])
    
    
    #Setting plot properties
    if ((dispout=='imagesc') or (dispout=='pcolor') or (dispout=='contour')):
    
        #Setting axis limits
        #plt.xlim([np.min(xgrid),np.max(xgrid)])
        #plt.ylim([np.min(ygrid),np.max(ygrid)])
    
        #Plotting colorbar
        plt.colorbar()
    
        #Setting label
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Wave Height')
        

    #--------------------------------------------------------------------------
    #Converting output array with size of (Mgrid,Ngrid,1) to array with size of (Mgrid,Ngrid)

    if len(xCenter)==1:
        Hgrid=np.reshape(Hgrid,(Mgrid,Ngrid)) #Reshaping array


    #--------------------------------------------------------------------------
    #Outputs
    return Hgrid

    #--------------------------------------------------------------------------
