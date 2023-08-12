def hurricanewindvelmax(Pc, CalcMethod='atkinson', VgToVCoeff=0.8, dispout='no'):
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

    scientimate.hurricanewindvelmax
    ===============================

    .. code:: python

        Vmax, Vgmax = scientimate.hurricanewindvelmax(Pc, CalcMethod='atkinson', VgToVCoeff=0.8, dispout='no')

    Description
    -----------

    | Calculate hurricane maximum wind velocity at the surface level
    | For Holland (1980) method use hurricanewindvelmaxh80 function
    | For Holland (2008) method use hurricanewindvelmaxh08 function

    Inputs
    ------

    Pc
        Hurricane central surface pressure in (Pa)
    Pn=101325
        | Ambient surface pressure (external pressure) in (Pa)
        | Standard atmosphere pressure is 101325 (Pa) 
        | Typical value: Pn=101500 (Pa) for the western North Pacific, Pn= 101000 (Pa) for the North Atlantic
        | (Batke et al., 2014)
    CalcMethod='atkinson'
        | Hurricane velocity calculation method 
        | 'atkinson': Hurricane maximum wind velocity at the surface level is calculated using Atkinson & Holliday (1977)
        | 'dvorak': Hurricane maximum wind velocity at the surface level is calculated Dvorak tabular pressure–wind model (Holland, 2010)
    VgToVCoeff=0.8
        | Coefficient to convert gradient wind velocity to wind velocity at 10 m above surface as: 
        | V=VgToVCoeff*Vg, if VgToVCoeff=1, then V=Vg
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    Vmax
        Maximum hurricane 1-min wind velocity at the surface level (at 10-m height) in (m/s)
    Vgmax
        | Maximum hurricane 1-min wind velocity at the gradient level in (m/s)
        | Vgmax is calculated from Vmax as Vgmax=Vmax/VgToVCoeff
        | Gradient wind velocity can be converted to wind velocity at 10 m above surface by V=VgToVCoeff*Vg
        | VgToVCoeff can be approximated as VgToVCoeff=0.8
        | For detail on converting gradient (mean boundary-layer) wind velocity to velocity 10 m above surface see
        | e.g. Graham & Numm (1959), Young & Vinoth (2013), Phadke et al. (2003), Powell et al. (2003), Valamanesh et al. (2016), Wei et al. (2017)
        | Shore Protection Manual (SPM), U.S. Army Corps of Engineers (1984)
        | Coastal Engineering Manual (CEM), U.S. Army Corps of Engineers (2015)

    Examples
    --------

    .. code:: python

        import scientimate as sm


        #EXAMPLE 1

        Pc=90200 #(Pa)
        Vmax,Vgmax=sm.hurricanewindvelmax(Pc,'atkinson',0.8,'no')


        #EXAMPLE 2

        #Hurricane Katrina centeral pressure (Pa)
        Pc=[100800,100700,100700,100600,100300,100000,99700,99400,98800,98400,98300,98700,\
            97900,96800,95900,95000,94200,94800,94100,93000,90900,90200,90500,91300,\
            92000,92300,92800,94800,96100,97800,98500,99000,99400,99600]

        Vmax,Vgmax=sm.hurricanewindvelmax(Pc,'atkinson',0.8,'yes')

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

    Harper, B. A., & Holland, G. J. (1999, January). 
    An updated parametric model of the tropical cyclone. 
    In Proc. 23rd Conf. Hurricanes and Tropical Meteorology.

    Holland, G. J., Belanger, J. I., & Fritz, A. (2010). 
    A revised model for radial profiles of hurricane winds. 
    Monthly Weather Review, 138(12), 4393-4401.

    Phadke, A. C., Martino, C. D., Cheung, K. F., & Houston, S. H. (2003). 
    Modeling of tropical cyclone winds and waves for emergency management. 
    Ocean Engineering, 30(4), 553-578.

    Powell, M. D., Vickery, P. J., & Reinhold, T. A. (2003). 
    Reduced drag coefficient for high wind speeds in tropical cyclones. 
    Nature, 422(6929), 279.

    U.S. Army Corps of Engineers (2015). 
    Coastal Engineering Manual. 
    Engineer Manual 1110-2-1100, Washington, D.C.: U.S. Army Corps of Engineers.

    Valamanesh, V., Myers, A. T., Arwade, S. R., Hajjar, J. F., Hines, E., & Pang, W. (2016). 
    Wind-wave prediction equations for probabilistic offshore hurricane hazard analysis. 
    Natural Hazards, 83(1), 541-562.

    Wei, K., Arwade, S. R., Myers, A. T., Valamanesh, V., & Pang, W. (2017). 
    Effect of wind and wave directionality on the structural performance of non‐operational offshore wind turbines supported by jackets during hurricanes. 
    Wind Energy, 20(2), 289-303.

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
    
    Pc=type2numpy(Pc)

    #--------------------------------------------------------------------------
    #Calculating hurricane maximum wind velocity at the surface level

    #Calculating hurricane maximum wind velocity at the surface level using Atkinson & Holliday (1977) 
    #It is suitable for the western North Pacific (Holland, 2010)
    if CalcMethod=='atkinson':
    
        Vmax=0.51*6.7*(1010-Pc*1e-2)**0.644 #hurricane maximum 1-min sustained wind velocity at 10m height in (m/s), Pc is in (Pa)
    
    #Calculating hurricane maximum wind velocity at the surface level for Dvorak tabular pressure–wind model (Holland, 2010)
    #It is suitable for the North Atlantic (Holland, 2010)
    elif CalcMethod=='dvorak':
    
        Vmax=3.92*(1015-Pc*1e-2)**0.644 #hurricane maximum 1-min sustained wind velocity at 10m height in (m/s), Pc is in (Pa)
    
    
    #Converting velocity 10 m above surface to gradient (mean boundary-layer) wind velocity
    #e.g. Graham & Numm (1959), Young & Vinoth (2013), Phadke et al. (2003), Powell et al. (2003), Valamanesh et al. (2016), Wei et al. (2017)
    #Shore Protection Manual (SPM), U.S. Army Corps of Engineers (1984)
    #Coastal Engineering Manual (CEM), U.S. Army Corps of Engineers (2015)
    Vgmax=Vmax/VgToVCoeff

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
        
        #Plotting data
        plt.scatter(Pc,Vmax,label='Vmax')
        plt.scatter(Pc,Vgmax,label='Vgmax')
            
        plt.xlabel('Pc (Pa)')
        plt.ylabel('Vmax, Vgmax, (m/s)')
        plt.legend()
        

    #--------------------------------------------------------------------------
    #Output
    return Vmax, Vgmax

    #--------------------------------------------------------------------------
