def fullydevwave(CalcMethod='pierson', windvel=0):
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

    scientimate.fullydevwave
    ========================

    .. code:: python

        EhatFullyDev, fphatFullyDev = scientimate.fullydevwave(CalcMethod='pierson', windvel=0)

    Description
    -----------

    Calculate a fully developed condition for wind wave growth

    Inputs
    ------

    CalcMethod='pierson'
        | Parametric wave model to be used for calculation 
        | 'pierson': Use method by Pierson and Moskowitz (1964)
        | 'spm': Use method by Shore Protection Manual (SPM), U.S. Army Corps of Engineers (1984)
        | 'cem': Use method by Coastal Engineering Manual (CEM),
        |     U.S. Army Corps of Engineers (2015)
    windvel=0
        | Wind velocity in (m/s)
        | Wind velocity should be measured (or represents velocity) at 10 m above surface
        | For 'cem' and 'spm' methods, wind velocity should be converted to duration of sustained wind by using gust factor

    Outputs
    -------

    EhatFullyDev
        Dimensionless wave energy for a fully developed condition, Ehat=g^2*m0/U10^4
    fphatFullyDev
        | Dimensionless peak wave frequency for a fully developed condition, fphat=fp*U10/9.81 
        | Note, g=9.81: gravitational acceleration
        |     U10: wind velocity
        |     fp: peak wave frequency
        |     m0: Zero-moment of water surface elevation power spectral density

    Examples
    --------

    .. code:: python

        import scientimate as sm
        EhatFullyDev,fphatFullyDev=sm.fullydevwave('pierson')

    References
    ----------

    Department of the Army, Waterways Experiment Station, Corps of Engineers, 
    and Coastal Engineering Research Center (1984), 
    Shore Protection Manual, Washington, 
    D.C., vol. 1, 4th ed., 532 pp.

    Pierson, W. J., & Moskowitz, L. (1964). 
    A proposed spectral form for fully developed wind seas based on the similarity theory of SA Kitaigorodskii. 
    Journal of geophysical research, 69(24), 5181-5190.

    U.S. Army Corps of Engineers (2015). 
    Coastal Engineering Manual. 
    Engineer Manual 1110-2-1100, Washington, D.C.: U.S. Army Corps of Engineers.

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

    #--------------------------------------------------------------------------
    #Calculating a fully developed condition

    #Calculating fully developed condition based on Pierson and Moskowitz (1964)
    if CalcMethod=='pierson':
    
        EhatFullyDev=3.64e-3 #Ehat for fully developed condition
        fphatFullyDev=0.133 #fphat for fully developed condition
    
    #Calculating fully developed condition based on Shore Protection Manual (SPM), U.S. Army Corps of Engineers (1984)
    elif CalcMethod=='spm':
    
        UA=0.71*windvel**1.23 #Wind stress factor or adjusted wind velcity
    
        Hm0hatFullyDev=2.433e-1 #Fully Developed Condition from SPM (1984) based on UA, Holthuijsen(2007), Hm0hat*=(gHm0)/(UA^2)=2.433*10^-1
        Hm0FullyDev=Hm0hatFullyDev*(UA**2)/9.81 #Significant wave height for fully Developed Condition
        m0FullyDev=(Hm0FullyDev**2)/16 #Zero-moment of surface elevation power spectral density for fully Developed Condition, Hm0=4*(m0)^0.5
        EhatFullyDev=9.81**2*m0FullyDev/(windvel**4) #Dimensionless wave energy for fully Developed Condition: Ehat=g^2*m0/U10^4
    
        TphatFullyDev=8.134 #Fully Developed Condition from SPM (1984), Holthuijsen(2007), Tphat*=(gTp)/UA=8.134
        TpFullyDev=TphatFullyDev*UA/9.81 #Peak wave period for fully Developed Condition
        fpFullyDev=1/TpFullyDev #Peak wave frequency for fully Developed Condition: fp=1/Tp
        fphatFullyDev=fpFullyDev*windvel/9.81 #Dimensionless peak wave frequency for fully Developed Condition: fphat=fp*U10/g
    
    #Calculating fully developed condition based Coastal Engineering Manual (CEM), U.S. Army Corps of Engineers (2015)
    elif CalcMethod=='cem':
    
        CD=0.001*(1.1+0.035*windvel) #Drag Coefficient
        Ustar=(CD*(windvel**2))**0.5 #Shear Velocity in (m/s)
    
        Hm0hatFullyDev=2.115e2 #Fully developed condition for Hm0hat=g*Hm0/U*^2
        Hm0FullyDev=Hm0hatFullyDev*(Ustar**2)/9.81 #Zero-moment wave height: Hm0hat=g*Hm0/U*^2
        m0FullyDev=(Hm0FullyDev**2)/16 #Zero-moment of surface elevation power spectral density, Hm0=4*(m0)^0.5
        EhatFullyDev=9.81**2*m0FullyDev/(windvel**4) #Dimensionless wave energy: Ehat=g^2*m0/U10^4
    
        TphatFullyDev=2.398e2 #Fully developed condition for Tphat=g*Tp/U*
        TpFullyDev=TphatFullyDev*Ustar/9.81 #Peak wave period: Tphat=g*Tp/U*
        fpFullyDev=1/TpFullyDev #Peak wave frequency: fp=1/Tp
        fphatFullyDev=fpFullyDev*windvel/9.81 #Dimensionless peak wave frequency: fphat=fp*U10/g
    
    #--------------------------------------------------------------------------
    #Outputs
    return EhatFullyDev, fphatFullyDev

    #--------------------------------------------------------------------------
