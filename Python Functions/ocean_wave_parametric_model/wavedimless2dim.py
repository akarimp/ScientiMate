def wavedimless2dim(windvel, dimlessValue, ValueType='waveheight'):
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

    scientimate.wavedimless2dim
    ===========================

    .. code:: python

        dimValue = scientimate.wavedimless2dim(windvel, dimlessValue, ValueType='waveheight')

    Description
    -----------

    Calculate dimensional numbers from dimensionless numbers

    Inputs
    ------

    windvel
        Wind velocity in (m/s)
    dimlessValue
        Dimensionless value to be converted to dimensional value
    ValueType='waveheight'
        | Type of the dimensionless value 
        | 'fecth': Dimensionless wind fetch (Fetchhat): Fetchhat=g*Fetch/U10^2
        | 'depth':Dimensionless depth (hhat): hhat=g*h/U10^2
        | 'frequency':Dimensionless wave frequency (fhat): fhat=f*U10/9.81
        | 'energy':Dimensionless wave energy (Ehat): Ehat=g^2*m0/U10^4
        | 'waveheight':Dimensionless wave height (Hhat): Hhat=g*H/U10^2
        | 'period':Dimensionless wave period (That): That=g*T/U10
        | 'wavelength':Dimensionless wave length (Lhat): Lhat=g*L/U10^2
        | 'wavenumber':Dimensionless wave number (khat): khat=k*U10^2/g
        | 'time':Dimensionless time (that): that=g*t/U10
        | 'length':Dimensionless length (xhat): xhat=g*x/U10^2

    Outputs
    -------

    dimValue
        | Dimensional value calculated from dimensionless value
        | For ValueType='fetch': dimValue is wind fetch in (m)
        | For ValueType='depth': dimValue is water depth in (m)
        | For ValueType='frequency': dimValue is wave frequency in (Hz)
        | For ValueType='energy': dimValue is zero-moment of surface elevation power spectral density in (m^2)
        | For ValueType='waveheight': dimValue is wave height in (m)
        | For ValueType='period': dimValue is wave period in (s)
        | For ValueType='wavelength': dimValue is wave length in (m)
        | For ValueType='wavenumber': dimValue is wave number in (1/m) or (Radian/m)
        | For ValueType='time': dimValue is time in (s)
        | For ValueType='length': dimValue is length in (m)
        | Note, g=9.81: gravitational acceleration
        |     U10: wind velocity
        |     Fetch: Wind fetch in (m)
        |     h: Water depth in (m)
        |     f: Wave frequency in (Hz)
        |     m0: Zero-moment of water surface elevation power spectral density in (m^2)
        |     H: Wave height in (m)
        |     T: Wave period in (s)
        |     L: Wave length in (m)
        |     k: Wave number in (1/m) or (Radian/m)
        |     t: Time in (s)
        |     x: Length in (m)

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        windvel=10*np.random.rand(100)
        Fetchhat=10000*np.random.rand(100)
        dimlessValue=Fetchhat.copy()
        dimValue=sm.wavedimless2dim(windvel,dimlessValue,'fetch')

    References
    ----------

    Sverdrup, H. U., & Munk, W. H. (1947). 
    Wind, sea, and swell: theory of relations for forecasting. 
    U.S. Navy Department, Hydrographic Office, Publication No. 601, 44 pp. 

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
    dimlessValue=type2numpy(dimlessValue)

    #--------------------------------------------------------------------------
    #Calculate dimensional numbers from dimensionless number

    #Calculating dimensionless Fetch
    if ValueType=='fetch':
        Fetchhat=dimlessValue
        Fetch=Fetchhat*(windvel**2)/9.81 #Wind fetch: Fetchhat=g*Fetch/U10^2
        dimValue=Fetch.copy()
    
    #Calculating dimensionless water depth
    elif ValueType=='depth':
        hhat=dimlessValue
        h=hhat*(windvel**2)/9.81 #Depth: hhat=g*h/U10^2
        dimValue=h.copy()
    
    #Calculating dimensionless wave frequency
    elif ValueType=='frequency':
        fhat=dimlessValue
        f=fhat*9.81/(windvel) #Wave frequency, fhat=f*U10/g
        dimValue=f.copy()
    
    #Calculating dimensionless wave energy
    elif ValueType=='energy':
        Ehat=dimlessValue
        m0=Ehat*(windvel**4)/(9.81**2) #Zero-Moment, Ehat=g^2*m0/U10^4
        dimValue=m0.copy()
    
    #Calculating dimensionless wave height
    elif ValueType=='waveheight':
        Hhat=dimlessValue
        H=Hhat*(windvel**2)/9.81 #Wave height: Hhat=g*H/U10^2
        dimValue=H.copy()
    
    #Calculating dimensionless wave period
    elif ValueType=='period':
        That=dimlessValue
        T=That*windvel/9.81 #Wave period: That=g*T/U10
        dimValue=T.copy()
    
    #Calculating dimensionless wave length
    elif ValueType=='wavelength':
        Lhat=dimlessValue
        L=Lhat*(windvel**2)/9.81 #Wave length: Lhat=g*L/U10^2
        dimValue=L.copy()
    
    #Calculating dimensionless wave number
    elif ValueType=='wavenumber':
        khat=dimlessValue
        k=khat*9.81/(windvel**2) #Wave number: khat=k*U10^2/g
        dimValue=k.copy()
    
    #Calculating dimensionless time
    elif ValueType=='time':
        that=dimlessValue
        t=that*windvel/9.81 #Time: that=g*t/U10
        dimValue=t.copy()
    
    #Calculating dimensionless length
    elif ValueType=='length':
        xhat=dimlessValue
        x=xhat*(windvel**2)/9.81 #Length: xhat=g*x/U10^2
        dimValue=x.copy()

    #--------------------------------------------------------------------------
    #Outputs
    return dimValue

    #--------------------------------------------------------------------------
