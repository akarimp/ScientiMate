def surfaceroughness(z, u, delta=None, dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2018-03-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.surfaceroughness
    ============================

    .. code:: python

        ustar, z0, d = scientimate.surfaceroughness(z, u, delta=None, dispout='no')

    Description
    -----------

    Calculate shear velocity and surface roughness from a given velocity profile using Karimpour et al. (2012) method

    Inputs
    ------

    z
        Distance from a surface (elevation, height) in (m)
    u
        Velocity at z in (m/s)
    delta=max(z)
        Boundary layer height in (m)
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    z0
        Surface roughness in (m)
    ustar
        Shear velocity (u*) in (m/s)
    d
        | Zero plane displacement distance in (m)
        | Note: Above values are for a logarithmic velocity profile as:
        |     u=(u*/K)*ln((z-d)/z0)

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        z=np.arange(0.1,1,0.05)
        u=2/0.4*np.log((z-0.003)/0.002)
        ustar,z0,d=sm.surfaceroughness(z,u,np.max(z),'yes')

    References
    ----------

    Karimpour, A., Kaye, N. B., & Baratian-Ghorghi, Z. (2012). 
    Modeling the neutrally stable atmospheric boundary layer for laboratory scale studies of the built environment. 
    Building and Environment, 49, 203-211.

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
    
    z=type2numpy(z)
    u=type2numpy(u)

    #--------------------------------------------------------------------------
    #Assign default values

    if delta is None: delta=np.max(z)

    #--------------------------------------------------------------------------
    #Calculating shear velocity (friction velocity)

    umaxminusu=np.max(u)-u #umax-u
    z_delta=z/delta
    
    z_delta_min=0.1 #Lower limit that data used in curve fitting
    z_delta_max=0.5 #Upper limit that data used in curve fitting
    
    umaxminusu=umaxminusu[(z_delta>=z_delta_min) & (z_delta<=z_delta_max)]
    z_delta=z_delta[(z_delta>=z_delta_min) & (z_delta<=z_delta_max)]
    
    x=z_delta.copy()
    y=umaxminusu.copy()
    p=np.polyfit(np.log(x),y,1)
    
    ustar=-0.4*p[0] #Shear velocity (friction velocity), u*

    #--------------------------------------------------------------------------
    #Calculating surface roughness and zero plane displacement distance

    exp_u_k_ustar=np.exp(u*(0.4/ustar))
    
    x=z.copy()
    y=exp_u_k_ustar.copy()
    p=np.polyfit(x,y,1)
    
    z0=1/p[0] #Surface roughness, z0
    d=np.abs(p[1]*z0) #Zero plane displacement distance

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
    
        ufit=ustar/0.4*np.log((z-d)/z0)
    
        plt.scatter(u,z,label='Data')
        plt.plot(ufit,z,label='Fitted Curve')
        plt.xlabel('Velocity')
        plt.ylabel('Elevation')
        plt.legend()


    #--------------------------------------------------------------------------
    #Outputs
    return ustar, z0, d

    #--------------------------------------------------------------------------
