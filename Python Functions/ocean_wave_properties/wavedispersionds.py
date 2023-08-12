def wavedispersionds(h, T, Uc=0):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2017-02-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.wavedispersionds
    ============================

    .. code:: python

        k, L, C, Cg = scientimate.wavedispersionds(h, T, Uc=0)

    Description
    -----------

    | Solve water wave dispersion relation with presence of current (Doppler shift)
    | Calculate wave number (k), wave length (L), wave celereity (C), and wave group velocity (Cg) using linear wave theory

    Inputs
    ------

    h
        Water depth in (m)
    T
        | Wave period in (s) 
        | If peak wave frequency (Tp) is used, calculated values represent peak wave 
    Uc=0
        | Current velocity in (m/s), for Doppler shift
        | Uc should have a same size as h
        | Note: inputs can be as a single value or a 1-D vertical array

    Outputs
    -------

    k
        Wave number in (radian/m)
    L
        Wave length in (m)
    C
        Wave celerity in (m/s)
    Cg
        Wave group celerity in (m/s)

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        k,L,C,Cg=sm.wavedispersionds(1,3,1)

        k,L,C,Cg=sm.wavedispersionds([1,1.1],[3,3.1],[1,1])

        k,L,C,Cg=sm.wavedispersionds(np.array([1,1.1]),np.array([3,3.1]),np.array([1,1]))

    References
    ----------


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
    
    h=type2numpy(h)
    T=type2numpy(T)
    Uc=type2numpy(Uc)

    #--------------------------------------------------------------------------
    #Calculating wave number (k)

    #Make Uc the same size as h
    if len(Uc)==1:
        Uc1=Uc
        Uc=np.zeros(len(h))
        Uc[:]=Uc1

    f=1/T #Wave frequency
    w=2*np.pi*f #Wave angular frequency, w=2*pi/T=2*pi*f
    
    #Deep water
    k0=(w**2)/9.81 #Deep water wave number
    k0h=k0*h
    
    
    #Calculating exact wave number (k)
    
    #Estimation of wave number (k) from Beji (2013)
    kh=k0h*(1+k0h**1.09*np.exp(-(1.55+1.3*k0h+0.216*k0h**2)))/np.sqrt(np.tanh(k0h)) #Calculating wave number from Beji (2013)
    kini=kh/h #Initial value for k (Wave number from Beji (2013))
    kini[w==0]=0

    #Calculating exact wave number (k)
    k=np.zeros(len(f)) #Pre-assigning array to make program run faster
    for i in range(0,len(f),1):
        if len(h)==1:
            fun = lambda x : (w[i]-x*Uc)**2-(9.81*x*np.tanh(x*h)) #(w-kU)**2=g.k.tanh(kd)
        else:
            fun = lambda x : (w[i]-x*Uc[i])**2-(9.81*x*np.tanh(x*h[i])) #(w-kU)**2=g.k.tanh(kd)

        k[i]=sp.optimize.fsolve(fun,kini[i]) #Wave number

    #--------------------------------------------------------------------------
    #Calculating wave length (L), wave celereity (C), and wave group velocity (Cg)

    L=2*np.pi/k #Wave length
    C=9.81*T/(2*np.pi)*np.tanh(k*h) #Wave celerity
    Cg=1/2*C*(1+(2*k*h)/np.sinh(2*k*h)) #Wave group celirity

    #--------------------------------------------------------------------------
    #Outputs
    return k, L, C, Cg

    #--------------------------------------------------------------------------
