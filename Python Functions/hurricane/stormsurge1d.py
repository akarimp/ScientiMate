def stormsurge1d(h0, U10, m=0, dispout='no'):
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

    scientimate.stormsurge1d
    ========================

    .. code:: python

        Eta, x, maxsurgeheight, L = scientimate.stormsurge1d(h0, U10, m=0, dispout='no')

    Description
    -----------

    Calculate one dimensional storm surge using Dean Dalrymple (1991) method

    Inputs
    ------

    h0
        Deep-water depth in (m), 
    U10
        10-min averaged wind velocity at 10 meter above a surface in (m/s)
    m=0
        | Bed slope
        | Note: m=h0/L, where L is a length of the continental shelf in (m)
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    Eta
        Surge height along a x axis in (m)
    x
        | Points on x axis in (m) 
        | x=0 located at a far-end of the water boundary
        | x=L located at a coastline
    maxsurgeheight
        Maximum surge height at a coastline (x=L) in (m)
    L
        | Length of the continental shelf in (m)
        | L=h0/m if m is not zero
        | L=50000 if m is zero (L=50000 meter for the hurricane Katrina)

    Examples
    --------

    .. code:: python

        import scientimate as sm

        h0=42 #Example: h0=42 m for the hurricane Katrina
        U10=40 #Example: U10=40 m/s for the hurricane Katrina 
        m=0.00084
        Eta,x,maxsurgeheight,L=sm.stormsurge1d(h0,U10,m,'yes')

    References
    ----------

    Dean, R. G., & Dalrymple, R. A. (1991). 
    Water wave mechanics for engineers and scientists (Vol. 2). 
    World Scientific Publishing Company.

    Wu, J. (1982). 
    Wind‐stress coefficients over sea surface from breeze to hurricane. 
    Journal of Geophysical Research: Oceans, 87(C12), 9704-9706.

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
    
    h0=type2numpy(h0)
    U10=type2numpy(U10)
    m=type2numpy(m)

    #--------------------------------------------------------------------------
    #Initial values

    n=1.3 #Shear stress coefficient
    Row=1028 #Seawater density in (kg/m3)
    Rhoa=1.204 #Air density in (kg/m3)
    
    #Calculating a length of the continental shelf
    if m==0:
        L=50000 #L=50000 m for the hurricane Katrina
    else:
        L=h0/m

    #--------------------------------------------------------------------------
    #Calculating wind drag coefficient based on Wu (1982), method used by ADCIRC

    CD=(0.8+0.065*U10)*1e-3 #Wind drag coefficient
    CD[U10<7.5]=1.2875*1e-3 #Wind drag coefficient
    #CD[CD>0.0034]=0.0034 #Wind drag coefficient
    Tauwind=Rhoa*CD*U10**2 #Wind shear stress in (N/m^2)
    A=n*Tauwind*L/(Row*9.81*h0**2)

    #--------------------------------------------------------------------------
    #Calculating the surge

    x=np.linspace(0,L,1000) #Points on a x axis
    
    Eta=np.zeros(len(x)) #Pre-assign array
    #Bed without slope
    if m==0:
        
        #Claculating Eta for a bed with a zero slope
        Eta=h0*(np.sqrt(1+2*A*x/L)-1)
        
    #Bed with slpoe
    else:
        
        h=-m*x+h0 #h=h0.*(1-x./L);
        
        #Claculating Eta for a bed with a slope
        Eta[0]=0
        
        for i in range(1,len(x),1):
            
            Etaini=Eta[i-1] #Initial value for Eta from previous step
            fun = lambda Eta1 : (x[i]/L-((1-(h[i]+Eta1)/h0)-A*np.log(((h[i]+Eta1)/h0-A)/(1-A)))) ##x/L=(1-(h+Eta)/h0)-A*log(((h+Eta)/h0-A)/(1-A))
            Eta[i]=sp.optimize.fsolve(fun,Etaini) #Surge height

    
    maxsurgeheight=Eta[-1]

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
    
        #Colors: https://htmlcolorcodes.com/
        
        if m==0:
            
            #Plotting a seabed
            x1=[0,L,L,L+L/10]
            y1=[-h0,-h0,0,0]
            plt.plot(x1,y1,color=np.array([160,64,0])/255,label='Seabed')
            
            #Plotting a water surface
            x2=[0,L]
            y2=[0,0]
            plt.plot(x2,y2,color=np.array([41,128,185])/255,label='Water Surface')
            
            #Plotting a storm surge
            plt.plot(x,Eta,color=np.array([231,76,60])/255,label='Storm Surge')
    
            plt.xlim(0,L+L/10)
            plt.xlabel('x')
            plt.ylabel('Elevation')
            plt.legend(loc='lower right')
            
        else:
            
            #Plotting a seabed
            x1=[0,L,L+L/10]
            y1=[-h0,m*L-h0,m*L-h0]
            plt.plot(x1,y1,color=np.array([160,64,0])/255,label='Seabed')
            
            #Plotting a water surface
            x2=[0,L]
            y2=[0,0]
            plt.plot(x2,y2,color=np.array([41,128,185])/255,label='Water Surface')
    
            #Plotting a storm surge
            plt.plot(x,Eta,color=np.array([231,76,60])/255,label='Storm Surge')
    
            plt.xlim(0,L+L/10)
            plt.xlabel('x')
            plt.ylabel('Elevation')
            plt.legend(loc='lower right')
    

    #--------------------------------------------------------------------------
    #Output
    return Eta, x, maxsurgeheight, L

    #--------------------------------------------------------------------------
