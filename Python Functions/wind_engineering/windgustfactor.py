def windgustfactor(t, t0=3600, CalcMethod='cem', Iu=0, dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2017-12-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.windgustfactor
    ==========================

    .. code:: python

        G = scientimate.windgustfactor(t, t0=3600, CalcMethod='cem', Iu=0, dispout='no')

    Description
    -----------

    Convert wind velocity of duration t0 to t

    Inputs
    ------

    t
        | Wind averaging duration to convert to (e.g. Sustained Wind Duration) in (second) 
        | Example: for 10-min wind averaging duration: t=600 (s)
    t0=3600
        | Wind averaging duration that data originally averaged at in (second) 
        | Example: for 60-min wind averaging duration: t0=3600 (s)
        | If t0 is not defined, then t0 is considered as t0=3600 s
        | t0 can only have one element
    CalcMethod='cem'
        | Wind gust factor calculation method 
        | 'durst': Durst (1960)
        | 'cem': Coastal Engineering Manual, U.S. Army Corps of Engineers (2015)
        | 'cook': Cook (1985) which is 3600 s gust factor based on Wieringa (1973) 600 s gust factor
        | 'krayer': Krayer & Marshall (1992)
        | 'cem' for t<=3600 s is similar to Dusrt (1960) relationship 
        | 'cem' for t>3600 s is Cook (1985) relationship with Iu=0.155  
        | 'cook' for Iu=0.175 is the commonly plotted curve  
    Iu=0
        Wind longitudinal turbulence intensity 
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    G
        | Wind gust factor
        | G=U(t)/U(t0)
        | G=(wind velocity avereged over t seconds)/(wind velocity avereged over t0 seconds)
        | If t0=3600 s then G is calculated respect with 1-hour avereged wind velocity
        | If t0=3600 s then G=(wind velocity avereged over t seconds)/(wind velocity avereged over 3600 seconds)
        | If t0=3600 s then G=U(t)/U(3600)

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        t=600
        G=sm.windgustfactor(t)

        t=range(1,3601,1)
        t0=3600
        G=sm.windgustfactor(t,t0,'durst',0,'yes')

        t=range(1,36001,1)
        t0=3600
        G=sm.windgustfactor(t,t0,'cem',0,'yes')

        t=np.arange(1,3601,1)
        t0=3600
        G=sm.windgustfactor(t,t0,'cook',0.175,'yes')

        t=range(1,3601,1)
        t0=3600
        G=sm.windgustfactor(t,t0,'krayer',0,'yes')

    References
    ----------

    American Society of Civil Engineers, ASCE-7. (2005). 
    Minimum design loads for buildings and other structures (Vol. 7). 
    Amer Society of Civil Engineers.

    Cook, N. J. (1985). 
    The designer's guide to wind loading on building structures. Part I: Background, damage survey, wind data, and structural classification. 
    Building Research Establishment, Watford.
    G_t_3600=1+0.42*Iu*log(3600/t)

    Durst, C. S. (1960). 
    Wind speeds over short periods of time. 
    Meteor. Mag, 89(1056), 181-187.
    G_t_3600=0.977+(0.64/(1+(t/38.14)**0.685)) #from ASCE-7 (2005)
    G_t_3600=0.96+(0.648/(1+(t/38)**0.638)) #from Krayer and Marshall (1992)

    Krayer, W. R., & Marshall, R. D. (1992). 
    Gust factors applied to hurricane winds. 
    Bulletin of the American Meteorological Society, 73(5), 613-618.
    G_t_3600=0.96+(0.839/(1+(t/36.27)**0.655))

    U.S. Army Corps of Engineers (2015). 
    Coastal Engineering Manual. 
    Engineer Manual 1110-2-1100, Washington, D.C.: U.S. Army Corps of Engineers.

    Wieringa, J. (1973). 
    Gust factors over open water and built-up country. 
    Boundary-Layer Meteorology, 3(4), 424-441.
    G_t_600=1+(1.42+0.3013*log((600/t)-4))*Iu

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
    
    t=type2numpy(t)

    #--------------------------------------------------------------------------
    #Calculating wind gust factor

    #Calculating wind gust factor from Durst (1960)
    if CalcMethod=='durst':
    
        #Calculating wind gust factor respect with 1-hour (3600 s) avereged wind velocity
        if t0==3600:
            
            G_t_3600=np.ones(len(t)) #Pre-assigning
            G_t_3600=0.977+(0.64/(1+(t/38.14)**0.685)) #Ratio of wind velocity averaged over duration t to 3600 s based on Durst (1960)
            G=G_t_3600.copy() #Wind gust factor as G_t_t0=U(t)/U(3600)
    
        #Converting wind velocity of duration t0 to t, e.g. t0=2-min averaged to t=10-min averaged
        elif t0!=3600:
    
            G_t_t0=np.ones(len(t)) #Pre-assigning
            G_t_t0=(0.977+(0.64/(1+(t/38.14)**0.685)))/(0.977+(0.64/(1+(t0/38.14)**0.685))) #Ratio of wind velocity averaged over duration t to t0 based on Durst (1960)
            G=G_t_t0.copy() #Wind gust factor as G_t_t0=U(t)/U(t0)
    
    #Calculating wind gust factor from Coastal Engineering Manual, U.S. Army Corps of Engineers (2015)
    elif CalcMethod=='cem':
    
        #Checking input times
    
        if t0<1: t0=1 #Shortest Sustained Wind Duration is 1 seceond
        t[t<1]=1 #Shortest Sustained Wind Duration is 1 seceond
        
        if t0>36000: t0=36000 #Longest Sustained Wind Duration is 36000 seceond (10 hr)
        t[t>36000]=36000 #Longest Sustained Wind Duration is 36000 seceond (10 hr)
    
        #Calculating wind gust factor respect with 1-hour (3600 s) avereged wind velocity
        if t0==3600:
            
            G_t_3600=np.ones(len(t)) #Pre-assigning
        
            G_t_3600[t<1]=(1.277+0.296*np.tanh(0.9*np.log10(45/1))) #Ratio of wind velocity averaged over duration t to 3600 s based on CEM (2015)
            G_t_3600[((t>=1) & (t<=3600))]=(1.277+0.296*np.tanh(0.9*np.log10(45/t[((t>=1) & (t<=3600))]))) #Ratio of wind velocity averaged over duration t to 3600 s based on CEM (2015)
            G_t_3600[((t>3600) & (t<36000))]=(1.5334-0.15*np.log10(t[((t>3600) & (t<36000))])) #Ratio of wind velocity averaged over duration t to 3600 s based on CEM (2015)
            G_t_3600[t>=36000]=(1.5334-0.15*np.log10(36000)) #Ratio of wind velocity averaged over duration t to 3600 s based on CEM (2015)
        
#            for i in range(0,len(t),1):
#                if t[i]<1:
#                    G_t_3600[i]=(1.277+0.296*np.tanh(0.9*np.log10(45/1))) #Ratio of wind velocity averaged over duration t to 3600 s based on CEM (2015)
#                elif ((t[i]>=1) and (t[i]<=3600)):
#                    G_t_3600[i]=(1.277+0.296*np.tanh(0.9*np.log10(45/t[i]))) #Ratio of wind velocity averaged over duration t to 3600 s based on CEM (2015)
#                elif ((t[i]>3600) and (t[i]<36000)):
#                    G_t_3600[i]=(1.5334-0.15*np.log10(t[i])) #Ratio of wind velocity averaged over duration t to 3600 s based on CEM (2015)
#                elif t[i]>=36000:
#                    G_t_3600[i]=(1.5334-0.15*np.log10(36000)) #Ratio of wind velocity averaged over duration t to 3600 s based on CEM (2015)
        
            G=G_t_3600.copy() #Wind gust factor as G_t_t0=U(t)/U(3600)
    
        #Converting wind velocity of duration t0 to t, e.g. t0=2-min averaged to t=10-min averaged
        elif t0!=3600:
        
            G_t_t0=np.ones(len(t)) #Pre-assigning
        
            for i in range(0,len(t),1):
                if ((t0<=3600) and (t[i]<=3600)):
                    G1=(1.277+0.296*np.tanh(0.9*np.log10(45/t[i])))/(1.277+0.296*np.tanh(0.9*np.log10(45/t0))) #Ratio of wind velocity averaged over duration t to t0 based on CEM (2015)
                    G_t_t0[i]=G1 #Wind gust factor as G_t_t0=U(t)/U(t0) 
                elif ((t0<=3600) and (t[i]>3600)):
                    G1=(1.277+0.296*np.tanh(0.9*np.log10(45/3600)))/(1.277+0.296*np.tanh(0.9*np.log10(45/t0))) #Ratio of wind velocity averaged over duration 3600 to t0 based on CEM (2015)
                    G2=(1.5334-0.15*np.log10(t[i]))/(1.5334-0.15*np.log10(3600)) #Ratio of wind velocity averaged over duration t to 3600 based on CEM (2015)
                    G_t_t0[i]=G1*G2 #Wind gust factor as G_t_t0=U(t)/U(t0) 
                elif ((t0>3600) and (t[i]<=3600)):
                    G1=(1.277+0.296*np.tanh(0.9*np.log10(45/t[i])))/(1.277+0.296*np.tanh(0.9*np.log10(45/3600))) #Ratio of wind velocity averaged over duration t to 3600 based on CEM (2015)
                    G2=(1.5334-0.15*np.log10(3600))/(1.5334-0.15*np.log10(t0)) #Ratio of wind velocity averaged over duration 3600 to t0 based on CEM (2015)
                    G_t_t0[i]=G1*G2 #Wind gust factor as G_t_t0=U(t)/U(t0) 
                elif ((t0>3600) and (t[i]>3600)):
                    G2=(1.5334-0.15*np.log10(t[i]))/(1.5334-0.15*np.log10(t0)) #Ratio of wind velocity averaged over duration t to t0 based on CEM (2015)
                    G_t_t0[i]=G2 #Wind gust factor as G_t_t0=U(t)/U(t0) 
        
            G=G_t_t0.copy() #Wind gust factor as G_t_t0=U(t)/U(t0)
    
    
    #Calculating wind gust factor from Cook (1985)
    elif CalcMethod=='cook':
    
        #Calculating wind gust factor respect with 1-hour (3600 s) avereged wind velocity
        if t0==3600:
            
            G_t_3600=np.ones(len(t)) #Pre-assigning
            G_t_3600=1+0.42*Iu*np.log(3600/t) #Ratio of wind velocity averaged over duration t to 3600 s based on Cook (1985)
            G=G_t_3600.copy() #Wind gust factor as G_t_t0=U(t)/U(3600)
    
        #Converting wind velocity of duration t0 to t, e.g. t0=2-min averaged to t=10-min averaged
        elif t0!=3600:
    
            G_t_t0=np.ones(len(t)) #Pre-assigning
            G_t_t0=(1+0.42*Iu*np.log(3600/t))/(1+0.42*Iu*np.log(3600/t0)) #Ratio of wind velocity averaged over duration t to t0 based on Cook (1985)
            G=G_t_t0.copy() #Wind gust factor as G_t_t0=U(t)/U(t0)
    
    #Calculating wind gust factor from Krayer & Marshall (1992)
    elif CalcMethod=='krayer':
    
        #Calculating wind gust factor respect with 1-hour (3600 s) avereged wind velocity
        if t0==3600:
            
            G_t_3600=np.ones(len(t)) #Pre-assigning
            G_t_3600=0.96+(0.839/(1+(t/36.27)**0.655)) #Ratio of wind velocity averaged over duration t to 3600 s based on Krayer & Marshall (1992)
            G=G_t_3600.copy() #Wind gust factor as G_t_t0=U(t)/U(3600)
    
        #Converting wind velocity of duration t0 to t, e.g. t0=2-min averaged to t=10-min averaged
        elif t0!=3600:
    
            G_t_t0=np.ones(len(t)) #Pre-assigning
            G_t_t0=(0.96+(0.839/(1+(t/36.27)**0.655)))/(0.96+(0.839/(1+(t0/36.27)**0.655))) #Ratio of wind velocity averaged over duration t to t0 based on Krayer & Marshall (1992)
            G=G_t_t0.copy() #Wind gust factor as G_t_t0=U(t)/U(t0)
    
    #--------------------------------------------------------------------------
    #Converting wind velocity averaged over t0 to wind velocity averaged over t

    #Windvel_t=Windvel_t0*G #Wind velocity averaged over duration of t

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
    
        #Plotting
        #plt.plot(t,G)
        plt.semilogx(t,G)
    
        plt.xlabel('t')
        plt.ylabel('Wind Gust Factor, G=U(t)/U(t0)')
    
    #--------------------------------------------------------------------------
    #Outputs
    return G

    #--------------------------------------------------------------------------
