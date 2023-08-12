def directionalpsd(Syy, f, Wavedir=0, calcMethod='mitsuyasu', Windvel=0, dtheta=15, dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2017-05-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.directionalpsd
    ==========================

    .. code:: python

        Syy2d, f2d, theta = scientimate.directionalpsd(Syy, f, Wavedir=0, calcMethod='mitsuyasu', Windvel=0, dtheta=15, dispout='no')

    Description
    -----------

    | Calculate wave directional spectrum using parametric directional spreading function
    | such as Mitsuyasu et al. (1975), Hasselmann et al. (1980) and Donelan et al. (1985)

    Inputs
    ------

    Syy
        One dimensional power spectral density (m^2/Hz)
    f
        Frequency (Hz)
    Wavedir=0
        Mean wave direction between 0 and 360 (Degree)
    CalcMethod='mitsuyasu'
        | Directional wave spectrum calculation method 
        | 'pierson': Pierson et al. (1955), 'cos2': D=1/pi*(cos((theta-theta_mean)/2))^2 
        | 'mitsuyasu': Mitsuyasu (1975), 'hasselmann': Hasselmann (1980), 'donelan': Donelan et al. (1985) 
        | 'flat': uniform distribution in all directions
    Windvel=0
        | Wind velocity at 10 meter above surface level in (m/s)
        | Wind velocity is required for Mitsuyasu (1975) and Hasselmann (1980) methods
        | For other methods use Windvel=0
    dtheta=15
        Direction interval at which directional spectrum calculated between 0 and 360 (Degree)
    dispout='no'
        | Define to display outputs or not
        | '2d': 2 dimensional plot, 'surface': Surface plot, 'polar': Polar plot, 'no': not display 

    Outputs
    -------

    Syy2d
        Directional wave power spectral density (m^2/Hz/Degree)
    f2d
        Directional frequency (Hz)
    theta
        Direction (Degree)

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        N=2**11 #Total number of points
        fs=8 #Sampling frequency
        df=fs/N #Frequency difference 
        f=np.arange(0,fs/2+df,df) #Frequency vector 
        f[0]=f[1]/2 #Assign temporarily non-zero value to fisrt element of f to prevent division by zero
        Syy=0.016*9.81**2/((2*np.pi)**4*(f**5))*np.exp(-1.25*(0.33/f)**4) #Calculating Spectrum 
        f[0]=0
        Syy[0]=0
        Syy2d,f2d,theta=sm.directionalpsd(Syy,f,45,'mitsuyasu',10,15,'polar')

    References
    ----------

    Banner, M. L. (1990). 
    Equilibrium spectra of wind waves. 
    Journal of Physical Oceanography, 20(7), 966-984.

    Donelan, M. A., Hamilton, J., & Hui, W. (1985). 
    Directional spectra of wind-generated waves. 
    Philosophical Transactions of the Royal Society of London A: Mathematical, Physical and Engineering Sciences, 315(1534), 509-562.

    Ewans, K. C. (1998). 
    Observations of the directional spectrum of fetch-limited waves. 
    Journal of Physical Oceanography, 28(3), 495-512.

    Goda, Y. (1999). 
    A comparative review on the functional forms of directional wave spectrum. 
    Coastal Engineering Journal, 41(01), 1-20.

    Hasselmann, D. E., Dunckel, M., & Ewing, J. A. (1980). 
    Directional wave spectra observed during JONSWAP 1973. 
    Journal of physical oceanography, 10(8), 1264-1280.

    Hwang, P. A., & Wang, D. W. (2001). 
    Directional distributions and mean square slopes in the equilibrium and saturation ranges of the wave spectrum. 
    Journal of physical oceanography, 31(5), 1346-1360.

    Mitsuyasu, H., Tasai, F., Suhara, T., Mizuno, S., Ohkusu, M., Honda, T., & Rikiishi, K. (1975). 
    Observations of the directional spectrum of ocean WavesUsing a cloverleaf buoy. 
    Journal of Physical Oceanography, 5(4), 750-760.

    Pierson, W. J., Neumann, G., & James, R. W. (1955). 
    Practical methods for observing and forecasting ocean waves by means of wave spectra and statistics. 
    Publication 603, U.S. Navy Hydrographic Office, 284 pp. 

    Sorensen, R. M. (2005). 
    Basic coastal engineering (Vol. 10). 
    Springer Science & Business Media.

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
    from scipy import special
    if dispout!='no':
        import matplotlib.pyplot as plt 
        from mpl_toolkits.mplot3d import Axes3D

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
    
    Syy=type2numpy(Syy)
    f=type2numpy(f)

    #--------------------------------------------------------------------------

    df=f[1]-f[0] #Frequency interval
    Wavedirrad=np.deg2rad(Wavedir) #Converting to radian

    #--------------------------------------------------------------------------
    #Calculating wave properties

    #Calculating spectral moments
    m0=np.sum(Syy*f**0*df)
    m1=np.sum(Syy*f**1*df)
    m2=np.sum(Syy*f**2*df)

    #Calculating wave properties
    Hm0=4*np.sqrt(m0) #Zero-Moment wave height
    Tm01=m0/m1 #mean period
    Tm02=(m0/m2)**0.5 #zero crossing period

    #Calculating peak period
    SyymaxIndx=np.argmax(Syy) #Locating an index of the spectrum peak
    fp=f[SyymaxIndx] #peak frequency
    Tp=1/fp #peak period

    #Calculating peak frequency from weighted integral (Young, 1995)
    #fp=(np.sum(Syy**5*f**1*df))/(np.sum(Syy**5*f**0*df)) #Peak frequency

    #--------------------------------------------------------------------------
    #Generating direction vector

    #thetadeg=np.linspace(0,360,(360+dtheta)/dtheta)
    thetadeg=np.arange(0,360+0.001,dtheta) #Direction vector from 0 to 360 Degree, equally spaced at dtheta
    thetarad=np.deg2rad(thetadeg) #Converting to radian

    #--------------------------------------------------------------------------
    #Directioal Spectrum Calculation (Ewans, 1998; Goda, 1999; Hwang and Wang, 2001; Sorensen, 2005)

    #Wave phase velocity (wave celerity) in a deep water at peak wave frequency
    Cp=9.81/(2*np.pi*fp)
    
    #Calculating directional spreading function
    if calcMethod=='pierson': #Pierson et al. (1955)
        
        #Directional spreading function: D(f,theta)
        #St. Denis and Pierson (1953): D(f,theta)=2/pi*(cos(theta-theta_mean))^2
    
        Gs=2/np.pi
    
        #Calculating directional spreading function
        f2d=np.zeros((len(f),len(thetadeg))) #Pre-assigning vector
        D=np.zeros((len(f),len(thetadeg))) #Pre-assigning vector
        for j in range(1,len(thetadeg),1):
            f2d[:,j]=f.copy() #Two dimensional frequency vector
            D[:,j]=(Gs*np.abs(np.cos(thetarad[j]-Wavedirrad))**(2))
            if ((Wavedirrad>0) and (Wavedirrad<=np.pi/2)): 
                if (np.abs((thetarad[j]-1*np.pi)-(Wavedirrad))<np.pi/2): 
                    D[:,j]=0 #D=0 for abs(theta)>pi/2
            elif ((Wavedirrad>np.pi/2) and (Wavedirrad<=3*np.pi/2)): 
                if (np.abs(thetarad[j]-Wavedirrad)>np.pi/2): 
                    D[:,j]=0 #D=0 for abs(theta)>pi/2
            elif ((Wavedirrad>3*np.pi/2) and (Wavedirrad<=2*np.pi)): 
                if (np.abs((thetarad[j]+1*np.pi)-(Wavedirrad))<np.pi/2): 
                    D[:,j]=0 #D=0 for abs(theta)>pi/2
    
    elif calcMethod=='cos2': #Mitsuyasu et al. (1975) for s=1
        
        #Directional spreading function: D(f,theta)
        #Mitsuyasu et al. (1975): D(f,theta)=G(s)*(cos((theta-theta_mean)/2))^2s
        #If s=1, then D(f,theta)=G(s)*(cos((theta-theta_mean)/2))^2
        #or, D(f,theta)=1/pi*(cos((theta-theta_mean)/2))^2
    
        #Calculating spreading parameter, s
        s=1
    
        Gs=2**(2*s-1)/(np.pi)*(sp.special.gamma(s+1))**2/(sp.special.gamma(2*s+1))
    
        #Calculating directional spreading function
        f2d=np.zeros((len(f),len(thetadeg))) #Pre-assigning vector
        D=np.zeros((len(f),len(thetadeg))) #Pre-assigning vector
        for j in range(1,len(thetadeg),1):
            f2d[:,j]=f.copy() #Two dimensional frequency vector
            D[:,j]=(Gs*np.abs(np.cos((thetarad[j]-Wavedirrad)/2))**(2*s))

    elif calcMethod=='mitsuyasu': #Mitsuyasu et al. (1975)
    
        Sp=11.5*(Windvel/Cp)**(-2.5)
    
        #Calculating spreading parameter, s
        s=np.zeros(len(f)) #Pre-assigning s vector
        s[f<fp]=Sp*(f[f<fp]/fp)**5
        s[f>=fp]=Sp*(f[f>=fp]/fp)**(-2.5)
    
        Gs=2**(2*s-1)/(np.pi)*(sp.special.gamma(s+1))**2/(sp.special.gamma(2*s+1))
    
        #Calculating directional spreading function
        #f2d=np.zeros((len(f),len(thetadeg))) #Pre-assigning vector
        #D=np.zeros((len(f),len(thetadeg))) #Pre-assigning vector
        #for i in range(0,len(f),1):
        #    for j in range(1,len(thetadeg),1):
        #        f2d[i,j]=f[i] #Two dimensional frequency vector
        #        D[i,j]=np.real(Gs[i]*np.abs(np.cos((thetarad[j]-Wavedirrad)/2))**(2*s[i]))
    
        #Calculating directional spreading function
        f2d=np.zeros((len(f),len(thetadeg))) #Pre-assigning vector
        D=np.zeros((len(f),len(thetadeg))) #Pre-assigning vector
        for j in range(1,len(thetadeg),1):
            f2d[:,j]=f.copy() #Two dimensional frequency vector
            D[:,j]=(Gs*np.abs(np.cos((thetarad[j]-Wavedirrad)/2))**(2*s))
    
    elif calcMethod=='hasselmann': #Hasselmann et al. (1980)
    
        mu=-2.33-1.45*(Windvel/Cp-1.17)
    
        #Calculating spreading parameter, s
        s=np.zeros(len(f)) #Pre-assigning s vector
        s[f<1.05*fp]=6.97*(f[f<1.05*fp]/fp)**4.06
        s[f>=1.05*fp]=9.77*(f[f>=1.05*fp]/fp)**mu
    
        Gs=2**(2*s-1)/(np.pi)*(sp.special.gamma(s+1))**2/(sp.special.gamma(2*s+1))
    
        #Calculating directional spreading function
        f2d=np.zeros((len(f),len(thetadeg))) #Pre-assigning vector
        D=np.zeros((len(f),len(thetadeg))) #Pre-assigning vector
        for j in range(1,len(thetadeg),1):
            f2d[:,j]=f.copy() #Two dimensional frequency vector
            D[:,j]=(Gs*np.abs(np.cos((thetarad[j]-Wavedirrad)/2))**(2*s))
    
    elif calcMethod=='donelan': #Donelan et al. (1985)
    
        #Calculating spreading parameter, Beta,  Donelan et al. (1985) and Banner (1990)
        Beta=np.zeros(len(f)) #Pre-assigning vector
        Beta[f<=0.56*fp]=1.24
        Beta[((f>0.56*fp) & (f<=0.95*fp))]=2.61*(f[((f>0.56*fp) & (f<=0.95*fp))]/fp)**1.3
        Beta[((f>0.95*fp) & (f<=1.6*fp))]=2.28*(f[((f>0.95*fp) & (f<=1.6*fp))]/fp)**-1.3
        Beta[f>1.6*fp]=10**(-0.4+0.8393*np.exp(-0.567*np.log((f[f>1.6*fp]/fp)**2)))
    
        #Calculating directional spreading function
        f2d=np.zeros((len(f),len(thetadeg))) #Pre-assigning vector
        D=np.zeros((len(f),len(thetadeg))) #Pre-assigning vector
        for j in range(1,len(thetadeg),1):
            f2d[:,j]=f.copy() #Two dimensional frequency vector
            D[:,j]=(0.5*Beta*(np.cosh(Beta*(thetarad[j]-Wavedirrad)))**(-2))
            if (((thetarad[j]-Wavedirrad)>np.pi) and ((thetarad[j]-Wavedirrad)<2*np.pi)): 
                D[:,j]=(0.5*Beta*(np.cosh(Beta*((thetarad[j]-2*np.pi)-Wavedirrad)))**(-2))
    
    elif calcMethod=='flat': #Uniform distribution in all directions 
        
        Gs=1/(2*np.pi)
    
        #Calculating directional spreading function
        f2d=np.zeros((len(f),len(thetadeg))) #Pre-assigning vector
        D=np.zeros((len(f),len(thetadeg))) #Pre-assigning vector
        for j in range(1,len(thetadeg),1):
            f2d[:,j]=f.copy() #Two dimensional frequency vector
            D[:,j]=Gs
    
    #Normalizing directional spreading function
    for i in range(0,len(f),1):
        D[i,:]=D[i,:]/(np.sum(D[i,:])*dtheta)
    
    Syy2d=np.zeros((len(f),len(thetadeg))) #Pre-assigning
    #Calculating directional power specral density
    for i in range(0,len(f),1):
        for j in range(0,len(thetadeg),1):
            Syy2d[i,j]=Syy[i]*D[i,j] #Directional spectrum

    #--------------------------------------------------------------------------
    #Calculating 1-D spectrum by integrating from directional spectrum

    Syy1d=np.zeros(len(f)) #Pre-assigning vector
    D1d=np.zeros(len(f)) #Pre-assigning vector
    for i in range(0,len(f),1):
        Syy1d[i]=np.sum(np.real(Syy2d[i,:])*dtheta)
        D1d[i]=np.sum(np.real(D[i,:])*dtheta)

    #--------------------------------------------------------------------------
    #Assigning output parameter

    theta=thetadeg.copy() #Direction vector

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='2d': #2 dimensional plot
    
        val=[Hm0,fp,Tp,Tm01,Tm02,m0,m1,m2]
        name=['Hm0','fp','Tp','Tm01','Tm02','m0','m1','m2']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            print('{0:10}= {1:0.10f}'.format(name[i],val[i])) 
    
        #Plotting
        plt.plot(f,Syy)
    
        plt.title('Power Spectral Density')
        plt.xlabel('Frequency(Hz)')
        plt.ylabel('Spectral Density(m^2/Hz)')
    
    elif dispout=='surface': #Surface plot
    
        val=[Hm0,fp,Tp,Tm01,Tm02,m0,m1,m2]
        name=['Hm0','fp','Tp','Tm01','Tm02','m0','m1','m2']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            print('{0:10}= {1:0.10f}'.format(name[i],val[i])) 
    
        thetaGrid,fGrid=np.meshgrid(thetadeg,f)

        fig=plt.figure()
        fig3d=fig.gca(projection='3d')
        #fig3d.view_init(elev=50,azim=-20)
        surf1=fig3d.plot_surface(fGrid,thetaGrid,Syy2d,cmap=plt.get_cmap(),edgecolor='none')
        fig.colorbar(surf1)
        fig3d.set_title('Directional Power Spectral Density')
        fig3d.set_xlabel('Frequency (Hz)')
        fig3d.set_ylabel('Direction (Degree)')
        fig3d.set_zlabel('Power Spectral Density (m^2/Hz/Degree)')
    
    elif dispout=='polar': #Polar plot
    
        val=[Hm0,fp,Tp,Tm01,Tm02,m0,m1,m2]
        name=['Hm0','fp','Tp','Tm01','Tm02','m0','m1','m2']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            print('{0:10}= {1:0.10f}'.format(name[i],val[i])) 
    
        thetaGrid,fGrid=np.meshgrid(thetadeg,f)
        ax=plt.subplot(111, projection='polar')
        count1=ax.contour(np.deg2rad(thetaGrid),fGrid,Syy2d,cmap=plt.get_cmap())

        plt.colorbar(count1)
        plt.title('Directional Power Spectral Density')
        #plt.xlabel('Direction (Degree)')
        #plt.ylabel('Frequency (Hz)')
        #plt.zlabel('Spectral Density (m^2/Hz/degree)')
    
    #--------------------------------------------------------------------------
    #Outputs
    return Syy2d, f2d, theta

    #--------------------------------------------------------------------------
