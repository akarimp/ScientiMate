def replacespike3dps(x,nrpeat=1,DespikeScaleFactor=1,interpMethod='linear',dispout='no'):
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

    scientimate.replacespike3dps
    ============================

    .. code:: python

        xDespiked, Indx = scientimate.replacespike3dps(x, nrpeat=1, DespikeScaleFactor=1, interpMethod='linear', dispout='yes')

    Description
    -----------

    Remove spikes in the time series based on 3D phase space method by Goring and Nikora (2002)

    Inputs
    ------

    x
        Input data
    nrpeat=1
        Number of time despiking procedure is repeating
    DespikeScaleFactor=1
        Scaling a despiking threshold
    interpMethod='linear'
        | Interpolation method for replacing spike points:
        | Matlab/Octave: 'linear', 'nearest', 'next', 'previous', 'pchip', 'cubic', 'spline'
        | Python/Scipy : 'linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic'
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    xDespiked
        Dispiked data
    Indx
        Index of despiked points

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np
        import scipy as sp
        from scipy import signal

        fs=128
        t=np.linspace(0,9.5,10*fs)
        x=np.sin(2*np.pi*0.3*t)+0.1*np.sin(2*np.pi*4*t)
        spikeloc=np.arange(10,len(t),100)
        x[np.int64(spikeloc+np.round(2*np.random.randn(len(spikeloc))))]=np.sign(np.random.randn(len(spikeloc)))
        x[220:225]=1.5
        x=x+5
        xDespiked,Indx=sm.replacespike3dps(x,2,1,'linear','yes')

        fs=2
        t=np.linspace(0,1023.5,1024*fs)
        x=sp.signal.detrend(0.5*np.cos(2*np.pi*0.2*t)+(-0.1+(0.1-(-0.1)))*np.random.rand(1024*fs))
        spikeloc=np.arange(10,len(t),100)
        x[np.int64(spikeloc+np.round(2*np.random.randn(len(spikeloc))))]=np.sign(np.random.randn(len(spikeloc)))
        xDespiked,Indx=sm.replacespike3dps(x,1,1,'linear','yes')

    References
    ----------

    Goring, D. G., & Nikora, V. I. (2002). 
    Despiking acoustic Doppler velocimeter data. 
    Journal of Hydraulic Engineering, 128(1), 117-126.

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
    
    x=type2numpy(x)

    #--------------------------------------------------------------------------

    N=len(x) #Length of the input data

    #Preserving an input data
    xInput=x.copy()
    xInputmean=np.mean(xInput) #Mean of input data

    #Removing a mean from input data
    x=x-xInputmean

    #--------------------------------------------------------------------------
    #Removing spikes

    for i in range(0,nrpeat,1):
    
        #Gradient
        x1stGrad=np.gradient(x) #First gradient of x
        x2ndGrad=np.gradient(x1stGrad) #Second gradient of x
    
        #Standard deviation
        xStd=np.std(x) #Standard deviation of x
        x1stGradStd=np.std(x1stGrad) #Standard deviation of x1stGrad
        x2ndGradStd=np.std(x2ndGrad) #Standard deviation of x2ndGrad
    
        #Rotation  angle
        #Theta=np.arctan2(np.sum(x*x2ndGrad),sum(x**2)) #Rotation  angle  of  the  principal  axis of x2ndGrad vs x
        Theta=np.arctan(np.sum(x*x2ndGrad)/np.sum(x**2)) #Rotation  angle  of  the  principal  axis of x2ndGrad vs x
         
        #Expected absolute maximum based on normal distributaion
        Lambda=np.sqrt(2*np.log(N))
    
        #Axis of ellipse
        a1=Lambda*xStd #Major axis of ellipse for x1stGrad vs x
        b1=Lambda*x1stGradStd #Minor axis of ellipse for x1stGrad vs x
    
        a3=Lambda*x1stGradStd #Major axis of ellipse for x2ndGradStd vs x1stGrad
        b3=Lambda*x2ndGradStd #Minor axis of ellipse for x2ndGradStd vs x1stGrad
    
        #Initial values for a2 and b2
        a2=a1.copy() #Major axis of ellipse for x2ndGradStd vs x
        b2=b1.copy() #Minor axis of ellipse for x2ndGradStd vs x
        a2prevstep=0
        b2prevstep=0
    
        #Loop with 5 iteration to define a2 and b2 
        for n in range(0,5,1):
            a2=np.sqrt(((Lambda*xStd)**2-(b2*np.sin(Theta))**2)/((np.cos(Theta))**2)) #Major axis of ellipse for x2ndGradStd vs x
            b2=np.sqrt(((Lambda*x2ndGradStd)**2-(a2*np.sin(Theta))**2)/((np.cos(Theta))**2)) #Minor axis of ellipse for x2ndGradStd vs x
            a2diff=a2-a2prevstep
            b2diff=b2-b2prevstep
            a2prevstep=a2.copy()
            b2prevstep=b2.copy()

    
        #Additional loop (if needed) with (100-5) iterations to define a2 and b2
        while ((n<=100) and ((np.abs(a2diff)>0.0001) or (np.abs(b2diff)>0.0001))):
            a2=np.sqrt(((Lambda*xStd)**2-(b2*np.sin(Theta))**2)/((np.cos(Theta))**2)) #Major axis of ellipse for x2ndGradStd vs x
            b2=np.sqrt(((Lambda*x2ndGradStd)**2-(a2*np.sin(Theta))**2)/((np.cos(Theta))**2)) #Minor axis of ellipse for x2ndGradStd vs x
            a2diff=a2-a2prevstep
            b2diff=b2-b2prevstep
            a2prevstep=a2.copy()
            b2prevstep=b2.copy()
            n=n+1
    
        #Calculating mean values
        xmean=np.mean(x)
        x1stGradmean=np.mean(x1stGrad)
        x2ndGradmean=np.mean(x2ndGrad)
    
        #Calculating distances based on ellipse equation
        distellipse12=((x-xmean)**2/a1**2+(x1stGrad-x1stGradmean)**2/b1**2)
        distellipse13=((x-xmean)**2/a2**2+(x2ndGrad-x2ndGradmean)**2/b2**2)
        distellipse23=((x1stGrad-x1stGradmean)**2/a3**2+(x2ndGrad-x2ndGradmean)**2/b3**2)
    
        #Defining spike points, if a point locatated outside ellipse it is considered as spike
        Indx1=np.int_((np.nonzero(distellipse12>1*DespikeScaleFactor))[0])
        Indx2=np.int_((np.nonzero(distellipse13>1*DespikeScaleFactor))[0])
        Indx3=np.int_((np.nonzero(distellipse23>1*DespikeScaleFactor))[0])
    
        #Fixing an indexing, method may pick 2 points right before and after the spike
        for m in range(1,len(Indx1),1):
            if Indx1[m]-Indx1[m-1]==2:
                Indxmean=(Indx1[m]+Indx1[m-1])/2
                Indx1[m]=Indxmean
                Indx1[m-1]=Indxmean

    
        for m in range(1,len(Indx2),1):
            if Indx2[m]-Indx2[m-1]==2:
                Indxmean=(Indx2[m]+Indx2[m-1])/2
                Indx2[m]=Indxmean
                Indx2[m-1]=Indxmean

    
        for m in range(1,len(Indx3),1):
            if Indx3[m]-Indx3[m-1]==2:
                Indxmean=(Indx3[m]+Indx3[m-1])/2
                Indx3[m]=Indxmean
                Indx3[m-1]=Indxmean

    
        #Setting spike points as NaN
        xDespiked=x.copy()
        xDespiked[Indx1]=np.nan
        xDespiked[Indx2]=np.nan
        xDespiked[Indx3]=np.nan
    
        #Locating all spike pints
        Indx4=np.int_((np.nonzero(np.isnan(xDespiked)==1))[0])
    
        #Locating all non-spike pints
        Indx5=np.int_((np.nonzero(np.isnan(xDespiked)!=1))[0])
    
        #Replacing spike points
        samples=np.arange(1,len(x)+1,1)
        if len(Indx4)!=0:
            #xDespiked[Indx4]=np.interp1(samples[Indx4],samples[Indx5],xDespiked[Indx5],interpMethod,'extrap')
            fun=sp.interpolate.interp1d(samples[Indx5],xDespiked[Indx5],kind=interpMethod,fill_value='extrapolated')
            xDespiked[Indx4]=fun(samples[Indx4])
    
        #Assigning a despiked data to x for next repeat
        x=xDespiked.copy()
    
    
    #Adding mean back to data
    xDespiked=xDespiked+xInputmean

    #--------------------------------------------------------------------------
    #Defining despiked points

    Indx=np.int_((np.nonzero(xInput!=xDespiked))[0])
    NumberDespikedPoint=len(Indx)
    PercentDespikedPoint=len(Indx)/len(x)*100

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':

        print('Total number of pints despiked = ',str(NumberDespikedPoint))
        print('Percent of despiked points     = ',str(round(PercentDespikedPoint,2)),' %')
        
        plt.subplot(2,1,1)
        plt.plot(samples,xInput)
        #plt.hold='True'
        plt.scatter(samples[Indx],xInput[Indx])
        plt.xlabel('Sample')
        plt.ylabel('Data')
        plt.title('Input Data')
    
        plt.subplot(2,1,2)
        plt.plot(samples,xDespiked)
        #plt.hold='True'
        plt.scatter(samples[Indx],xDespiked[Indx])
        plt.xlabel('Sample')
        plt.ylabel('Data')
        plt.title('Despiked Data')
    

    #--------------------------------------------------------------------------
    #Outputs
    return xDespiked, Indx

    #--------------------------------------------------------------------------
