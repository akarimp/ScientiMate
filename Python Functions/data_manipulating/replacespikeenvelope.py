def replacespikeenvelope(x, lowbound, upbound, nrpeat=1, interpMethod='linear', dispout='no'):
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

    scientimate.replacespikeenvelope
    ================================

    .. code:: python

        xDespiked, Indx = scientimate.replacespikeenvelope(x, lowbound, upbound, nrpeat=1, interpMethod='linear', dispout='no')

    Description
    -----------

    Remove spikes in the time series that are outside a defined envelope

    Inputs
    ------

    x
        Input data
    lowbound
        Lower boundary of the data, all data points should be larger than that
    upbound
        Upper boundary of the data, all data points should be smaller than that
    nrpeat=1
        Number of time despiking procedure is repeating
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
        x[220:225]=1.5
        x=x+5
        xDespiked,Indx=sm.replacespikeenvelope(x,3.9,6.1,1,'linear','yes')

        fs=2
        t=np.linspace(0,1023.5,1024*fs)
        x=sp.signal.detrend(0.5*np.cos(2*np.pi*0.2*t)+(-0.1+(0.1-(-0.1)))*np.random.rand(1024*fs))
        spikeloc=np.arange(10,len(t),100)
        x[np.int64(spikeloc+np.round(2*np.random.randn(len(spikeloc))))]=np.sign(np.random.randn(len(spikeloc)))
        xDespiked,Indx=sm.replacespikeenvelope(x,-0.6,0.6,1,'linear','yes')

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

    #Preserving an input data
    xInput=x.copy()

    #--------------------------------------------------------------------------
    #Removing spikes

    for i in range(0,nrpeat,1):
        
        #Checking the defined boundary
        Indx1=np.int_((np.nonzero((x>lowbound) & (x<upbound)))[0])
        if len(Indx1)==0:
            print('Warning-----------------------------------------')
            print('No data points are between defined boundaries')
            print('Data are not de-spiked')
            print('------------------------------------------------')
            lowbound=np.min(x)
            upbound=np.max(x)
        
        #Defining spike points, if a point locatated outside boundary limits is considered as spike
        Indx2=np.int_((np.nonzero(x<lowbound))[0])
        Indx3=np.int_((np.nonzero(x>upbound))[0])
    
        #Setting spike points as NaN
        xDespiked=x.copy()
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
