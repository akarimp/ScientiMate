def replacemissing1d(x, what2replace='both', interpMethod='linear', dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2017-02-01/2020-02-01                  +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.replacemissing1d
    ============================

    .. code:: python

        xReplaced, NaN_Indx = scientimate.replacemissing1d(x, what2replace='both', interpMethod='linear', dispout='no')

    Description
    -----------

    Replace missing data points in 1d data such as time series

    Inputs
    ------

    x
        Input data
    what2replace='both'
        | What needs to be replaced
        | 'NaN': replacing NaN data points
        | 'Inf': replacing Inf data points
        | 'both': replacing NaN and Inf data points
        | Number: replacing data points equal to Number
    interpMethod='linear'
        | Interpolation method for replacing spike points:
        | Matlab/Octave: 'linear', 'nearest', 'next', 'previous', 'pchip', 'cubic', 'spline'
        | Python/Scipy : 'linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic'
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    xReplaced
        Replaced data
    NaN_Indx
        Logical index of replaced points

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
        x[np.int64(spikeloc+np.round(2*np.random.randn(len(spikeloc))))]=np.nan
        x=x+5
        x[220:225]=np.nan
        xReplaced,NaN_Indx=sm.replacemissing1d(x,'NaN','linear','yes')

        fs=2
        t=np.linspace(0,1023.5,1024*fs)
        x=sp.signal.detrend(0.5*np.cos(2*np.pi*0.2*t)+(-0.1+(0.1-(-0.1)))*np.random.rand(1024*fs))
        spikeloc=np.arange(10,len(t),100)
        x=x+1
        x[np.int64(spikeloc+np.round(2*np.random.randn(len(spikeloc))))]=0
        xReplaced,NaN_Indx=sm.replacemissing1d(x,0,'linear','yes')

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
    #Replacing missing data points

    #Assigning empty to NaN_Indx
    NaN_Indx=[]
    
    #Defining missing data points
    if what2replace=='NaN':
        #Define NaN data points
        NaN_Indx=(np.isnan(x)==1)
        Valid_Indx=(np.isnan(x)!=1)
    elif what2replace=='Inf':
        #Define Inf data points
        NaN_Indx=(np.isinf(x)==1)
        Valid_Indx=(np.isinf(x)!=1)
    elif what2replace=='both':
        #Define NaN and Inf data points
        NaN_Indx=((np.isnan(x)==1) | (np.isinf(x)==1))
        Valid_Indx=((np.isnan(x)!=1) | (np.isinf(x)!=1))
    elif ((isinstance(what2replace,int)==True) or (isinstance(what2replace,float)==True)):
        #Define NaN and Inf data points
        NaN_Indx=(x==what2replace)
        Valid_Indx=(x!=what2replace)

    #Assigning x to xReplaced
    xReplaced=x.copy()

    #Replacing missing points
    samples=np.arange(1,len(x)+1,1)
    if True in NaN_Indx:
        #xReplaced[NaN_Indx]=np.interp(samples[NaN_Indx==True],samples[Valid_Indx==True],xReplaced[Valid_Indx==True])
        fun=sp.interpolate.interp1d(samples[Valid_Indx==True],xReplaced[Valid_Indx==True],kind=interpMethod,bounds_error=False,fill_value='extrapolate')
        xReplaced[NaN_Indx==True]=fun(samples[NaN_Indx==True])
        
    #Assigning a despiked data to x for next repeat
    #x=xReplaced.copy()


    #--------------------------------------------------------------------------
    #Defining replaced points

    if True in NaN_Indx:
        NumberReplacedPoint=np.sum(NaN_Indx)
        PercentReplacedPoint=np.sum(NaN_Indx)/np.size(x)*100
    else:
        NumberReplacedPoint=0
        PercentReplacedPoint=0

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':

        print('Total number of pints replaced = ',str(NumberReplacedPoint))
        print('Percent of replaced points     = ',str(round(PercentReplacedPoint,2)),' %')
        
        plt.subplot(2,1,1)
        plt.plot(samples,xInput)
        plt.xlabel('Sample')
        plt.ylabel('Data')
        plt.title('Input Data')
    
        plt.subplot(2,1,2)
        plt.plot(samples,xReplaced)
        plt.scatter(samples[NaN_Indx],xReplaced[NaN_Indx])
        plt.xlabel('Sample')
        plt.ylabel('Data')
        plt.title('Replaced Data')
    

    #--------------------------------------------------------------------------
    #Outputs
    return xReplaced, NaN_Indx

    #--------------------------------------------------------------------------
