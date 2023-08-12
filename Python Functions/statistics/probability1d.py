def probability1d(x, binedge=None, dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2017-06-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.probability1d
    =========================

    .. code:: python

        fdensity, fdensitycumulative, bincenter, xmean, xstd = scientimate.probability1d(x, binedge=None, dispout='no')

    Description
    -----------

    Calculate 1D probability density distribution for a given dataset

    Inputs
    ------

    x
        Input data 
    binedge
        | Bin edges  
        | length(binedge)=number of bin +1   
        | If there are N bins in histogram/distribution, then values in binedge are as:   
        | 1st value: left edge of 1st bin, 2nd value: left edge of 2nd bin, ...   
        | (N)th value: left edge of last bin, (N+1)th value: right edge of last bin   
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    fdensity
        Probability density distribution
    fdensitycumulative
        Cumulative probability density distribution
    bincenter
        Bin center
    xmean
        Mean value of input data
    xstd
        Standard deviation of input data

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        x=(-0.1+(0.1-(-0.1)))*np.random.randn(1024*2)
        binedge=np.linspace(min(x),max(x),11)
        fdensity,fdensitycumulative,bincenter,xmean,xstd=sm.probability1d(x,binedge,'yes')

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
    #Assign default values

    if binedge is not None: binedge=type2numpy(binedge)
    if binedge is None: binedge=np.linspace(np.min(x),np.max(x),11)

    #--------------------------------------------------------------------------
    #Removing NaN and Inf values

    #Total number of elements in dataset
    NTotal=len(x)
    
    #Finding NaN and Inf in x dataset
    Indx1=np.int_((np.nonzero((np.isnan(x)!=1) & (np.isinf(x)!=1)))[0])
    x1=x[Indx1]
    
    #Datasets without NaN and Inf
    x=x1.copy()

    #Total number of elements in dataset after removal of NaN and Inf
    N=len(x)
    
    #Total number of NaN and Inf elements in dataset
    NNaNInf=NTotal-N

    #--------------------------------------------------------------------------
    #Calculating probability density distribution

    binwidth=np.diff(binedge) #Bin width, length(binwidth) is one element less than length(binedge) 
    bincenter=binedge[0:-1]+binwidth/2 #Bin center
    
    NPoints,binedgeout=np.histogram(x,binedge) #Calculating number of data points in each bin
    histarea=NPoints*binwidth #Area under each element of histogram
    fdensity=NPoints/(np.sum(histarea)) #Calculating probability density distribution
    #Note: np.sum(fdensity*binwidth)=1
    
    #--------------------------------------------------------------------------
    #Calculating cumulative probability density distribution

    fdensitycumulative=np.zeros(len(bincenter)) #Pre-assigning vector
    for i in range(0,len(bincenter),1):
        fdensitycumulative[i]=np.sum(fdensity[0:i+1]*binwidth[0:i+1])

    #--------------------------------------------------------------------------
    #Calculating statistical properties

    xmean=np.mean(x) #Mean value of the input data
    xstd=np.std(x,ddof=1) #Standard deviation of the input data

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
    
        val=[xmean,xstd,NTotal,N,NNaNInf]
        name=['mean','st. dev.','N Input','N Used','N NaN Inf']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            #print('{0:10}= {1:0.10f}'.format(name[i],val[i])) 
            print('     {0:15}     {1:+0.10g}'.format(name[i],val[i])) 

        #Plotting
        plt.subplot(2,1,1)
        plt.bar(bincenter,fdensity,width=0.8*binwidth)
        plt.title('Probability Density Distribution')
        plt.ylabel('Probability Density, f')
        plt.xlabel('Data')
        
        plt.subplot(2,1,2)
        plt.plot(bincenter,fdensitycumulative)
        plt.title('Cumulative Probability Density Distribution')
        plt.ylabel('Cumulative Probability Density, f')
        plt.xlabel('Data')
    
    #--------------------------------------------------------------------------
    #Outputs
    return fdensity, fdensitycumulative, bincenter, xmean, xstd

    #--------------------------------------------------------------------------
