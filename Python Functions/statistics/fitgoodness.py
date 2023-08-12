def fitgoodness(x, y, dispout='no'):
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

    scientimate.fitgoodness
    =======================

    .. code:: python

        r, R2, RMSE, MAE, SI, NSE, d, Bias, NMBias, RE = scientimate.fitgoodness(x, y, dispout='no')

    Description
    -----------

    Calculate goodness of fit parameters

    Inputs
    ------

    x
        Dataset with true (exact or expected) values, such as theoretical values
    y
        | Dataset that needs to be evaluated, such as model results or estimated values
        | Accuracy of y dataset is evaluated against x dataset
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    r
        Pearson correlation coefficient
    R2
        Coefficient of determination
    RMSE
        Root mean square error
    MAE
        Mean absolute error
    SI
        Scatter index
    NSE
        Nash Sutcliffe efficiency coefficient
    d
        Index of agreement
    Bias
        Bias
    NMBias
        Normalized mean bias
    RE
        Relative error

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        x=(-0.1+(0.1-(-0.1)))*np.random.randn(1024*2)
        y=x+(-0.01+(0.01-(-0.01)))*np.random.randn(1024*2)
        r,R2,RMSE,MAE,SI,NSE,d,Bias,NMBias,RE=sm.fitgoodness(x,y,'yes')

        x=[1,2,3,4,5,6,7,8,9,10]
        y=[1.1,1.98,3.3,4.2,4.8,5.95,7.5,7.7,8.99,10.5]
        r,R2,RMSE,MAE,SI,NSE,d,Bias,NMBias,RE=sm.fitgoodness(x,y,'yes')

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
    y=type2numpy(y)

    #--------------------------------------------------------------------------
    #Removing NaN and Inf values

    #Total number of elements in dataset
    NTotal=len(x)
    
    #Finding NaN and Inf in x dataset
    Indx1=np.int_((np.nonzero((np.isnan(x)!=1) & (np.isinf(x)!=1)))[0])
    x1=x[Indx1]
    y1=y[Indx1]
    
    #Finding NaN and Inf in y dataset
    Indx2=np.int_((np.nonzero((np.isnan(y1)!=1) & (np.isinf(y1)!=1)))[0])
    x2=x1[Indx2]
    y2=y1[Indx2]
    
    #Datasets without NaN and Inf
    x=x2.copy()
    y=y2.copy()
    
    #Total number of elements in dataset after removal of NaN and Inf
    N=len(x)
    
    #Total number of NaN and Inf elements in dataset
    NNaNInf=NTotal-N

    #--------------------------------------------------------------------------
    #Calculating statistical properties for x and y data

    xmean=np.mean(x) #Mean value of x data
    xstd=np.std(x,ddof=1) #Standard deviation of x data
    
    ymean=np.mean(y) #Mean value of y data
    ystd=np.std(y,ddof=1) #Standard deviation of y data

    #--------------------------------------------------------------------------
    #Calculating goodness of fit parameters (quantitavie evaluation of model) 

    #Pearson's correlation coefficient
    r=(np.sum((x-xmean)*(y-ymean)))/(np.sqrt((np.sum((x-xmean)**2))*(np.sum((y-ymean)**2)))) #Pearson's correlation coefficient 
    
    #Coefficient of determination
    #R2=1-(np.sum((x-y)**2))/(np.sum((x-xmean)**2)) #Coefficient of determination
    R2=r**2 #Coefficient of determination
    
    #Root-mean-square error
    RMSE=np.sqrt((np.sum((y-x)**2))/N) #Root-mean-square error
    
    #Mean absolute error 
    MAE=(np.sum(np.abs(y-x)))/N #Mean absolute error
    
    #Scatter index
    SI=RMSE/((np.sum(x))/N) #Scatter index
    
    #Nash Sutcliffe efficiency coefficient
    NSE=1-((np.sum((x-y)**2))/(np.sum((x-xmean)**2))) #Nash Sutcliffe efficiency coefficient
    
    #Index of agreement 
    d=1-(np.sum((x-y)**2))/(np.sum((np.abs(y-xmean)+np.abs(x-xmean))**2)) #Index of agreement
    
    #Bias
    #Bias=((sum(y(:,1)))/N-(sum(x(:,1)))/N); #Bias
    Bias=(ymean-xmean) #Bias
    
    #Normalized mean bias
    #NMBias=((np.sum(y))/N-(np.sum(x))/N)/((np.sum(x))/N) #Normalized mean bias
    NMBias=(ymean-xmean)/(xmean) #Normalized mean bias
    
    #Relative error
    RE=(y-x)/x #Relative error

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':

        val=[R2,r,RMSE,MAE,SI,NSE,d,Bias,NMBias,xmean,xstd,ymean,ystd,NTotal,N,NNaNInf]
        name=['R2','r','RMSE','MAE','SI','NSE','d','Bias','NMBias','x mean','x st. dev.','y mean','y st. dev.','N Input','N Used','N NaN Inf']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            print('{0:10}= {1:0.10f}'.format(name[i],val[i])) 
    
        #Plotting
        plt.scatter(x,y,label='Data')
        plt.plot([np.min(x),np.max(x)],[np.min(x),np.max(x)],label='1:1 Line')
    
        plt.xlabel('x Data')
        plt.ylabel('y Data')
        plt.legend()

    #--------------------------------------------------------------------------
    #Outputs
    return r, R2, RMSE, MAE, SI, NSE, d, Bias, NMBias, RE

    #--------------------------------------------------------------------------
