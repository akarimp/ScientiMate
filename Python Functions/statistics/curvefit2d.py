def curvefit2d(x, y, MathExpression, coefIniGuess, fitmethod='fmin', dispout='no'):
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

    scientimate.curvefit2d
    ======================

    .. code:: python

        c, yPredicted = scientimate.curvefit2d(x, y, MathExpression, coefIniGuess, fitmethod='fmin', dispout='no')

    Description
    -----------

    Fit curve to 2 dimensinal input dataset

    Inputs
    ------

    x
        x data
    y
        y data
    MathExpression
        | Right hand side of the relation y=f(x) as 'f(x)'
        | Example: y=c[0]*x**2+c[1]*x+c[2] then MathExpression='c[0]*x**2+c[1]*x+c[2]'
        | Example: y=c[0]*exp(x)+c[1] then MathExpression='c[0]*np.exp(x)+c[1]'
        | Example: y=c[0]*sin(x)+c[1] then MathExpression='c[0]*np.sin(x)+c[1]'
        | Desired coefficients to be found should be as : c[0], c[1],...,c[n]
    coefIniGuess
        | Initial guess for desired coefficients to be found
        | coefIniGuess=[guess0,guess1,...,guessn]
        | guess0 is initial guess for c[0],...,guessn is initial guess for c[n] 
    fitmethod
        | Fitting method: 
        | 'lsq': curve-fitting using nonlinear least-squares  
        | 'fmin': curve-fitting by minimizing a sum of squared errors
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    c
        Desired coefficients to be found
    yPredicted
        Predicted value from fitted curve

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        x=np.linspace(0,10,100)
        y=0.5*x**2+2*x+10+(-10+(10-(-10)))*np.random.rand(100)
        c,yPredicted=sm.curvefit2d(x,y,'c[0]*x**2+c[1]*x+c[2]',[1,1,1],'fmin','yes')

        x=np.linspace(0,10,100)
        y=5*x**2+7+(-200+(200-(-200)))*np.random.rand(100)
        c,yPredicted=sm.curvefit2d(x,y,'c[0]*x**c[1]+c[2]',[1,2,10],'lsq','yes')

        x=np.linspace(0,10,100)
        y=0.5*np.exp(x)+100+(-200+(200-(-200)))*np.random.rand(100)
        c,yPredicted=sm.curvefit2d(x,y,'c[0]*np.exp(x)+c[1]',[1,1],'fmin','yes')

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
    coefIniGuess=type2numpy(coefIniGuess)

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
    #Fitting curve 

    if fitmethod=='lsq':
        
        #Constructing mathematical expression for defining function
        MathExpressionforFun='lambda x,*c :'+'('+MathExpression+')'
        
        #Define Lambda Anonymous Functions
        #fun=lambda x,*c : (c[0]*x**2+c[1]*x+c[2])
        fun = eval(MathExpressionforFun)
    
        #Calculating coefficient for input function
        c,_ = sp.optimize.curve_fit(fun,x,y,coefIniGuess) 
    
    elif fitmethod=='fmin': 
    
        #Constructing mathematical expression for defining function
        MathExpressionforFun='lambda c,x,y :'+'(sum(('+MathExpression+'-y)**2))' #Sum of squared error

        #Define Anonymous Functions
        # fun = lambda x :(np.sum((x[0]*Windvelall**x[1]+x[2])-Turbdespikeall)**2)
        fun = eval(MathExpressionforFun)
    
        #Calculating coefficient for input function
        #c = sp.optimize.fmin(fun,coefIniGuess,args=(x,y))
        optimize_result = sp.optimize.minimize(fun,coefIniGuess,args=(x,y))
        c=optimize_result.x #Solution of optimization

    #--------------------------------------------------------------------------
    #Calculating Predicted value from fitted curve

    MathExpressionforyPredicted='lambda x, *c :'+'('+MathExpression+')'
    fun = eval(MathExpressionforyPredicted)
    yPredicted=np.zeros(len(x))
    for i in range(0,len(x)):
        yPredicted[i]=fun(x[i],*c)

    #--------------------------------------------------------------------------
    #Calculating statistical properties for y and yPredicted data

    ymean=np.mean(y) #Mean value of y data
    ystd=np.std(y,ddof=1) #Standard deviation of y data

    yPredictedmean=np.mean(yPredicted) #Mean value of yPredicted data
    yPredictedstd=np.std(yPredicted,ddof=1) #Standard deviation of yPredicted data

    #--------------------------------------------------------------------------
    #Calculating goodness of fit parameters (quantitavie evaluation of model) 

    #Pearson's correlation coefficient
    r=(np.sum((y-ymean)*(yPredicted-yPredictedmean)))/(np.sqrt((np.sum((y-ymean)**2))*(np.sum((yPredicted-yPredictedmean)**2)))) #Pearson's correlation coefficient 

    #Coefficient of determination
    #R2=1-(np.sum((y-yPredicted)**2))/(np.sum((y-ymean)**2)) #Coefficient of determination
    R2=r**2 #Coefficient of determination

    #Root-mean-square error
    RMSE=np.sqrt((np.sum((yPredicted-y)**2))/N) #Root-mean-square error

    #Mean absolute error 
    MAE=(np.sum(np.abs(yPredicted-y)))/N #Mean absolute error

    #Scatter index
    SI=RMSE/((np.sum(y))/N) #Scatter index

    #Nash Sutcliffe efficiency coefficient
    NSE=1-((np.sum((y-yPredicted)**2))/(np.sum((y-ymean)**2))) #Nash Sutcliffe efficiency coefficient

    #Index of agreement 
    d=1-(np.sum((y-yPredicted)**2))/(np.sum((np.abs(yPredicted-ymean)+np.abs(y-ymean))**2)) #Index of agreement

    #Bias
    #Bias=((np.sum(yPredicted))/N-(np.sum(y))/N) #Bias
    Bias=(yPredictedmean-ymean) #Bias

    #Normalized mean bias
    #NMBias=((np.sum(yPredicted))/N-(np.sum(y))/N)/((np.sum(y))/N) #Normalized mean bias
    NMBias=(yPredictedmean-ymean)/(ymean) #Normalized mean bias

    #Relative error
    RE=(yPredicted-y)/y #Relative error

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
    
        for i in range(0,len(c)):
            val=[c[i]]
            name=['c'+str(i)]
            print('{0:10}= {1:0.10f}'.format(name[0],val[0])) 

        val=[RMSE,R2,NTotal,N,NNaNInf]
        name=['RMSE','R2','N Input','N Used','N NaN Inf']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            print('{0:10}= {1:0.10f}'.format(name[i],val[i])) 

        #Plotting
        plt.scatter(x,y,label='Data')
        #plt.hold='True'
        plt.plot(x,yPredicted,label='Fitted Curve')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.legend()

    #--------------------------------------------------------------------------
    #Outputs
    return c, yPredicted

    #--------------------------------------------------------------------------
