def curvefit3d(x, y, z, MathExpression, coefIniGuess, fitmethod='fmin', dispout='no'):
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

    scientimate.curvefit3d
    ======================

    .. code:: python

        c, zPredicted = scientimate.curvefit3d(x, y, z, MathExpression, coefIniGuess, fitmethod='fmin', dispout='no')

    Description
    -----------

    Fit curve to 3 dimensinal input dataset

    Inputs
    ------

    x
        x data
    y
        y data
    z
        z data
    MathExpression
        | Right hand side of the relation z=f(x,y) as 'f(x,y)'
        | Example: z=c[0]*x**2+c[1]*y+c[2] then MathExpression='c[0]*x**2+c[1]*y+c[2]'
        | Example: z=c[0]*exp(x)+c[1]*sin(y) then MathExpression='c[0]*np.exp(x)+c[1]*np.sin(y)'
        | Desired coefficients to be found should be as : c[0], c[1],...,c[n]
    coefIniGuess
        | Initial guess for desired coefficients to be found
        | coefIniGuess=[guess0,guess1,...,guessn]
        | guess0 is initial guess for c[0],...,guessn is initial guess for c(n) 
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
    zPredicted
        Predicted value from fitted curve

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        xx,yy=np.meshgrid(np.linspace(0,10,100),np.linspace(0,10,100))
        zz=xx**2+yy**2+10+(-10+(10-(-10)))*np.random.rand(100,100)
        x=xx.reshape(-1) #Flatten array
        y=yy.reshape(-1) #Flatten array
        z=zz.reshape(-1) #Flatten array
        c,zPredicted=sm.curvefit3d(x,y,z,'c[0]*x**2+c[1]*y**2',[1,1],'fmin','yes')

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
    
    x=type2numpy(x)
    y=type2numpy(y)
    z=type2numpy(z)
    coefIniGuess=type2numpy(coefIniGuess)

    #--------------------------------------------------------------------------
    #Removing NaN and Inf values

    #Total number of elements in dataset
    NTotal=len(x)
    
    #Finding NaN and Inf in x dataset
    Indx1=np.int_((np.nonzero((np.isnan(x)!=1) & (np.isinf(x)!=1)))[0])
    #Indx1=~(np.isnan(x) | np.isinf(x))
    x1=x[Indx1]
    y1=y[Indx1]
    z1=z[Indx1]

    #Finding NaN and Inf in y dataset
    Indx2=np.int_((np.nonzero((np.isnan(y1)!=1) & (np.isinf(y1)!=1)))[0])
    #Indx2=~(np.isnan(y1) | np.isinf(1))
    x2=x1[Indx2]
    y2=y1[Indx2]
    z2=z1[Indx2]

    #Finding NaN and Inf in z dataset
    Indx3=np.int_((np.nonzero((np.isnan(z2)!=1) & (np.isinf(z2)!=1)))[0])
    #Indx3=~(np.isnan(z2) | np.isinf(z2))
    x3=x2[Indx3]
    y3=y2[Indx3]
    z3=z2[Indx3]

    #Datasets without NaN and Inf
    if fitmethod=='fmin':
        x=x3.copy()
        y=y3.copy()
        z=z3.copy()

    #Total number of elements in dataset after removal of NaN and Inf
    N=len(x)
    
    #Total number of NaN and Inf elements in dataset
    NNaNInf=NTotal-N

    #--------------------------------------------------------------------------
    #Fitting curve 

    if fitmethod=='lsq':
        
        #Constructing mathematical expression for defining function
        MathExpressionforFun='lambda x,*c :'+'('+MathExpression+')'
        
        #Define Anonymous Functions
        #fun=lambda x,*c : (c[0]*x**2+c[1]*y+c[2])
        fun = eval(MathExpressionforFun)
    
        #Calculating coefficient for input function
        c,_ = sp.optimize.curve_fit(fun,x,z,coefIniGuess) 
    
    elif fitmethod=='fmin': 
    
        #Constructing mathematical expression for defining function
        MathExpressionforFun='lambda c,x,y,z :'+'(sum(('+MathExpression+'-z)**2))' #Sum of squared error
    
        #Define Anonymous Functions
        # fun = lambda x :(np.sum((x[0]*Windvelall**x[1]+x[2])-Turbdespikeall)**2)
        fun = eval(MathExpressionforFun)
    
        #Calculating coefficient for input function
        #c = sp.optimize.fmin(fun,coefIniGuess,args=(x,y,z))
        optimize_result = sp.optimize.minimize(fun,coefIniGuess,args=(x,y,z))
        c=optimize_result.x #Solution of optimization

    #--------------------------------------------------------------------------
    #Calculating Predicted value from fitted curve

    MathExpressionforyPredicted='lambda x,y,*c :'+'('+MathExpression+')'
    fun = eval(MathExpressionforyPredicted)
    zPredicted=np.zeros(len(x))
    for i in range(0,len(x)):
        zPredicted[i]=fun(x[i],y[i],*c)

    #--------------------------------------------------------------------------
    #Calculating statistical properties for z and zPredicted data

    zmean=np.mean(z) #Mean value of z data
    zstd=np.std(z,ddof=1) #Standard deviation of z data

    zPredictedmean=np.mean(zPredicted) #Mean value of zPredicted data
    zPredictedstd=np.std(zPredicted,ddof=1) #Standard deviation of zPredicted data

    #--------------------------------------------------------------------------
    #Calculating goodness of fit parameters (quantitavie evaluation of model) 

    #Pearson's correlation coefficient
    r=(np.sum((z-zmean)*(zPredicted-zPredictedmean)))/(np.sqrt((np.sum((z-zmean)**2))*(np.sum((zPredicted-zPredictedmean)**2)))) #Pearson's correlation coefficient 

    #Coefficient of determination
    #R2=1-(np.sum((z-zPredicted)**2))/(np.sum((z-zmean)**2)) #Coefficient of determination
    R2=r**2 #Coefficient of determination

    #Root-mean-square error
    RMSE=np.sqrt((np.sum((zPredicted-z)**2))/N) #Root-mean-square error

    #Mean absolute error 
    MAE=(np.sum(np.abs(zPredicted-z)))/N #Mean absolute error

    #Scatter index
    SI=RMSE/((np.sum(z))/N) #Scatter index

    #Nash Sutcliffe efficiency coefficient
    NSE=1-((np.sum((z-zPredicted)**2))/(np.sum((z-zmean)**2))) #Nash Sutcliffe efficiency coefficient

    #Index of agreement 
    d=1-(np.sum((z-zPredicted)**2))/(np.sum((np.abs(zPredicted-zmean)+np.abs(z-zmean))**2)) #Index of agreement

    #Bias
    #Bias=((np.sum(zPredicted))/N-(np.sum(z))/N) #Bias
    Bias=(zPredictedmean-zmean) #Bias

    #Normalized mean bias
    #NMBias=((np.sum(zPredicted))/N-(np.sum(z))/N)/((np.sum(z))/N) #Normalized mean bias
    NMBias=(zPredictedmean-zmean)/(zmean) #Normalized mean bias

    #Relative error
    RE=(zPredicted-z)/z #Relative error

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
        xx,yy=np.meshgrid(np.linspace(np.min(x),np.max(x),25),np.linspace(np.min(y),np.max(y),25))
        zPredictedforplot=np.zeros((len(xx[:,0]),len(yy[:,0])))
        for i in range(0,len(xx[:,0])):
            for j in range(0,len(yy[:,0])):
                zPredictedforplot[i,j]=fun(xx[i,j],yy[i,j],*c)
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1, projection='3d')
        ax.scatter(x,y,z,label='Data')
        #plt.hold='True'
        ax.plot_surface(xx,yy,zPredictedforplot,cmap=plt.cm.get_cmap(),label='Fitted Curve')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        #ax.legend()

    #--------------------------------------------------------------------------
    #Outputs
    return c, zPredicted

    #--------------------------------------------------------------------------
