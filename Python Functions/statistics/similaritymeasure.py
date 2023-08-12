def similaritymeasure(x, y, CalcMethod='euclidean', dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2020-01-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.similaritymeasure
    =============================

    .. code:: python

        d = scientimate.similaritymeasure(x, y, CalcMethod='euclidean', dispout='no')

    Description
    -----------

    Measure similarity between two arrays

    Inputs
    ------

    x
        First array, its similarity is measured against y array
    y
        Second array, its similarity is measured against x array
    CalcMethod='euclidean'
        | Similarity  calculation method
        | 'euclidean': Euclidean distance
        | 'manhattan': Manhattan distance
        | 'minkowski': Minkowski distance (power=3)
        | 'cosine': Cosine distance
        | 'pearson': Pearson's correlation coefficient
        | 'spearman': spearman's correlation coefficient
        | 'norm': Absolute difference of norm
        | 'covariance': Covariance
        | 'inv_covariance': Euclidean distance of inverse covariance
        | 'histogram': Mean of absolute difference of histogram
        | 't-test': Two-sample t-test statistic
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    d
        Arrays similarity measure

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        x=[0,2,4,6]
        y=[2,3,5,7]
        d=sm.similaritymeasure(x,y,'euclidean','yes')

        x=[1,2,3,4,5,6,7,8,9,10]
        y=[1.1,1.98,3.3,4.2,4.8,5.95,7.5,7.7,8.99,10.5]
        d=sm.similaritymeasure(x,y,'pearson','yes')

    References
    ----------
    Kianimajd, A., Ruano, M. G., Carvalho, P., Henriques, J., Rocha, T., Paredes, S., & Ruano, A. E. (2017).
    Comparison of different methods of measuring similarity in physiologic time series.
    IFAC-PapersOnLine, 50(1), 11005-11010.

    * https://en.wikipedia.org/wiki/Similarity_measure
    * https://en.wikipedia.org/wiki/Goodness_of_fit
    * https://dataaspirant.com/2015/04/11/five-most-popular-similarity-measures-implementation-in-python/
    * https://towardsdatascience.com/similarity-measures-e3dbd4e58660
    * https://www.mathworks.com/matlabcentral/answers/377944-how-to-calculate-a-percentage-of-similarity-between-two-arrays
    * https://en.wikipedia.org/wiki/Template_matching
    * https://www.mathworks.com/help/images/ref/normxcorr2.html
    * https://www.mathworks.com/help/stats/hypothesis-tests-1.html
    
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
    from scipy import stats
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

    #Finding NaN and Inf in x dataset
    x[((np.isnan(x)==1) | (np.isinf(x)==1))]=0
    
    #Finding NaN and Inf in y dataset
    y[((np.isnan(x)==1) | (np.isinf(x)==1))]=0
    
    #--------------------------------------------------------------------------
    #Calculating statistical properties for x and y data

    xmean=np.mean(x) #Mean value of x data
    xstd=np.std(x,ddof=1) #Standard deviation of x data
    
    ymean=np.mean(y) #Mean value of y data
    ystd=np.std(y,ddof=1) #Standard deviation of y data

    #--------------------------------------------------------------------------
    #Calculate similarity

    #Calculate similarity using Euclidean distance
    if CalcMethod=='euclidean':
        #Or d=np.sqrt(np.sum((x-y)**2))
        p=2
        d=(np.sum((np.abs(x-y))**p))**(1/p)
    
    #Calculate similarity using Manhattan distance
    elif CalcMethod=='manhattan':
        #Or d=np.sum(np.abs(x-y))
        p=1
        d=(np.sum((np.abs(x-y))**p))**(1/p)
    
    #Calculate similarity using Minkowski distance
    elif CalcMethod=='minkowski':
        p=3
        d=(np.sum((np.abs(x-y))**p))**(1/p)
    
    #Calculate similarity using Cosine distance
    elif CalcMethod=='cosine':
        d=(np.sum(x*y))/(np.sqrt(np.sum(x**2))*np.sqrt(np.sum(y**2)))

    #Pearson's correlation coefficient
    elif CalcMethod=='pearson':
        #Or d=(np.sum((x-xmean)*(y-ymean)))/(np.sqrt((np.sum((x-xmean)**2))*(np.sum((y-ymean)**2)))) #Pearson's correlation coefficient
        d, p_val = sp.stats.pearsonr(np.reshape(x,(-1)),np.reshape(y,(-1))) #Pearson's correlation coefficient

    #Spearman's Spearman coefficient
    elif CalcMethod=='spearman':
        d, p_val = sp.stats.spearmanr(np.reshape(x,(-1)), np.reshape(y,(-1))) #Spearman's correlation coefficient

    #Absolute difference of norm
    #https://www.mathworks.com/matlabcentral/answers/19752-metrics-for-matrices-similarity
    elif CalcMethod=='norm':
        d=np.abs(np.linalg.norm(np.reshape(x,(-1)))-np.linalg.norm(np.reshape(y,(-1))))

    #Covariance
    elif CalcMethod=='covariance':
        cov_matrix=np.cov(np.reshape(x,(-1)),np.reshape(y,(-1)))
        d=cov_matrix[0,1]

    #Euclidean distance of inverse covariance
    elif CalcMethod=='inv_covariance':
        cov_x=np.cov(x)
        cov_y=np.cov(y)
        if ((np.ndim(x)>1) and (np.ndim(x)>1)):
            cov_x_inv=np.linalg.inv(cov_x)
            cov_y_inv=np.linalg.inv(cov_y)
            d=np.sqrt(np.sum((cov_x_inv-cov_y_inv)**2)) #Euclidean distance
        elif ((np.size(cov_x)==1) and (np.size(cov_y)==1)):
            d=np.sqrt(np.sum((1/cov_x-1/cov_y)**2)) #Euclidean distance
        else:
            print('Cannot calculate inverse covariance, inputs should have the same dimensions')

    #Mean of absolute difference of histogram
    elif CalcMethod=='histogram':
        d=np.mean((np.abs(np.histogram(np.reshape(x,(-1)))[0]-np.histogram(np.reshape(y,(-1)))[0])))

    #Two-sample t-test statistic
    elif CalcMethod=='t-test':
        d, p_val = sp.stats.ttest_ind(np.reshape(x,(-1)), np.reshape(y,(-1)), equal_var = False)


    #Normalize similarity value
    #d=d/np.max(d)  #Normalize to [0, 1]

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':

        val=[d]
        name=['Similarity measure']
        for i in range(0,len(val)):
            #print ('\n',name[i],val[i],sep='     ',end='')
            print('{0:10}= {1:0.10f}'.format(name[i],val[i])) 
    
        plt.subplot(1,2,1)
        x_reshp=np.reshape(x,(-1,1))
        y_reshp=np.reshape(y,(-1,1))
        x_y_data=np.concatenate((x_reshp,y_reshp), axis=1)
        plt.hist(x_y_data,label=['x','y'])
        plt.xlabel('Bins')
        plt.ylabel('Frequency')
        plt.legend()


        plt.subplot(1,2,2)
        plt.scatter(x,y,label='Data')
        plt.plot([np.min(x),np.max(x)],[np.min(x),np.max(x)],label='1:1 Line')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.legend()


    #--------------------------------------------------------------------------
    #Outputs
    return d

    #--------------------------------------------------------------------------
