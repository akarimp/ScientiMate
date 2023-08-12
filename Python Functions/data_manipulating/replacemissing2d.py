def replacemissing2d(x, what2replace='both', gridsize_x=1, gridsize_y=1, interpMethod='nearest', dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2020-02-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    scientimate.replacemissing2d
    ============================

    .. code:: python
    
        xReplaced, NaN_Indx = scientimate.replacemissing2d(x, what2replace='both', gridsize_x=1, gridsize_y=1, interpMethod='nearest', dispout='no')
    
    Description
    -----------
    
    Replace missing data points in 2d array
    
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
    gridsize_x=1
        | Grid size (distance between grid points) in x direction
        | Leave gridsize_x=1 if you do not have it
    gridsize_y=1
        | Grid size (distance between grid points) in y direction
        | Leave gridsize_y=1 if you do not have it
    interpMethod='nearest'
        | Interpolation method
        | 'linear': Use default or 'linear' method to interpolate
        | 'nearest': Use nearest neighbor method to interpolate
        | 'knn': Use nearest neighbor method to interpolate (Use 'knn' for large array)
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
        from numpy import random
    
        x=[[1,0,3],[2,5,np.nan],[3,np.nan,1],[5,7,2]]
        xReplaced, NaN_Indx = sm.replacemissing2d(x, 'NaN', 1, 1, 'nearest', 'yes')
    
        rng = np.random.default_rng()
        xgrid=rng.normal(size=(100,50))
        xgrid(rng.integers(0,99,(20,1)),rng.integers(0,49,(10,1)))=np.nan
        xReplaced, NaN_Indx = sm.replacemissing2d(xgrid, 'NaN', 1, 1, 'knn', 'yes')
    
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
        import matplotlib as mpl
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
    
    #Create 2d i and j positions
    M,N=np.shape(x)
    samples_ii,samples_jj=np.meshgrid(range(0,N,1),range(0,M,1))
    samples_ii=gridsize_x*samples_ii #Map position to have correct grid length
    samples_jj=gridsize_y*samples_jj #Map position to have correct grid length
    
    #Replacing missing points
    if interpMethod=='knn':
    
        if True in NaN_Indx:

            NaN_row_col=np.argwhere(NaN_Indx==True) #Find row and columns that have NaN values
            NaN_row=NaN_row_col[:,0]
            NaN_col=NaN_row_col[:,1]
            
            #for i=1:M
                #for j=1:N
            for i in range(0,len(NaN_row),1):
                Indx_ii=NaN_row[i] #Index of NaN point in x direction
                Indx_jj=NaN_col[i] #Index of NaN point in y direction
                if NaN_Indx[Indx_ii,Indx_jj]==True:
                    Dist=np.sqrt((samples_ii-samples_ii[Indx_ii,Indx_jj])**2+(samples_jj-samples_jj[Indx_ii,Indx_jj])**2) #Calculating distance of (x,y) to (xPoint,yPoint)
                    Dist[NaN_Indx==True]=np.max(Dist) #make sure a NaN point is not replace by another NaN point (NaN point distance to itself is zero)
                    MinDist_Indx=np.argmin(np.reshape(Dist,-1)) #Finding minimum distances
                    x_1d=np.reshape(xReplaced,-1)
                    xReplaced[Indx_ii,Indx_jj]=x_1d[MinDist_Indx] #Value of the nearest point to (Indx_ii,Indx_jj)

    
    else:
    
        if True in NaN_Indx:
            xReplaced[NaN_Indx==True]=sp.interpolate.griddata((samples_ii[Valid_Indx==True],samples_jj[Valid_Indx==True]),xReplaced[Valid_Indx==True],(samples_ii[NaN_Indx==True],samples_jj[NaN_Indx==True]),method=interpMethod)
    
            #Replacing NaN data point resulted from default method with ones from nearest method
            if np.sum(np.isnan(xReplaced))>0:
                xReplacednearest=sp.interpolate.griddata((samples_ii,samples_jj),xReplaced,(samples_ii,samples_jj),method='nearest')
                xReplaced[np.isnan(xReplaced)==True]=xReplacednearest[np.isnan(xReplaced)==True] 


    #Assigning a replaced data to x for next repeat
    x=xReplaced.copy()
    
    
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
        plt.imshow(xInput,interpolation=None,cmap=mpl.cm.get_cmap(),aspect='auto',origin='lower')
        plt.title('Input Data')
    
        plt.subplot(2,1,2)
        plt.imshow(xReplaced,interpolation=None,cmap=mpl.cm.get_cmap(),aspect='auto',origin='lower')
        plt.scatter(samples_ii[NaN_Indx]/gridsize_x,samples_jj[NaN_Indx]/gridsize_y,c='r')
        plt.title('Replaced Data')
        
    
    #--------------------------------------------------------------------------
    #Outputs
    return xReplaced, NaN_Indx

    #--------------------------------------------------------------------------
