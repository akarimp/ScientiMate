def downsamplexy(x, y, RetainRatio=0.5):
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
    
    scientimate.downsamplexy
    ========================
    
    .. code:: python
    
        x_ds, y_ds = scientimate.downsamplexy(x, y, RetainRatio)
    
    Description
    -----------
    
    Downsample x and y data and retain given ratio
    
    Inputs
    ------
    
    x
        x data
    y
        y data
    RetainRatio=0.5
        | Define percentage of data to retain, value between 0 and 1
        | Example: RetainRatio=0.8 means 80% of data are retained
    
    Outputs
    -------
    
    x_ds
        Downsample x data
    y_ds
        Downsample y data
    
    Examples
    --------
    
    .. code:: python
    
        import scientimate as sm
        import numpy as np
        from numpy import random

        rng = np.random.default_rng()
        x=10*rng.random((1000,1))
        y=x**2
        x_ds, y_ds=sm.downsamplexy(x, y, 0.7)
    
        rng = np.random.default_rng()
        xgrid=(-90-(-91))*rng.random((1000,500))+(-91)
        ygrid=xgrid**2
        x_ds, y_ds=sm.downsamplexy(xgrid, ygrid, 0.3)
    
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
    #if dispout=='yes':
    #    import matplotlib.pyplot as plt 

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
    #Down sampling the data
    
    #Check RetainRatio
    if ((RetainRatio<=0) or (RetainRatio>1)): RetainRatio=1
    
    if np.ndim(x)==1:
    
        #Data size
        Sz1,=np.shape(x)
    
        #Defining index number to down sample data
        dsIndx_dim1=np.int64(np.linspace(0,Sz1-1,int(Sz1*RetainRatio)))
    
        #Retaining down sampled data
        x_ds=x[dsIndx_dim1] #x data
        y_ds=y[dsIndx_dim1] #y data
    
    elif np.ndim(x)==2:
    
        #Data size
        Sz1,Sz2=np.shape(x)

        if (Sz1==1 or Sz2==1):

            #Defining index number to down sample data
            dsIndx_dim1=np.int64(np.linspace(0,Sz1-1,int(Sz1*RetainRatio)))
        
            #Retaining down sampled data
            x_ds=x[dsIndx_dim1] #x data
            y_ds=y[dsIndx_dim1] #y data

        else:

            #Defining index number to down sample data
            dsIndx_dim1=np.int64(np.linspace(0,Sz1-1,int(Sz1*RetainRatio)))
            dsIndx_dim2=np.int64(np.linspace(0,Sz2-1,int(Sz2*RetainRatio)))
        
            #Retaining down sampled data in 1st direction
            x_ds_1=x[dsIndx_dim1,:] #x data in 1st dim
            y_ds_1=y[dsIndx_dim1,:] #y data in 1st dim
        
            #Retaining down sampled data in 2nd direction
            x_ds=x_ds_1[:,dsIndx_dim2] #x data in 2nd dim
            y_ds=y_ds_1[:,dsIndx_dim2] #y data in 2nd dim
    

    elif np.ndim(x)==3:
    
        #Data size
        Sz1,Sz2,Sz3=np.shape(x)

        #Defining index number to down sample data
        dsIndx_dim1=np.int64(np.linspace(0,Sz1-1,int(Sz1*RetainRatio)))
        dsIndx_dim2=np.int64(np.linspace(0,Sz2-1,int(Sz2*RetainRatio)))
        dsIndx_dim3=np.int64(np.linspace(0,Sz3-1,int(Sz3*RetainRatio)))

        #Retaining down sampled data in 1st direction
        x_ds_1=x[dsIndx_dim1,:,:] #x data in 1st dim
        y_ds_1=y[dsIndx_dim1,:,:] #y data in 1st dim
    
        #Retaining down sampled data in 2nd direction
        x_ds_2=x_ds_1[:,dsIndx_dim2,:] #x data in 2nd dim
        y_ds_2=y_ds_1[:,dsIndx_dim2,:] #y data in 2nd dim
    
        #Retaining down sampled data in 3rd direction
        x_ds=x_ds_2[:,:,dsIndx_dim3] #x data in 3rd dim
        y_ds=y_ds_2[:,:,dsIndx_dim3] #y data in 3rd dim
    
    elif np.ndim(x)==4:
    
        #Data size
        Sz1,Sz2,Sz3,Sz4=np.shape(x)

        #Defining index number to down sample data
        dsIndx_dim1=np.int64(np.linspace(0,Sz1-1,int(Sz1*RetainRatio)))
        dsIndx_dim2=np.int64(np.linspace(0,Sz2-1,int(Sz2*RetainRatio)))
        dsIndx_dim3=np.int64(np.linspace(0,Sz3-1,int(Sz3*RetainRatio)))
        dsIndx_dim4=np.int64(np.linspace(0,Sz4-1,int(Sz4*RetainRatio)))

        #Retaining down sampled data in 1st direction
        x_ds_1=x[dsIndx_dim1,:,:,:] #x data in 1st dim
        y_ds_1=y[dsIndx_dim1,:,:,:] #y data in 1st dim
    
        #Retaining down sampled data in 2nd direction
        x_ds_2=x_ds_1[:,dsIndx_dim2,:,:] #x data in 2nd dim
        y_ds_2=y_ds_1[:,dsIndx_dim2,:,:] #y data in 2nd dim
    
        #Retaining down sampled data in 3rd direction
        x_ds_3=x_ds_2[:,:,dsIndx_dim3,:] #x data in 3rd dim
        y_ds_3=y_ds_2[:,:,dsIndx_dim3,:] #y data in 3rd dim

        #Retaining down sampled data in 4th direction
        x_ds=x_ds_3[:,:,:,dsIndx_dim4] #x data in 4st dim
        y_ds=y_ds_3[:,:,:,dsIndx_dim4] #y data in 4st dim


    #--------------------------------------------------------------------------
    #Outputs
    return x_ds, y_ds

    #--------------------------------------------------------------------------
