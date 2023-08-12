def downsamplex(x, RetainRatio=0.5):
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
    
    scientimate.downsamplex
    =======================
    
    .. code:: python
    
        x_ds = scientimate.downsamplex(x, RetainRatio)
    
    Description
    -----------
    
    Downsample x data and retain given ratio
    
    Inputs
    ------
    
    x
        x data
    RetainRatio=0.5
        | Define percentage of data to retain, value between 0 and 1
        | Example: RetainRatio=0.8 means 80% of data are retained
    
    Outputs
    -------
    
    x_ds
        Downsample x data
    
    Examples
    --------
    
    .. code:: python
    
        import scientimate as sm
        import numpy as np
        from numpy import random

        rng = np.random.default_rng()
        x=10*rng.random((1000,1))
        x_ds=sm.downsamplex(x, 0.7)
    
        rng = np.random.default_rng()
        xgrid=(-90-(-91))*rng.random((1000,500))+(-91)
        x_ds=sm.downsamplex(xgrid, 0.3)

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
    

    elif np.ndim(x)==2:
    
        #Data size
        Sz1,Sz2=np.shape(x)
    
        if (Sz1==1 or Sz2==1):

            #Defining index number to down sample data
            dsIndx_dim1=np.int64(np.linspace(0,Sz1-1,int(Sz1*RetainRatio)))
    
            #Retaining down sampled data
            x_ds=x[dsIndx_dim1] #x data


        else:
            
            #Defining index number to down sample data
            dsIndx_dim1=np.int64(np.linspace(0,Sz1-1,int(Sz1*RetainRatio)))
            dsIndx_dim2=np.int64(np.linspace(0,Sz2-1,int(Sz2*RetainRatio)))
        
            #Retaining down sampled data in 1st direction
            x_ds_1=x[dsIndx_dim1,:] #x data
        
            #Retaining down sampled data in 2nd direction
            x_ds=x_ds_1[:,dsIndx_dim2] #y data
            

    elif np.ndim(x)==3:
    
        #Data size
        Sz1,Sz2,Sz3=np.shape(x)

        #Defining index number to down sample data
        dsIndx_dim1=np.int64(np.linspace(0,Sz1-1,int(Sz1*RetainRatio)))
        dsIndx_dim2=np.int64(np.linspace(0,Sz2-1,int(Sz2*RetainRatio)))
        dsIndx_dim3=np.int64(np.linspace(0,Sz3-1,int(Sz3*RetainRatio)))
    
        #Retaining down sampled data in 1st direction
        x_ds_1=x[dsIndx_dim1,:,:] #x data
    
        #Retaining down sampled data in 2nd direction
        x_ds=x_ds_1[:,dsIndx_dim2,:] #y data

        #Retaining down sampled data in 3rd direction
        x_ds=x_ds_1[:,:,dsIndx_dim3] #z data

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
    
        #Retaining down sampled data in 2nd direction
        x_ds=x_ds_1[:,dsIndx_dim2,:,:] #x data in 2nd dim

        #Retaining down sampled data in 3rd direction
        x_ds=x_ds_1[:,:,dsIndx_dim3,:] #x data in 3rd dim

        #Retaining down sampled data in 4th direction
        x_ds=x_ds_1[:,:,:,dsIndx_dim4] #x data in 4st dim


    #--------------------------------------------------------------------------
    #Outputs
    return x_ds

    #--------------------------------------------------------------------------
