def downsamplexyz(x, y, z, RetainRatio=0.5):
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
    
    scientimate.downsamplexyz
    =========================
    
    .. code:: python
    
        x_ds, y_ds, z_ds = scientimate.downsamplexyz(x, y, z, RetainRatio)
    
    Description
    -----------
    
    Downsample x, y,  and z data and retain given ratio
    
    Inputs
    ------
    
    x
        x data
    y
        y data
    z
        z data
    RetainRatio=0.5
        | Define percentage of data to retain, value between 0 and 1
        | Example: RetainRatio=0.8 means 80% of data are retained
    
    Outputs
    -------
    
    x_ds
        Downsample x data
    y_ds
        Downsample y data
    z_ds
        Downsample z data
    
    Examples
    --------
    
    .. code:: python
    
        import scientimate as sm
        import numpy as np
        from numpy import random

        rng = np.random.default_rng()
        x=10*rng.random((1000,1))
        y=10*rng.random((1000,1))
        z=x**2+y**2
        x_ds, y_ds, z_ds=sm.downsamplexyz(x, y, z, 0.7)
    
        rng = np.random.default_rng()
        x=(-90-(-91))*rng.random((1000,1))+(-91)
        y=(31-(30))*rng.random((1000,1))+(30)
        xgrid,ygrid=np.meshgrid(np.linspace(min(x),max(x),1000),np.linspace(min(y),max(y),500))
        zgrid=xgrid**2+ygrid**2
        x_ds, y_ds, z_ds=sm.downsamplexyz(xgrid, ygrid, zgrid, 0.3)
    
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
    z=type2numpy(z)

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
        z_ds=z[dsIndx_dim1] #z data
    
    elif np.ndim(x)==2:
    
        #Data size
        Sz1,Sz2=np.shape(x)
    
        if (Sz1==1 or Sz2==1):

            #Defining index number to down sample data
            dsIndx_dim1=np.int64(np.linspace(0,Sz1-1,int(Sz1*RetainRatio)))
        
            #Retaining down sampled data
            x_ds=x[dsIndx_dim1] #x data
            y_ds=y[dsIndx_dim1] #y data
            z_ds=z[dsIndx_dim1] #z data

        else:

            #Defining index number to down sample data
            dsIndx_dim1=np.int64(np.linspace(0,Sz1-1,int(Sz1*RetainRatio)))
            dsIndx_dim2=np.int64(np.linspace(0,Sz2-1,int(Sz2*RetainRatio)))
        
            #Retaining down sampled data in 1st direction
            x_ds_1=x[dsIndx_dim1,:] #x data
            y_ds_1=y[dsIndx_dim1,:] #y data
            z_ds_1=z[dsIndx_dim1,:] #z data
        
            #Retaining down sampled data in 2nd direction
            x_ds=x_ds_1[:,dsIndx_dim2] #x data
            y_ds=y_ds_1[:,dsIndx_dim2] #y data
            z_ds=z_ds_1[:,dsIndx_dim2] #y data


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
        z_ds_1=z[dsIndx_dim1,:,:] #z data in 1st dim
    
        #Retaining down sampled data in 2nd direction
        x_ds_2=x_ds_1[:,dsIndx_dim2,:] #x data in 2nd dim
        y_ds_2=y_ds_1[:,dsIndx_dim2,:] #y data in 2nd dim
        z_ds_2=z_ds_1[:,dsIndx_dim2,:] #z data in 2nd dim
    
        #Retaining down sampled data in 3rd direction
        x_ds=x_ds_2[:,:,dsIndx_dim3] #x data in 3rd dim
        y_ds=y_ds_2[:,:,dsIndx_dim3] #y data in 3rd dim
        z_ds=z_ds_2[:,:,dsIndx_dim3] #z data in 3rd dim
    
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
        z_ds_1=z[dsIndx_dim1,:,:,:] #z data in 1st dim
    
        #Retaining down sampled data in 2nd direction
        x_ds_2=x_ds_1[:,dsIndx_dim2,:,:] #x data in 2nd dim
        y_ds_2=y_ds_1[:,dsIndx_dim2,:,:] #y data in 2nd dim
        z_ds_2=z_ds_1[:,dsIndx_dim2,:,:] #z data in 2nd dim
    
        #Retaining down sampled data in 3rd direction
        x_ds_3=x_ds_2[:,:,dsIndx_dim3,:] #x data in 3rd dim
        y_ds_3=y_ds_2[:,:,dsIndx_dim3,:] #y data in 3rd dim
        z_ds_3=z_ds_2[:,:,dsIndx_dim3,:] #z data in 3rd dim

        #Retaining down sampled data in 4th direction
        x_ds=x_ds_3[:,:,:,dsIndx_dim4] #x data in 4st dim
        y_ds=y_ds_3[:,:,:,dsIndx_dim4] #y data in 4st dim
        z_ds=z_ds_3[:,:,:,dsIndx_dim4] #z data in 4st dim


    #--------------------------------------------------------------------------
    #Outputs
    return x_ds, y_ds, z_ds

    #--------------------------------------------------------------------------
