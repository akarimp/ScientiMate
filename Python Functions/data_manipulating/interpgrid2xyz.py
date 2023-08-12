def interpgrid2xyz(xgrid, ygrid, zgrid, x, y, dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2017-10-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.interpgrid2xyz
    ==========================

    .. code:: python

        z = scientimate.interpgrid2xyz(xgrid, ygrid, zgrid, x, y, dispout='no')

    Description
    -----------

    Interpolate 2d gridded data on given scatter point(s) using nearest neighbor method 

    Inputs
    ------

    xgrid
        x data as a [M*N] array
    ygrid
        y data as a [M*N] array
    zgrid
        z data as z(x,y) as a [M*N] array
    x
        x of point that nearest point to that is desired to be found
    y
        y of point that nearest point to that is desired to be found
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    z
        Value of interpolated data at (x,y) as z(x,y)

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        xgrid,ygrid=np.meshgrid(np.linspace(0,10,100),np.linspace(0,10,100))
        zgrid=ygrid*np.sin(xgrid)-xgrid*np.cos(ygrid)
        x=10*np.random.rand(100,1)
        y=10*np.random.rand(100,1)
        z=sm.interpgrid2xyz(xgrid,ygrid,zgrid,x,y,'yes')

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
    
    xgrid=type2numpy(xgrid)
    ygrid=type2numpy(ygrid)
    zgrid=type2numpy(zgrid)
    x=type2numpy(x)
    y=type2numpy(y)

    #--------------------------------------------------------------------------
    #Calculating z(x,y) values

    #Interpolating data into grid using nearest neighbor method
    #First method
    fun=sp.interpolate.NearestNDInterpolator((np.reshape(xgrid,np.size(xgrid)),np.reshape(ygrid,np.size(ygrid))),np.reshape(zgrid,np.size(zgrid)))
    z=fun(x,y)

    #Second method
    #Mapping values of x and y to associated index
    #M,N=np.shape(xgrid)
    #xmap=np.interp(x,np.array([np.nanmin(xgrid),np.nanmax(xgrid)]),np.array([0,N-1]))
    #ymap=np.interp(y,np.array([np.nanmin(ygrid),np.nanmax(ygrid)]),np.array([0,M-1]))
    #xmap[xmap<0]=0
    #xmap[xmap>N-1]=N-1
    #ymap[ymap<0]=0
    #ymap[ymap>M-1]=M-1
    #z=zgrid[np.int64(ymap),np.int64(xmap)]

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
        
        #Plotting data
        fig=plt.figure()
        ax=fig.gca(projection='3d')
        surf=ax.plot_surface(xgrid,ygrid,zgrid,cmap=plt.get_cmap())
        ax.scatter(x,y,z,color='r',label='Input Data')
        
        fig.colorbar(surf)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        fig.legend()


    #--------------------------------------------------------------------------
    #Outputs
    return z

    #--------------------------------------------------------------------------
