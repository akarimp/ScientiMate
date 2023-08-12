def gridgenerator(xmin, xmax, ymin, ymax, gridsize=100, gridsizetype='points', dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2017-12-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.gridgenerator
    =========================

    .. code:: python

        xgrid, ygrid = scientimate.gridgenerator(xmin, xmax, ymin, ymax, gridsize=100, gridsizetype='points', dispout='no')

    Description
    -----------

    Generate 2d x-y grid

    Inputs
    ------

    xmin
        Minimum x of the domain to be generated
    xmax
        Maximum x of the domain to be generated
    ymin
        Minimum y of the domain to be generated
    ymax
        Maximum y of the domain to be generated
    gridsize=100
        | Grid size in x (longitude) and y (latitude) directions to interpolate elevation data on them
        |     if gridsizetype='length' then gridsize is a distance between grid points
        |     if gridsizetype='points' then gridsize is number of grid points in each direction
    gridsizetype='points'
        | Grid size type 
        |     'number': gridsize is considered as number of grid points in each direction
        |     'length': gridsize is considered as length between grid points
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    xgrid
        x of the defined mesh
    ygrid
        y of the defined mesh

    Examples
    --------

    .. code:: python

        import scientimate as sm
        xgrid,ygrid=sm.gridgenerator(0,100,0,100,100,'points','yes')
        xgrid,ygrid=sm.gridgenerator(0,100,0,100,20,'length','yes')

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
    
    #x=type2numpy(x)

    #--------------------------------------------------------------------------
    #Generating grid

    #Calculating grid size in each direction using grid length
    if gridsizetype=='length':
    
        #Checking x and y range to generate correct number of grid points
        if ((xmax-xmin)%gridsize)!=0:
            xmax=xmax-((xmax-xmin)%gridsize)
    
        if ((ymax-ymin)%gridsize)!=0:
            ymax=ymax-((ymax-ymin)%gridsize)
    
        gridsizex=gridsize #Grid size in lat (x) direction
        gridsizey=gridsize #Grid size in lon (y) direction
    
        ngridx=((xmax-xmin)/gridsizex)+1 #Number of grid points in each direction
        ngridy=((ymax-ymin)/gridsizey)+1 #Number of grid points in each direction

        #Generating grid
        xgrid,ygrid=np.meshgrid(np.arange(xmin,xmax+gridsizex,gridsizex),np.arange(ymin,ymax+gridsizey,gridsizey))

    #Calculating grid size in each direction using number grid points 
    elif gridsizetype=='points':
    
        ngridx=gridsize #Number of grid points in each direction
        ngridy=gridsize #Number of grid points in each direction

        gridsizex=(xmax-xmin)/(ngridx-1) #Grid size in x direction
        gridsizey=(ymax-ymin)/(ngridy-1) #Grid size in y direction
        
        #Generating grid
        xgrid,ygrid=np.meshgrid(np.linspace(xmin,xmax,ngridx),np.linspace(ymin,ymax,ngridy))

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
        
        M,N=np.shape(xgrid)
    
        #Plotting data
        for i in range(0,N,1):
            plt.plot([xgrid[0,i],xgrid[0,i]],[np.min(ygrid),np.max(ygrid)])
    

        for i in range(0,M,1):
            plt.plot([np.min(xgrid),np.max(xgrid)],[ygrid[i,0],ygrid[i,0]])

            
        plt.xlabel('x')
        plt.ylabel('y')
        

    #--------------------------------------------------------------------------
    #Outputs
    return xgrid, ygrid

    #--------------------------------------------------------------------------
