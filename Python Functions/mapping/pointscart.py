def pointscart(x1, y1, x2, y2, nmidpoints=1, dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2017-08-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.pointscart
    ======================

    .. code:: python

        xpoints, ypoints, distxy = scientimate.pointscart(x1, y1, x2, y2, nmidpoints=1, dispout='no')

    Description
    -----------

    Generate points between point (x1,y1) and (x2,y2) on cartesian coordinate

    Inputs
    ------

    x1
        x of start point (first point)
    y1
        y of start point (first point)
    x2
        x of end point (last point) 
    y2
        y of end point (last point) 
    nmidpoints=1
        | Number of middle points generated between first point and last point
        | nmidpoints=1;: 1 point between first point and last point is generated
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    xpoint
        | x of points from (x1,y1) and (x2,y2)
        | x1 and x2 are included, total number of points are nmidpoints+2
    ypoint
        | y of points from (x1,y1) and (x2,y2)
        | y1 and y2 are included, total number of points are nmidpoints+2
    distxy
        Distance at each points of (xpoint,ypoint) from the first point of (x1,y1)

    Examples
    --------

    .. code:: python

        import scientimate as sm
        x1=-3
        y1=-3
        x2=3
        y2=3
        xpoints,ypoints,distxy=sm.pointscart(x1,y1,x2,y2,10,'yes')

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

    x1=type2numpy(x1)
    y1=type2numpy(y1)
    x2=type2numpy(x2)
    y2=type2numpy(y2)

    #--------------------------------------------------------------------------
    #Generating points from (x1,y1) to (x2,y2)

    if (nmidpoints<1): nmidpoints=1
    
    xpoints=np.linspace(x1,x2,nmidpoints+2)
    ypoints=np.linspace(y1,y2,nmidpoints+2)

    #--------------------------------------------------------------------------
    #Calculating distance from (x1,y1) to (x2,y2)

    #Distance at each points of (xpoint,ypoint) from the first point of (x1,y1)
    distxy=np.sqrt((xpoints-xpoints[0])**2+(ypoints-ypoints[0])**2) 

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
        
        #Plotting data
        plt.scatter(x1,y1,label='Start Point')
        plt.scatter(x2,y2,label='End Point')
        plt.scatter(xpoints[1:-1],ypoints[1:-1])
        
        plt.xlabel('x')
        plt.ylabel('y')
        plt.legend()

    #--------------------------------------------------------------------------
    #Outputs
    return xpoints, ypoints, distxy

    #--------------------------------------------------------------------------
