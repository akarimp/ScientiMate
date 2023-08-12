def endpointcart(x1, y1, linelength=1, lineangle=0, dispout='no'):
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

    scientimate.endpointcart
    ========================

    .. code:: python

        x2, y2, lineslope = scientimate.endpointcart(x1, y1, linelength=1, lineangle=0, dispout='no')

    Description
    -----------

    Find an end point of the straight line segment from its starting point (x1,y1) and its angle on cartesian coordinate

    Inputs
    ------

    x1
        x of start point (first point)
    y1
        y of start point (first point)
    linelength=1
        Length of a line segment
    lineangle=0
        Angle of a line segment from start point (first point) toward end point (last point) in (Degree)
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    x2
        x of end point (last point) 
    y2
        y of end point (last point) 
    lineslope
        | Slope of a line segment from start point (first point) toward end point (last point) in (Degree)
        | lineslope=tan(deg2rad(lineangle))

    Examples
    --------

    .. code:: python

        import scientimate as sm

        x1=3
        y1=3
        linelength=5
        lineangle=45
        x2,y2,lineslope=sm.endpointcart(x1,y1,linelength,lineangle,'yes')

        x1=[0,0,0,0,0,0,0,0]
        y1=[0,0,0,0,0,0,0,0]
        linelength=5
        lineangle=[0,45,90,135,180,225,270,315]
        x2,y2,lineslope=sm.endpointcart(x1,y1,linelength,lineangle,'no')

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
    linelength=type2numpy(linelength)
    lineangle=type2numpy(lineangle)

    #--------------------------------------------------------------------------
    #Calculating end point (x2,y2) from start point (x1,y1)

    #x at an end of a line segment starts from (x1,y1) with length of linelength and angle of lineangle
    x2=x1+linelength*np.cos(np.deg2rad(lineangle)) 
    
    #y at an end of a line segment starts from (x1,y1) with length of linelength and angle of lineangle
    y2=y1+linelength*np.sin(np.deg2rad(lineangle))
    
    # Slope of a line segment from (x1,y1) toward (x2,y2) in (Degree)
    lineslope=np.tan(np.deg2rad(lineangle))

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
        
        #Plotting data
        plt.scatter(x1,y1,label='Start Point')
        plt.scatter(x2,y2,label='End Point')
        for i in range(0,len(x1),1):
            plt.plot([x1[i],x2[i]],[y1[i],y2[i]])
        
        plt.xlabel('x')
        plt.ylabel('y')
        plt.legend()
        

    #--------------------------------------------------------------------------
    #Outputs
    return x2, y2, lineslope

    #--------------------------------------------------------------------------
