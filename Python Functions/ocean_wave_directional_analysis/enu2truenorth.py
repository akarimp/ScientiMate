def enu2truenorth(directionENU, dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2017-05-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.enu2truenorth
    =========================

    .. code:: python

        directionTN = scientimate.enu2truenorth(directionENU, dispout='no')

    Description
    -----------

    Convert mathematical direction (angle) from ENU (East North Up) coordinate system to compass direction with respect to true north

    Inputs
    ------

    directionENU
        | Direction (angle) in ENU (East North Up) coordinate system between 0 and 360 (Degree)
        | If coordinate system is ENU, then x is East and y is North  
    dispout='no'
        | Define to display outputs or not ('yes': display, 'no': not display)
        | Note: inputs can be as a single value or a 1-D vertical array

    Outputs
    -------

    directionTN
        | Direction (angle) in compass direction with respect to true north (Degree)
        | In true north coordinate system, wave comes from as:
        | 0 degree: from north, 90 degree: from east, 180 degree: from south, 270 degree: from west  

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np
        
        directionTN=sm.enu2truenorth(90,'yes')
        
        directionTN=sm.enu2truenorth([15,30,45,60,90],'no')
        
        directionTN=sm.enu2truenorth(np.array([15,30,45,60,90]),'no')

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
    
    directionENU=type2numpy(directionENU)

    #--------------------------------------------------------------------------
    #Converting mathematical angle direction to compass direction with respect to true north

    directionTN=270-directionENU #Converting mathematical direction to compass direction
    directionTN[directionTN<0]=directionTN[directionTN<0]+360
    directionTN[directionTN>360]=directionTN[directionTN>360]-360

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':

        def cart2pol(x,y):
            rho=np.sqrt(x**2+y**2)
            phi = np.arctan2(y,x)
            return(phi,rho)

        def pol2cart(phi,rho):
            x=rho*np.cos(phi)
            y=rho*np.sin(phi)
            return(x,y)

        ENUx,ENUy=pol2cart(np.deg2rad(directionENU),1)
        
        ax=plt.subplot(111, projection='polar')
        arr1=ax.plot([np.deg2rad(directionENU),np.deg2rad(directionENU)],[0,1],label='ENU')
        arr1=ax.plot([np.deg2rad(directionTN),np.deg2rad(directionTN)],[0,1],label='True North')
        #arr1=ax.arrow(np.deg2rad(directionENU),0,0,1,label='ENU')
        #arr1=ax.arrow(np.deg2rad(directionTN),0,0,1,label='True North')
    
        plt.legend()


    #--------------------------------------------------------------------------
    #Outputs
    return directionTN
    
    #--------------------------------------------------------------------------
