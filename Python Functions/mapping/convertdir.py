def convertdir(dirin, CalcMethod='metetotrig'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2017-07-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.convertdir
    =======================

    .. code:: python

        dirout = scientimate.convertdir(dirin, CalcMethod='metetotrig')

    Description
    -----------

    Convert direction from one system to another one

    Inputs
    ------

    dirin
        Direction to be converted in (Degree)
    CalcMethod='metetotrig'
        | Direction conversion method 
        | 'metetotrig': Converting meteorological direction to trigonometric direction, trig=-mete+270;
        | 'trigtomete': Converting trigonometric direction to meteorological direction, mete=270-trig;
        | 'aztomete': Converting azimuth (bearing) to meteorological direction, mete=az+180;
        | 'metetoaz': Converting meteorological direction to azimuth (bearing), az=mete-180;
        | 'aztotrig': Converting azimuth (bearing) to trigonometric direction, trig=-az+90;
        | 'trigtoaz': Converting trigonometric direction to azimuth (bearing), az=90-trig;
        | meteorological direction:
        |     0 (degree): from North, 90 (degree): from East, 180 (degree): from South, 270 (degree): from West 
        | azimuth (bearing) direction which is measured clockwise from the north:
        |     0 (degree): toward North, 90 (degree): toward East, 180 (degree): toward South, 270 (degree): toward West 

    Outputs
    -------

    dirout
        Converted direction in (Degree)

    Examples
    --------

    .. code:: python

        import scientimate as sm

        metedir=[0,90,180,270,360] 
        dirout=sm.convertdir(metedir,'metetotrig')

        azimuth=[0,90,180,270,360] 
        dirout=sm.convertdir(azimuth,'aztomete')

        azimuth=[0,90,180,270,360] 
        dirout=sm.convertdir(azimuth,'aztotrig')

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
    
    dirin=type2numpy(dirin)

    #--------------------------------------------------------------------------
    #Converting direction 

    #Converting meteorological direction to trigonometric direction
    if CalcMethod=='metetotrig':
        
        dirout=-dirin+270
    
    #Converting trigonometric direction to meteorological direction
    elif CalcMethod=='trigtomete':
    
        dirout=270-dirin
    
    #Converting azimuth (bearing) to meteorological direction
    elif CalcMethod=='aztomete':
    
        dirout=dirin+180
    
    #Converting meteorological direction to azimuth (bearing)
    elif CalcMethod=='metetoaz':
    
        dirout=dirin-180
    
    #Converting azimuth (bearing) to trigonometric direction
    elif CalcMethod=='aztotrig':
    
        dirout=-dirin+90
    
    #Converting trigonometric direction to azimuth (bearing)
    elif CalcMethod=='trigtoaz':
    
        dirout=90-dirin
    
    
    #Add 360 to all numbers to have them all positive
    #Use mod(360) to take care of the ones larger than 360
    dirout=((dirout+360)%360) 

    #--------------------------------------------------------------------------
    #Outputs
    return dirout

    #--------------------------------------------------------------------------
