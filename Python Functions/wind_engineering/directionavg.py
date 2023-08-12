def directionavg(direction, NPointsAvg=None, NPointsInterval=None, dispout='no'):
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

    scientimate.directionavg
    ========================

    .. code:: python

        diravg = scientimate.directionavg(direction, NPointsAvg=None, NPointsInterval=None, dispout='no')

    Description
    -----------

    Average direction

    Inputs
    ------

    direction
        Direction time series data in (degree)
    NPointsAvg=length(direction(:,1))
        Number of data points from start of each section (interval) to be averaged
    NPointsInterval=length(direction(:,1))
        Number of points that each section (interval) has
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    diravg
        Averaged direction data in (degree)

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        direction=45*np.random.rand(10)
        diravg=sm.directionavg(direction)

        direction=225*np.random.rand(5*60) #One data point every minute for 5 hours
        diravg=sm.directionavg(direction,10,60,'yes')

    References
    ----------

    Yamartino, R. J. (1984). 
    A comparison of several "single-pass" estimators of the standard deviation of wind direction. 
    Journal of Climate and Applied Meteorology, 23(9), 1362-1366.

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
    
    direction=type2numpy(direction)

    #--------------------------------------------------------------------------
    #Assign default values

    if NPointsAvg is None: NPointsAvg=len(direction)
    if NPointsInterval is None: NPointsInterval=len(direction)

    #--------------------------------------------------------------------------
    #Calculating number of bursts in the input file

    burst=int(len(direction)/NPointsInterval) #Number of burst in input file
    sample=NPointsInterval #Number of data points in 1 burst
    if (NPointsAvg>NPointsInterval): NPointsAvg=NPointsInterval #Checking number of points to be averaged

    #--------------------------------------------------------------------------
    #Direction averaging using Yamartino (1984) method

    #Notes:
    #Use arctan2 instead of arctan to get the averaged wind direction between -180 and 180
    #Add 360 to all numbers to have them all positive
    #Use mod(360) to take care of the ones larger than 360
    
    DirSimpleAvg=np.zeros(burst) #Pre-assigning array
    diravg=np.zeros(burst) #Pre-assigning array
    for i in range(0,burst,1):
        j1=i*sample
        j2=j1+(NPointsAvg-1)
        DirSimpleAvg[i]=np.nanmean(direction[j1:j2+1])
        x=np.cos(np.deg2rad(direction[j1:j2]))
        y=np.sin(np.deg2rad(direction[j1:j2]))
        xmean=np.nanmean(x)
        ymean=np.nanmean(y)
        Dir=np.rad2deg(np.arctan2(ymean,xmean))
        diravg[i]=((Dir+360)%360)

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':

        dt1=np.arange(0,len(direction),1) #Delat_t for all data points 
        dt2=np.arange(0,(burst)*NPointsInterval,NPointsInterval) #Delat_t for averaged data points 
    
        plt.plot(dt1,direction,label='Input Data')
        plt.plot(dt2,diravg,label='Avg Data')
        plt.xlabel('dt (s)')
        plt.ylabel('Direction (Degree)')
        plt.title('Direction')
        plt.legend()

    #--------------------------------------------------------------------------
    #Outputs
    return diravg

    #--------------------------------------------------------------------------
