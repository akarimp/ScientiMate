def windavg(windvel, winddir, NPointsAvg=None, NPointsInterval=None, dispout='no'):
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

    scientimate.windavg
    ===================

    .. code:: python

        windvelavg, winddiravg = scientimate.windavg(windvel, winddir, NPointsAvg=None, NPointsInterval=None, dispout='no')

    Description
    -----------

    Average wind velocity and wind direction

    Inputs
    ------

    windvel
        Wind velocity time series data
    winddir
        Wind direction time series data in (degree)
    NPointsAvg=length(windvel(:,1))
        Number of data points from start of each section (interval) to be averaged
    NPointsInterval=length(windvel(:,1))
        Number of points that each section (interval) has
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    windvelavg
        Averaged wind velocity data
    winddiravg
        Averaged wind direction data in (degree)

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        windvel=10*np.random.rand(10)
        winddir=45*np.random.rand(10)
        windvelavg,winddiravg=sm.windavg(windvel,winddir)

        windvel=10*np.random.rand(5*60) #One data point every minute for 5 hours
        winddir=225*np.random.rand(5*60) #One data point every minute for 5 hours
        windvelavg,winddiravg=sm.windavg(windvel,winddir,10,60,'yes')

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
    
    windvel=type2numpy(windvel)
    winddir=type2numpy(winddir)

    #--------------------------------------------------------------------------
    #Assign default values

    if NPointsAvg is None: NPointsAvg=len(windvel)
    if NPointsInterval is None: NPointsInterval=len(windvel)

    #--------------------------------------------------------------------------
    #Calculating number of bursts in the input file

    burst=int(len(windvel)/NPointsInterval) #Number of burst in input file
    sample=NPointsInterval #Number of data points in 1 burst
    if (NPointsAvg>NPointsInterval): NPointsAvg=NPointsInterval #Checking number of points to be averaged
    
    #--------------------------------------------------------------------------
    #Wind velocity averaging 

    windvelavg=np.zeros(burst) #Pre-assigning array
    for i in range(0,burst,1):
        j1=i*sample
        j2=j1+(NPointsAvg-1)
        windvelavg[i]=np.nanmean(windvel[j1:j2+1])

    #--------------------------------------------------------------------------
    #Wind direction averaging using Yamartino (1984) method

    #Notes:
    #Use arctan2 instead of arctan to get the averaged wind direction between -180 and 180
    #Add 360 to all numbers to have them all positive
    #Use mod(360) to take care of the ones larger than 360
    
    DirSimpleAvg=np.zeros(burst) #Pre-assigning array
    winddiravg=np.zeros(burst) #Pre-assigning array
    for i in range(0,burst,1):
        j1=i*sample
        j2=j1+(NPointsAvg-1)
        DirSimpleAvg[i]=np.nanmean(winddir[j1:j2+1])
        x=np.cos(np.deg2rad(winddir[j1:j2]))
        y=np.sin(np.deg2rad(winddir[j1:j2]))
        xmean=np.nanmean(x)
        ymean=np.nanmean(y)
        Dir=np.rad2deg(np.arctan2(ymean,xmean))
        winddiravg[i]=((Dir+360)%360)

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':

        dt1=np.arange(0,len(windvel),1) #Delat_t for all data points 
        dt2=np.arange(0,(burst)*NPointsInterval,NPointsInterval) #Delat_t for averaged data points 
    
        plt.subplot(2,1,1)
        plt.plot(dt1,windvel,label='Input Data')
        plt.plot(dt2,windvelavg,label='Avg Data')
        plt.xlabel('dt (s)')
        plt.ylabel('Velocity')
        plt.title('Velocity')
        plt.legend()
    
        plt.subplot(2,1,2)
        plt.plot(dt1,winddir,label='Input Data')
        plt.plot(dt2,winddiravg,label='Avg Data')
        plt.xlabel('dt (s)')
        plt.ylabel('Direction (Degree)')
        plt.title('Direction')
        plt.legend()
    
    #--------------------------------------------------------------------------
    #Outputs
    return windvelavg, winddiravg

    #--------------------------------------------------------------------------
