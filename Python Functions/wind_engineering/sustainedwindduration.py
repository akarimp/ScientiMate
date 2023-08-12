def sustainedwindduration(windvel, winddir, MaxwindvelVar=2.5, MaxwinddirVar=15, WindTimeInterval=3600, window_length=7, dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2021-01-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.sustainedwindduration
    =================================

    .. code:: python

        SustWindDur = scientimate.sustainedwindduration(windvel, winddir, MaxwindvelVar=2.5, MaxwinddirVar=15, WindTimeInterval=3600, window_length=7, dispout='no')

    Description
    -----------

    Calculate the sustained wind duration

    Inputs
    ------

    windvel
        Wind velocity time series data in (m/s)
    winddir
        Wind direction time series data in (Degree)
    MaxwindvelVar=2.5
        | Maximum allowed wind velocity variation around a mean to allow wind to be considered sustained in (m/s)
        | Coastal Engineering Manual (2015) suggests 2.5 m/s
    MaxwinddirVar=15
        | Maximum allowed wind direction variation around a mean to allow wind to be considered sustained in (Degree)
        | Coastal Engineering Manual (2015) suggests 15 degree
    WindTimeInterval=3600
        | Time interval (time step) between two consecutive wind measurements (data points) in (second)
        | Example: for 60-min measurement interval: WindTimeInterval=3500 (s)
    window_length=7
        Number of data points in sliding window for calculating moving average values, must be odd
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    SustWindDur
        Sustained Wind Duration in (second)

    Examples
    --------

    .. code:: python

        import scientimate as sm

        #Data from https://tidesandcurrents.noaa.gov for Grand Isle, LA, USA (8761724), for June 1st, 2017, reported hourly
        windvel=[3,4.7,4.9,5.3,3.3,3.4,3.3,3.8,3.1,2,1.3,1.2,1.5,3.2,2.9,3,2.9,3.7,3.7,3.1,3.4,2.6,2.5,2.5] #24 Hour wind velocity
        winddir=[78,86,88,107,131,151,163,163,153,150,148,105,105,75,95,103,97,103,108,111,124,183,171,113] #24 Hour wind direction
        SustWindDur=sm.sustainedwindduration(windvel,winddir,2.5,15,3600,7,'yes')

    References
    ----------

    U.S. Army Corps of Engineers (2015). 
    Coastal Engineering Manual. 
    Engineer Manual 1110-2-1100, Washington, D.C.: U.S. Army Corps of Engineers.

    Yamartino, R. J. (1984). 
    A comparison of several "single-pass" estimators of the standard deviation of wind direction. 
    Journal of Climate and Applied Meteorology, 23(9), 1362-1366.

    Change point detection:
    * https://en.wikipedia.org/wiki/Change_detection
    * https://en.wikipedia.org/wiki/Step_detection
    * https://www.mathworks.com/help/signal/ref/findchangepts.html
    * https://www.mathworks.com/help/matlab/ref/ischange.html
    * https://github.com/deepcharles/ruptures

    .. License and Disclaimer
    .. --------------------
    ..
    .. Copyright (c) 2021 Arash Karimpour
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
    #Calculate moving average wind

    #Input size
    if windvel.ndim==1:
        M,=np.shape(windvel)
    elif windvel.ndim==2:
        M,N=np.shape(windvel)

    #Recatngular normal window
    window_length_half=int((window_length-1)/2)
    rect_window_norm=np.ones(window_length)/window_length

    #Moving average wind velocity
    windvel_moveavg = np.convolve(windvel, rect_window_norm, mode='same')
    windvel_moveavg[0:window_length_half]=np.nanmean(windvel[0:window_length_half])
    windvel_moveavg[-window_length_half:]=np.nanmean(windvel[-window_length_half:])

    #Moving average wind velocity difference
    windvel_diff = np.diff(windvel)
    #windvel_grad = np.gradient(windvel,2)

    #Find large variations and replace them with original data
    index_windvel_diff_large=np.nonzero(np.abs(windvel_diff)>=30)[0]
    if np.size(index_windvel_diff_large)!=0:
        for i in range(0,len(index_windvel_diff_large),1):
            index_to_recover=np.arange(index_windvel_diff_large[i]-window_length_half,index_windvel_diff_large[i]+window_length_half,1)
            index_to_recover[index_to_recover<0]=0
            index_to_recover[index_to_recover>M-1]=M-1
            windvel_moveavg[index_to_recover]=windvel[index_to_recover]

    #Moving average wind direction
    winddir_moveavg = np.convolve(winddir, rect_window_norm, mode='same')
    winddir_moveavg[0:window_length_half]=np.nanmean(winddir[0:window_length_half])
    winddir_moveavg[-window_length_half:]=np.nanmean(winddir[-window_length_half:])

    #Moving average wind direction difference
    winddir_diff = np.diff(winddir)
    #winddir_grad = np.gradient(winddir,2)

    #Find large variations and replace them with original data
    index_winddir_diff_large=np.nonzero(np.abs(winddir_diff)>=30)[0]
    if np.size(index_winddir_diff_large)!=0:
        for i in range(0,len(index_winddir_diff_large),1):
            index_to_recover=np.arange(index_winddir_diff_large[i]-window_length_half,index_winddir_diff_large[i]+window_length_half,1)
            index_to_recover[index_to_recover<0]=0
            index_to_recover[index_to_recover>M-1]=M-1
            winddir_moveavg[index_to_recover]=winddir[index_to_recover]


    #--------------------------------------------------------------------------
    #Calculate sustained wind duration

    #Pre-assigning array
    delta_windvel=np.zeros(M) #Pre-assigning array
    delta_winddir=np.zeros(M) #Pre-assigning array
    SustWind_Identifier=np.zeros(M) #Pre-assigning array
    SustWindDur_Counter=np.zeros(M) #Pre-assigning array
    SustWindDur=np.zeros(M) #Pre-assigning array

    #Assigning first element
    SustWind_Identifier[0]=1 #Sustained wind index (1: Sustained, 0: Un-sustained)
    SustWindDur_Counter[0]=SustWind_Identifier[0] #Sustained wind duration index
    SustWindDur[0]=SustWindDur_Counter[0]*WindTimeInterval # Sustained Wind Duration

    #Calculating the sustained wind duration
    k=0 #Sustained wind index counter
    for i in range(1,M,1):
        
        delta_windvel[i]=windvel[i]-windvel_moveavg[i-1]
        delta_winddir[i]=winddir[i]-winddir_moveavg[i-1]
        
        #Check if winddir_moveavg is right before 360 degree and winddir is right after 360 (0) degree
        if ((winddir_moveavg[i-1]>(360-MaxwinddirVar)) and (winddir_moveavg[i-1]<=360 and winddir_moveavg[i]>=0) and (winddir_moveavg[i]<MaxwinddirVar)):
            delta_winddir[i]=(winddir[i]+360)-winddir_moveavg[i-1]
        
        #Check if winddir_moveavg is right after 360 (0) degree and winddir is right before 360 degree
        if ((winddir_moveavg[i-1]>=0) and (winddir_moveavg[i-1]<MaxwinddirVar) and (winddir_moveavg[i]>(360-MaxwinddirVar)) and (winddir_moveavg[i]<=360)):
            delta_winddir[i]=(360-winddir[i])-winddir_moveavg[i-1]
        
        #Calculate sustained wind duration counter, each counter means duration equal to 1 WindTimeInterval
        if ((np.abs(delta_windvel[i])<=MaxwindvelVar) and (np.abs(delta_winddir[i])<=MaxwinddirVar)): #Check for Sustained Wind (1: Sustained, 0: Un-sustained)
            SustWind_Identifier[i]=1 #Check for Sustained Wind (1: Sustained, 0: Un-sustained)
            SustWindDur_Counter[i]=SustWindDur_Counter[i-1]+SustWind_Identifier[i] #Sustained Wind Duration Unit
        else:
            SustWind_Identifier[i]=0 #Check for Sustained Wind (1: Sustained, 0: Un-sustained)
            SustWindDur_Counter[i]=1
            k=i-1


    SustWindDur=SustWindDur_Counter*WindTimeInterval #Sustained Wind Duration
    SustWindDur[SustWindDur_Counter==0]=1*WindTimeInterval #Shortest Sustained Wind Duration is WindTimeInterval seceond
    SustWindDur[SustWindDur<WindTimeInterval]=WindTimeInterval #Shortest Sustained Wind Duration is WindTimeInterval seceond

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
    
        Time=np.arange(1,len(windvel)+1)*WindTimeInterval
    
        plt.subplot(3,1,1)
        plt.plot(Time/3600,windvel,label='Wind Velcity')
        plt.plot(Time/3600,windvel_moveavg,label='Moving Average')
        plt.xlabel('Time (Hour)')
        plt.ylabel('Velocity (m/s)')
        plt.title('Velocity')
        plt.legend()
    
        plt.subplot(3,1,2)
        plt.plot(Time/3600,winddir,label='Wind Direction')
        plt.plot(Time/3600,winddir_moveavg,label='Moving Average')
        plt.xlabel('Time (Hour)')
        plt.ylabel('Direction (Degree)')
        plt.title('Direction')
        plt.legend()
    
        plt.subplot(3,1,3)
        plt.plot(Time/3600,SustWindDur/3600)
        plt.xlabel('Time (Hour)')
        plt.ylabel('Sustained Wind Duration (Hour)')
        plt.title('Sustained Wind Duration')
    

    #--------------------------------------------------------------------------
    #Outputs
    return SustWindDur

    #--------------------------------------------------------------------------
