def smoothwind(windvel=[], winddir=[], window_length=5, dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2021-02-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    scientimate.smoothwind
    ======================

    .. code:: python

        windvel_smooth, winddir_smooth = scientimate.smoothwind(windvel=[], winddir=[], window_length=5, dispout='no')

    Description
    -----------
    
    Smooth wind data using moving average window
    
    Inputs
    ------
    
    windvel=[]
        | Wind velocity time series data
        | Leave wind velocity empty if only winddir is available
    winddir=[]
        | Wind direction time series data
        | Leave wind direction empty if only windvel is available
    window_length=5
        Number of data points in sliding window for calculating moving average values, must be odd
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)
    
    Outputs
    -------
    
    windvel_smooth
        Smoothed wind velocity
    winddir_smooth
        | Smoothed wind direction
        | Wind directions with large variations are not smoothed
    
    Examples
    --------
    
    .. code:: python

        import scientimate as sm

        #Data from https://tidesandcurrents.noaa.gov for Grand Isle, LA, USA (8761724), for June 1st, 2017, reported hourly
        windvel=[3,4.7,4.9,5.3,3.3,3.4,3.3,3.8,3.1,2,1.3,1.2,1.5,3.2,2.9,3,2.9,3.7,3.7,3.1,3.4,2.6,2.5,2.5] #24 Hour wind velocity
        winddir=[78,86,88,107,131,151,163,163,153,150,148,105,105,75,95,103,97,103,108,111,124,183,171,113] #24 Hour wind direction
        windvel_smooth,winddir_smooth=sm.smoothwind(windvel,winddir,5,'yes')
    
    References
    ----------
    
    U.S. Army Corps of Engineers (2015). 
    Coastal Engineering Manual. 
    Engineer Manual 1110-2-1100, Washington, D.C.: U.S. Army Corps of Engineers.
    
    Yamartino, R. J. (1984). 
    A comparison of several "single-pass" estimators of the standard deviation of wind direction. 
    Journal of Climate and Applied Meteorology, 23(9), 1362-1366.
    
    .. License & Disclaimer
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
    #Define rectangular moving window
    
    #Input size
    if np.size(windvel)!=0:
        if windvel.ndim==1:
            M,=np.shape(windvel)
        elif windvel.ndim==2:
            M,N=np.shape(windvel)
    elif np.size(winddir)!=0:
        if winddir.ndim==1:
            M,=np.shape(winddir)
        elif winddir.ndim==2:
            M,N=np.shape(winddir)
    
    #Recatngular normal window
    window_length_half=int((window_length-1)/2)
    rect_window_norm=np.ones(window_length)/window_length
    
    #--------------------------------------------------------------------------
    #Calculate moving average velocity
    
    #Moving average wind velocity
    if np.size(windvel)!=0:
    
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
    
    else:
        windvel_moveavg=[]
    
    
    #--------------------------------------------------------------------------
    #Calculate moving average direction
    
    #Moving average wind direction
    if np.size(winddir)!=0:
    
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
    
    else:
        winddir_moveavg=[]
    
    
    #--------------------------------------------------------------------------
    #Calculate sustained wind duration
    
    windvel_smooth=windvel_moveavg.copy() #Smoothed wind velocity
    winddir_smooth=winddir_moveavg.copy() #Smoothed wind direction
    
    #--------------------------------------------------------------------------
    #Displaying results
    
    if dispout=='yes':
    
        plt.subplot(2,1,1)
        if np.size(windvel)!=0:
            plt.plot(windvel,label='Wind Velcity')
            plt.plot(windvel_smooth,label='Smoothed Wind Velcity')
            plt.xlabel('Time')
            plt.ylabel('Velocity')
            plt.title('Velocity')
            plt.legend()

    
        plt.subplot(2,1,2)
        if np.size(winddir)!=0:
            plt.plot(winddir,label='Wind Direction')
            plt.plot(winddir_smooth,label='Smoothed Wind Direction')
            plt.xlabel('Time')
            plt.ylabel('Direction')
            plt.title('Direction')
            plt.legend()


    #--------------------------------------------------------------------------
    #Outputs
    return windvel_smooth, winddir_smooth

    #--------------------------------------------------------------------------
    