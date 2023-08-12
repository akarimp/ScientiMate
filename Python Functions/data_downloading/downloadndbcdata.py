def downloadndbcdata(station_ID, data_type='last5', data_year=2020, data_month=1, filelocation=None):
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
    
    scientimate.downloadndbcdata
    ============================
    
    .. code:: python
    
        station_url, station_real_time_url, station_hist_data_url = scientimate.downloadndbcdata(station_ID, data_type='last5', data_year=2020, data_month=1, filelocation=None)
    
    Description
    -----------
    
    | Download meteorological data from NOAA National Data Buoy Center
    | https://www.ndbc.noaa.gov/
    | To extract data from file, use readtimeseriesfile
    
    Inputs
    ------
    
    station_ID
        | Station ID between ' ' mark, example: '42035'
        | Find stationdata at: https://www.ndbc.noaa.gov/
    data_type='last5'
        | Define which data to download
        | 'last5' : download data of last 5 days
        | 'last45' : download data of last 45 days
        | 'current_year' : download 1 month of data from current year
        | 'historical' : download 1 year of historical data
    data_year=2020
        | Define a year that its data to be downloaded
        | It is required for data_type='current_year' and data_type='historical'
    data_month=1
        | Define a month of a year that its data to be downloaded
        | It is one of the: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12
        | It is required for data_type='current_year'
    filelocation=os.getcwd()
        Location of data file to be downloaded (saved) between ' ' mark, example: r'C:\'
    
    Outputs
    -------
    
    station_url
        Station web address
    station_real_time_url
        Station real time data web address
    station_hist_data_url
        Station historical data web address
    
    Examples
    --------
    
    .. code:: python
    
        import scientimate as sm
        import os

        station_ID='42035'
        data_type='last5'
        station_url, station_real_time_url, station_hist_data_url = sm.downloadndbcdata(station_ID, data_type)

        station_ID='42035'
        data_type='last45'
        station_url, station_real_time_url, station_hist_data_url = sm.downloadndbcdata(station_ID, data_type)

        station_ID='42035'
        data_type='current_year'
        data_year=2020 #Assume current year is 2020
        data_month=1
        filelocation=os.getcwd()
        station_url, station_real_time_url, station_hist_data_url = sm.downloadndbcdata(station_ID, data_type, data_year, data_month, filelocation)

        station_ID='42035'
        data_type='historical'
        data_year=2019
        station_url, station_real_time_url, station_hist_data_url = sm.downloadndbcdata(station_ID, data_type, data_year)

    
    References
    ----------
    
    * https://www.ndbc.noaa.gov
    * https://www.ndbc.noaa.gov/measdes.shtml
    * https://www.ndbc.noaa.gov/docs/ndbc_web_data_guide.pdf
    * https://www.ndbc.noaa.gov/stations.shtml
    * https://www.ndbc.noaa.gov/wstat.shtml
    * https://www.ndbc.noaa.gov/activestations.xml
    * https://www.ndbc.noaa.gov/rt_data_access.shtml
    * https://github.com/nickc1/buoypy
    * https://github.com/fitnr/buoyant
    * https://github.com/GenSci/NDBC
    * https://mhkit-software.github.io/MHKiT/
    * https://github.com/MHKiT-Software/MHKiT-Python
    * https://datalab.marine.rutgers.edu/
    * https://github.com/seagrinch/ooilab
    * https://pecos.readthedocs.io
    
    Measurement Descriptions (https://www.ndbc.noaa.gov/measdes.shtml):
    
    WDIR
        Wind direction (the direction the wind is coming from in degrees clockwise from true N) during the same period used for WSPD. See Wind Averaging Methods
    WSPD
        Wind speed (m/s) averaged over an eight-minute period for buoys and a two-minute period for land stations. Reported Hourly. See Wind Averaging Methods.
    GST 
        Peak 5 or 8 second gust speed (m/s) measured during the eight-minute or two-minute period. The 5 or 8 second period can be determined by payload, See the Sensor Reporting, Sampling, and Accuracy section.
    WVHT
        Significant wave height (meters) is calculated as the average of the highest one-third of all of the wave heights during the 20-minute sampling period. See the Wave Measurements section.
    DPD 
        Dominant wave period (seconds) is the period with the maximum wave energy. See the Wave Measurements section.
    APD 
        Average wave period (seconds) of all waves during the 20-minute period. See the Wave Measurements section.
    MWD 
        The direction from which the waves at the dominant period (DPD) are coming. The units are degrees from true North, increasing clockwise, with North as 0 (zero) degrees and East as 90 degrees. See the Wave Measurements section.
    PRES
        Sea level pressure (hPa). For C-MAN sites and Great Lakes buoys, the recorded pressure is reduced to sea level using the method described in NWS Technical Procedures Bulletin 291 (11/14/80). ( labeled BAR in Historical files)
    ATMP
        Air temperature (Celsius). For sensor heights on buoys, see Hull Descriptions. For sensor heights at C-MAN stations, see C-MAN Sensor Locations
    WTMP
        Sea surface temperature (Celsius). For buoys the depth is referenced to the hull's waterline. For fixed platforms it varies with tide, but is referenced to, or near Mean Lower Low Water (MLLW).
    DEWP
        Dewpoint temperature taken at the same height as the air temperature measurement.
    VIS 
        Station visibility (nautical miles). Note that buoy stations are limited to reports from 0 to 1.6 nmi.
    PTDY
        Pressure Tendency is the direction (plus or minus) and the amount of pressure change (hPa)for a three hour period ending at the time of observation. (not in Historical files)
    TIDE
        The water level in feet above or below Mean Lower Low Water (MLLW).
    
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

    import os
    import urllib.request

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
    #Assign default values

    if filelocation is None: filelocation=os.getcwd()

    #--------------------------------------------------------------------------
    #Change folder to data folder
    Current_Folder=os.getcwd() #Current path
    os.chdir(filelocation)

    #--------------------------------------------------------------------------
    #Station URLs
    station_url = 'https://www.ndbc.noaa.gov/station_page.php?station=' + station_ID
    station_real_time_url = 'https://www.ndbc.noaa.gov/station_realtime.php?station=' + station_ID
    station_hist_data_url = 'https://www.ndbc.noaa.gov/station_history.php?station=' + station_ID
    
    #--------------------------------------------------------------------------
    #Download data
    
    #Download last 5 days data
    if data_type=='last5':
    
        data_url = \
            'https://www.ndbc.noaa.gov/data/5day2/' + \
            station_ID + '_5day.txt'
    
        data_filename = station_ID + '_5day.txt'

    
    #Download last 45 days data
    elif data_type=='last45':
    
        data_url = \
            'https://www.ndbc.noaa.gov/data/realtime2/' + \
            station_ID + '.txt'
    
        data_filename = station_ID + '.txt'

    
    #Download current year data
    elif data_type=='current_year':
    
        month_id = ['1','2','3','4','5','6','7','8','9','a','b','c']
        month_name = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        data_url = \
            'https://www.ndbc.noaa.gov/view_text_file.php?' + \
            'filename=' + station_ID + month_id[data_month-1] + str(data_year) + '.txt.gz' + '&' + \
            'dir=data/stdmet/' + month_name[data_month-1] + '/'
    
        data_filename = station_ID + month_id[data_month-1] + str(data_year) + '.txt'

    
    #Download historical data
    elif data_type=='historical':
    
        data_url = \
            'https://www.ndbc.noaa.gov/view_text_file.php?' + \
            'filename=' + station_ID + 'h' + str(data_year) + '.txt.gz' + '&' + \
            'dir=data/historical/stdmet/'
    
        data_filename = station_ID + 'h' + str(data_year) + '.txt'


    #Save data
    with urllib.request.urlopen(data_url) as f:
        data_bytes = f.read() #'bytes' object
        data_str = data_bytes.decode('utf-8') #Convert to str
        with open(data_filename,'w') as file:
            file.write(data_str)

    #--------------------------------------------------------------------------
    #Return to current folder
    os.chdir(Current_Folder)

    #--------------------------------------------------------------------------
    #Outputs
    return station_url, station_real_time_url, station_hist_data_url

    #--------------------------------------------------------------------------
    