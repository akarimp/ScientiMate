def downloadtidecurrentdata(station_ID, begin_date, end_date, export_parameter='water_level', interval='h', units='metric', datum='NAVD', time_zone='gmt', filelocation=None):
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
    
    scientimate.downloadtidecurrentdata
    ===================================
    
    .. code:: python
    
        station_information_json, station_datum_json, station_sensors_json, station_url = scientimate.downloadtidecurrentdata(station_ID, begin_date, end_date, export_parameter='water_level', interval='h', units='metric', datum='NAVD', time_zone='gmt', filelocation=None)
    
    Description
    -----------
    
    | Download meteorological data from NOAA’s Center for Operational Oceanographic Products and Services (CO-OPS)
    | https://tidesandcurrents.noaa.gov/
    | To extract data from file, use readtimeseriesfile
    
    Inputs
    ------
    
    station_ID
        | Station ID between ' ' mark, example: '8761724'
        | Find stationdata at: https://tidesandcurrents.noaa.gov/
    begin_date
        | Start data to download data
        | Date format: 'yyyyMMdd', example: '20210101'
    end_date
        | End data to download data
        | Date format: 'yyyyMMdd', example: '20210131'
    export_parameter='water_level'
        | Data to be exported (measured at the station, depending on availability)
        | 'water_level': water levels
        | 'air_temperature': Air temperature
        | 'water_temperature': Water temperature
        | 'wind': Wind speed, direction, and gusts
        | 'air_pressure': Barometric pressure
        | 'air_gap': Air Gap (distance between a bridge and the water's surface) at the station.
        | 'conductivity': The water's conductivity
        | 'visibility': Visibility from the station's visibility sensor. A measure of atmospheric clarity.
        | 'humidity': Relative humidity
        | 'salinity': Salinity and specific gravity
        | 'hourly_height': Verified hourly height water level
        | 'high_low': Verified high/low water level
        | 'daily_mean': Verified daily mean water level
        | 'monthly_mean': Verified monthly mean water level
        | 'one_minute_water_level': One minute water level
        | 'predictions': 6 minute predictions water level*
        | 'datums': datums data for the stations
        | 'currents': Currents data
        | 'currents_predictions': Currents predictions data
    interval='h'
        | Interval of reported data
        | '6': 6-min interval data are returned
        | 'h': Hourly data are returned
    units='metric'
        | Data unit system
        | 'metric': Metric (Celsius, meters, cm/s) units
        | 'english': English (fahrenheit, feet, knots) units
    datum='NAVD'
        | Elevation datum
        | 'CRD': Columbia River Datum
        | 'IGLD': International Great Lakes Datum
        | 'LWD': Great Lakes Low Water Datum (Chart Datum)
        | 'MHHW': Mean Higher High Water
        | 'MHW': Mean High Water
        | 'MTL': Mean Tide Level
        | 'MSL': Mean Sea Level
        | 'MLW': Mean Low Water
        | 'MLLW': Mean Lower Low Water
        | 'NAVD': North American Vertical Datum
        | 'STND': Station Datum
    time_zone='gmt'
        | Time zone
        | 'gmt': Greenwich Mean Time
        | 'lst': Local Standard Time. The time local to the requested station
        | 'lst_ldt': Local Standard/Local Daylight Time. The time local to the requested station
    filelocation=os.getcwd()
        Location of data file to be downloaded (saved) between ' ' mark, example: r'C:\'

    Outputs
    -------
    
    station_information_json
        Station information
    station_datum_json
        Station datum information
    station_sensors_json
        Station sensor information
    station_url
        Station web address
    
    Examples
    --------
    
    .. code:: python
    
        import scientimate as sm
        import os

        station_ID='8761724'
        begin_date='20210101'
        end_date='20210131'
        export_parameter='water_level'
        interval='h'
        units='metric'
        datum='STND'
        time_zone='gmt'
        filelocation=os.getcwd()
        station_information_json, station_datum_json, station_sensors_json, station_url = sm.downloadtidecurrentdata(station_ID, begin_date, end_date, export_parameter, interval, units, datum, time_zone, filelocation)
    
    References
    ----------
    
    * https://tidesandcurrents.noaa.gov/
    * https://tidesandcurrents.noaa.gov/web_services_info.html
    * https://tidesandcurrents.noaa.gov/api-helper/url-generator.html
    * https://api.tidesandcurrents.noaa.gov/api/prod/
    * https://tidesandcurrents.noaa.gov/datum_options.htmlv
    * https://github.com/GClunies/noaa_coops
    
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
    import json
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
    station_url = 'https://tidesandcurrents.noaa.gov/stationhome.html?id=' + station_ID
    station_water_level_url = 'https://tidesandcurrents.noaa.gov/waterlevels.html?id=' + station_ID
    station_meteorological_url = 'https://tidesandcurrents.noaa.gov/met.html?id=' + station_ID
    station_physical_oceanography_url = 'https://tidesandcurrents.noaa.gov/physocean.html?id=' + station_ID

    #--------------------------------------------------------------------------
    #Station information in json format
    station_information_json_url = 'https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/' + station_ID + '/details.json?units=metric'
    station_datum_json_url = 'https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/' + station_ID + '/datums.json?units=metric'
    station_sensors_json_url = 'https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/' + station_ID + '/sensors.json?units=metric'

    with urllib.request.urlopen(station_information_json_url) as f:
        station_information_json = json.load(f)
    
    with urllib.request.urlopen(station_datum_json_url) as f:
        station_datum_json = json.load(f)

    with urllib.request.urlopen(station_sensors_json_url) as f:
        station_sensors_json = json.load(f)

    #--------------------------------------------------------------------------
    #Download data
    
    output_format = 'csv'

    data_url = \
        'https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?' + \
        'begin_date=' + begin_date + '&' + \
        'end_date=' + end_date + '&' + \
        'station=' + station_ID + '&' + \
        'product=' + export_parameter + '&' + \
        'datum=' + datum + '&' + \
        'time_zone=' + time_zone + '&' + \
        'units=' + units + '&' + \
        'format=' + output_format + '&' + \
        'interval=' + interval

    data_filename = 'CO-OPS_' + station_ID + '_' + export_parameter + '_' + begin_date + '_' + end_date + '.csv'

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
    return station_information_json, station_datum_json, station_sensors_json, station_url

    #--------------------------------------------------------------------------
