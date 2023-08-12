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
