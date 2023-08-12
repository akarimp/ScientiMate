.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2021-07-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

scientimate.readtimeseriesfile
==============================

.. code:: python

    data_DataFrame, data_array = scientimate.readtimeseriesfile(filename, file_directory=None, \
        column_delimiter=None, header_row=0, skip_rows=None, index_column=False, \
        nan_values=None, interpolate_nan='no', drop_duplicate='no', \
        date_columns=[0], date_format=None, date_interval=None, \
        output_start_date=None, output_end_date=None, \
        output_date_interval=None, output_columns='all', \
        dispout='no')

Description
-----------

Read and extract time-series data from ASCII, text and Comma Separated Values (CSV) file

Inputs
------

filename
    | Name of data file between ' ' mark, example: 'data.csv'
    | Acceptable file types are:
    | ASCII file, example: 'data.xyz'
    | Text file, example: 'data.txt'
    | Comma Separated Values (CSV) file, example: 'data.csv'
file_directory=None
    Location of data file between ' ' mark, example: r'C:'
column_delimiter=None
    | Character (delimiter) that separates column from each other between ' ' mark, example: ','
    | For column_delimiter=None, it uses default delimiter
    | For ASCII or text file, if data are separated by a single space, use ' ' 
    | For ASCII or text file, if data are separated by tab, use '\t' 
    | For ASCII or text file, if data are separated by comma, use ',' 
    | For CSV (comma-separated values) file, use ','
    | Use '\s+' for white space longer than 1
header_row=0
    Index of a row that contains column headers
skip_rows=None
    List of row indexes to be skipped
index_column=False
    Index of a column that contains row labels
nan_values=None
    List of values to be considered as NaN
interpolate_nan='no'
    | Interpolate and replace NaN values 
    | 'yes': it replaces NaN values
    | 'no': it does not replace NaN values
drop_duplicate='no'
    | Drop duplicate values 
    | 'yes': it drops duplicate values
    | 'no': it does not drop duplicate values
date_columns=[0]
    List of column indexes that contain date data
date_format=None
    Format of date
    Example: '%Y-%m-%d %H:%M:%S'
date_interval=None
    Inteval between two consecutive dates
    Example:

        | '30S': 30 seconds interval
        | '40min': 40 minutes interval
        | '40T': 40 minutes interval
        | '12H': 12 hours interval
        | '1D': 1 day interval

output_start_date=None
    Data will be extracted from output_start_date to output_end_date
    If None, then it starts from begining of data
    Example: '2019-01-01 00:00:00'
output_end_date=None
    Data will be extracted from output_start_date to output_end_date
    If None, then it ends at end of data
    Example: '2019-12-31 23:00:00'
output_date_interval=None
    Inteval that data are extracted at
    Example:

        | '30S': 30 seconds interval
        | '40min': 40 minutes interval
        | '40T': 40 minutes interval
        | '12H': 12 hours interval
        | '1D': 1 day interval

output_columns='all'
    List of column indexes to be extracted
    'all': all columns will be extracted
dispout='no'
    | Display outputs
    | 'yes': it displays outputs
    | 'no': it does not display outputs

Outputs
-------

data_DataFrame
    date extracted from input data file as Pandas DataFrame
data_array
    date extracted from input data file as NumPy array

Examples
--------

.. code:: python

    import scientimate as sm
    import os

    #First example
    #Download NDBC data
    station_ID='42035'
    data_type='historical'
    data_year=2019
    data_month=None
    file_directory=os.getcwd()
    station_url, station_real_time_url, station_hist_data_url = sm.downloadndbcdata(station_ID, data_type, data_year, data_month, file_directory)

    #Read NDBC data
    filename = station_ID+'h'+str(data_year)+'.txt'
    file_directory = os.getcwd() #Example: r'C:'
    column_delimiter = '\s+'
    header_row = 0
    skip_rows = [1]
    nan_values = [99.00,999,999.0,99.0]
    interpolate_nan = 'yes'
    drop_duplicate = 'yes'
    date_columns = [0,1,2,3,4]
    date_format = '%Y %m %d %H %M'
    date_interval = '10min'
    output_start_date = '2019-02-01 00:00:00'
    output_end_date = '2019-03-01 00:00:00'
    output_date_interval = '10min'
    output_columns = [5,6,8,9,10,17]
    dispout = 'yes'
    data_DataFrame, data_array=sm.readtimeseriesfile(filename, file_directory=file_directory, \
        column_delimiter=column_delimiter, header_row=header_row, skip_rows=skip_rows, index_column=False, \
        nan_values=nan_values, interpolate_nan=interpolate_nan, drop_duplicate=drop_duplicate, \
        date_columns=date_columns, date_format=date_format, date_interval=date_interval, \
        output_start_date=output_start_date, output_end_date=output_end_date, \
        output_date_interval=output_date_interval, output_columns=output_columns, \
        dispout=dispout)


    #Second example
    #Download Tide and Current data
    station_ID='8761724'
    begin_date='20190201'
    end_date='20190301'
    export_parameter='water_level'
    interval='h'
    units='metric'
    datum='STND'
    time_zone='gmt'
    file_directory=os.getcwd()
    station_information_json, station_datum_json, station_sensors_json, station_url = sm.downloadtidecurrentdata(station_ID, begin_date, end_date, export_parameter, interval, units, datum, time_zone, file_directory)

    #Read Tide and Current data
    filename = 'CO-OPS_'+station_ID+'_'+export_parameter+'_'+begin_date+'_'+end_date+'.csv'
    file_directory = os.getcwd() #Example: r'C:'
    column_delimiter = ','
    header_row = 0
    skip_rows = None
    nan_values = ['NaN']
    interpolate_nan = 'yes'
    drop_duplicate = 'yes'
    date_columns = [0]
    date_format = '%Y-%m-%d %H:%M'
    date_interval = '6min'
    output_start_date = '2019-02-01 00:00:00'
    output_end_date = '2019-03-01 00:00:00'
    output_date_interval = '6min'
    output_columns = [1]
    dispout = 'yes'
    data_DataFrame, data_array=sm.readtimeseriesfile(filename, file_directory=file_directory, \
        column_delimiter=column_delimiter, header_row=header_row, skip_rows=skip_rows, index_column=False, \
        nan_values=nan_values, interpolate_nan=interpolate_nan, drop_duplicate=drop_duplicate, \
        date_columns=date_columns, date_format=date_format, date_interval=date_interval, \
        output_start_date=output_start_date, output_end_date=output_end_date, \
        output_date_interval=output_date_interval, output_columns=output_columns, \
        dispout=dispout)

References
----------

* https://cheatography.com/davechild/cheat-sheets/regular-expressions/
* https://docs.python.org/3/library/datetime.html#strftime-and-strptime-behavior
* https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#timeseries-offset-aliases

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
