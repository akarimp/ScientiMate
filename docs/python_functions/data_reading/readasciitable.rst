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

scientimate.readasciitable
==========================

.. code:: python

    data_DataFrame, data_array = scientimate.readasciitable(file_name, file_directory=None, \
        column_delimiter=None, \
        skip_beginning_rows=0, header_row=None, unit_row='no', \
        output_start_row=0, output_end_row=-1, \
        output_start_column=0, output_end_column=-1, output_columns=None, \
        drop_duplicate='no', find_missing_rows='no', interpolate_nan='no', nan_values=None, \
        index_column=False, \
        dispout='no')

Description
-----------

Read and extract data from ASCII, text, Comma Separated Values (CSV), and spreadsheet (xlsx, xls, ods) file

Inputs
------

file_name
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
skip_beginning_rows=0
    Number of rows to be skipped at the beginning of the file
header_row=None
    | Index of a row that contains column headers
    | Example: header_row=0 means 1st line is a header row
unit_row='no'
    Define if there is a unit row after a header row or not
    | 'yes': there is a unit row after header row
    | 'no': there is not a unit row after header row
output_start_row=0
    | Index of first row to read
    | Rows will be extracted from output_start_row to output_end_row
output_end_row=-1
    | Index of last row to read
    | Rows will be extracted from output_start_row to output_end_row
output_start_column=0
    | Index of first column to read
    | Columns will be extracted from output_start_column to output_end_column
output_end_column=-1
    | Index of last column to read
    | Columns will be extracted from output_start_column to output_end_column
output_columns=None
    | List of column indexes to be extracted
    | None: columns will be extracted from output_start_column to output_end_column
    | If output_columns is not None then output_start_column and output_end_column are ignored
    | Example: output_columns=[0,1] means 1st and 2nd columns are returned
drop_duplicate='no'
    | Drop duplicate values 
    | 'yes': it drops duplicate values
    | 'no': it does not drop duplicate values
find_missing_rows='no'
    | Find missing rows (indexes) and insert them as empty rows 
    | 'yes': it find missing dates
    | 'no': it does not find missing dates
interpolate_nan='no'
    | Interpolate and replace NaN values 
    | 'yes': it replaces NaN values
    | 'no': it does not replace NaN values
nan_values=None
    List of values to be considered as NaN
index_column=False
    Index of a column that contains row labels
dispout='no'
    | Display outputs
    | 'yes': it displays outputs
    | 'no': it does not display outputs

Outputs
-------

data_DataFrame
    data extracted from input data file as Pandas DataFrame
data_array
    data extracted from input data file as NumPy array

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
    file_name = station_ID+'h'+str(data_year)+'.txt'
    file_directory = os.getcwd() #Example: r'C:'
    column_delimiter = '\s+'
    skip_beginning_rows = 0
    header_row = 0
    unit_row = 'yes'
    output_start_row = 0
    output_end_row = -1
    output_start_column = 0
    output_end_column = -1
    output_columns = [5,6,8,9,10,17]
    drop_duplicate = 'yes'
    find_missing_rows = 'yes'
    interpolate_nan = 'yes'
    nan_values = [99.00,999,999.0,99.0]
    index_column = False
    dispout = 'yes'
    data_DataFrame, data_array=sm.readasciitable(file_name, file_directory=file_directory, \
        column_delimiter=column_delimiter, \
        skip_beginning_rows=skip_beginning_rows, header_row=header_row, unit_row=unit_row, \
        output_start_row=output_start_row, output_end_row=output_end_row, \
        output_start_column=output_start_column, output_end_column=output_end_column, output_columns=output_columns, \
        drop_duplicate=drop_duplicate, find_missing_rows=find_missing_rows, interpolate_nan=interpolate_nan, nan_values=nan_values, \
        index_column=index_column, \
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
    file_name = 'CO-OPS_'+station_ID+'_'+export_parameter+'_'+begin_date+'_'+end_date+'.csv'
    file_directory = os.getcwd() #Example: r'C:'
    column_delimiter = ','
    skip_beginning_rows = 0
    header_row = 0
    unit_row = 'no'
    output_start_row = 0
    output_end_row = -1
    output_start_column = 0
    output_end_column = -1
    output_columns = [1]
    drop_duplicate = 'yes'
    find_missing_rows = 'yes'
    interpolate_nan = 'yes'
    nan_values = ['NaN']
    index_column = False
    dispout = 'yes'
    data_DataFrame, data_array=sm.readasciitable(file_name, file_directory=file_directory, \
        column_delimiter=column_delimiter, \
        skip_beginning_rows=skip_beginning_rows, header_row=header_row, unit_row=unit_row, \
        output_start_row=output_start_row, output_end_row=output_end_row, \
        output_start_column=output_start_column, output_end_column=output_end_column, output_columns=output_columns, \
        drop_duplicate=drop_duplicate, find_missing_rows=find_missing_rows, interpolate_nan=interpolate_nan, nan_values=nan_values, \
        index_column=index_column, \
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
