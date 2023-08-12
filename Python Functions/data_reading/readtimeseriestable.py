def readtimeseriestable(file_name, file_directory=None, \
    column_delimiter=None, \
    skip_beginning_rows=0, header_row=None, unit_row='no', \
    datetime_column=[0], datetime_format=None, date_interval=None, \
    output_start_row=0, output_end_row=-1, \
    output_start_column=0, output_end_column=-1, output_columns=None, \
    output_start_date=None, output_end_date=None, output_date_interval=None, \
    drop_duplicate='no', find_missing_rows='no', interpolate_nan='no', nan_values=None, \
    index_column=False, \
    dispout='no'):
    """
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

    scientimate.readtimeseriestable
    ===============================

    .. code:: python

        data_DataFrame, data_array = scientimate.readtimeseriestable(file_name, file_directory=None, \
            column_delimiter=None, \
            skip_beginning_rows=0, header_row=0, unit_row='no', \
            datetime_column=[0], datetime_format=None, date_interval=None, \
            output_start_row=0, output_end_row=-1, \
            output_start_column=0, output_end_column=-1, output_columns=None, \
            output_start_date=None, output_end_date=None, output_date_interval=None, \
            drop_duplicate='no', find_missing_rows='no', interpolate_nan='no', nan_values=None, \
            index_column=False, \
            dispout='no')

    Description
    -----------

    Read and extract time-series data from ASCII, text, Comma Separated Values (CSV), and spreadsheet (xlsx, xls, ods) file

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
    datetime_column=[0]
        List of column indexes that contain date and time data
    datetime_format=None
        | Date and time format
        | If datetime_format=None, then date format infer automatically
        | Example: '%Y-%m-%d %H:%M:%S'
    date_interval=None
        Inteval between two consecutive dates
        Example:

            | '30S': 30 seconds interval
            | '40min': 40 minutes interval
            | '40T': 40 minutes interval
            | '12H': 12 hours interval
            | '1D': 1 day interval

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
    output_start_date=None
        | Data will be extracted from output_start_date to output_end_date
        | If None, then it starts from begining of data
        | Example: '2019-01-01 00:00:00'
    output_end_date=None
        | Data will be extracted from output_start_date to output_end_date
        | If None, then it ends at end of data
        | Example: '2019-12-31 23:00:00'
    output_date_interval=None
        Inteval that data are extracted at
        Example:

            | '30S': 30 seconds interval
            | '40min': 40 minutes interval
            | '40T': 40 minutes interval
            | '12H': 12 hours interval
            | '1D': 1 day interval

    drop_duplicate='no'
        | Drop duplicate values 
        | 'yes': it drops duplicate values
        | 'no': it does not drop duplicate values
    find_missing_rows='no'
        | Find missing rows (dates) and insert them as empty rows 
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
        datetime_column = [0,1,2,3,4]
        datetime_format = '%Y %m %d %H %M'
        date_interval = '10min'
        output_start_row = 0
        output_end_row = -1
        output_start_column = 0
        output_end_column = -1
        output_columns = [5,6,8,9,10,17]
        output_start_date = '2019-02-01 00:00:00'
        output_end_date = '2019-03-01 00:00:00'
        output_date_interval = '10min'
        drop_duplicate = 'yes'
        find_missing_rows = 'yes'
        interpolate_nan = 'yes'
        nan_values = [99.00,999,999.0,99.0]
        index_column = False
        dispout = 'yes'
        data_DataFrame, data_array=sm.readtimeseriestable(file_name, file_directory=file_directory, \
            column_delimiter=column_delimiter, \
            skip_beginning_rows=skip_beginning_rows, header_row=header_row, unit_row=unit_row, \
            datetime_column=datetime_column, datetime_format=datetime_format, date_interval=date_interval, \
            output_start_row=output_start_row, output_end_row=output_end_row, \
            output_start_column=output_start_column, output_end_column=output_end_column, output_columns=output_columns, \
            output_start_date=output_start_date, output_end_date=output_end_date, output_date_interval=output_date_interval, \
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
        datetime_column = [0]
        datetime_format = '%Y-%m-%d %H:%M'
        date_interval = '6min'
        output_start_row = 0
        output_end_row = -1
        output_start_column = 0
        output_end_column = -1
        output_columns = [1]
        output_start_date = '2019-02-01 00:00:00'
        output_end_date = '2019-03-01 00:00:00'
        output_date_interval = '6min'
        drop_duplicate = 'yes'
        find_missing_rows = 'yes'
        interpolate_nan = 'yes'
        nan_values = ['NaN']
        index_column = False
        dispout = 'yes'
        data_DataFrame, data_array=sm.readtimeseriestable(file_name, file_directory=file_directory, \
            column_delimiter=column_delimiter, \
            skip_beginning_rows=skip_beginning_rows, header_row=header_row, unit_row=unit_row, \
            datetime_column=datetime_column, datetime_format=datetime_format, date_interval=date_interval, \
            output_start_row=output_start_row, output_end_row=output_end_row, \
            output_start_column=output_start_column, output_end_column=output_end_column, output_columns=output_columns, \
            output_start_date=output_start_date, output_end_date=output_end_date, output_date_interval=output_date_interval, \
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
    """

    #--------------------------------------------------------------------------
    #CODE
    #--------------------------------------------------------------------------
    #Import required packages

    import os
    import numpy as np
    import scipy as sp
    import pandas as pd
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
    
    #x=type2numpy(x)

    #--------------------------------------------------------------------------
    #Check inputs

    if file_directory is None: file_directory=os.getcwd()
    if (skip_beginning_rows<=0): skip_beginning_rows=None
    #if (output_start_row<0): output_start_row=0
    #if (output_start_column<0): output_start_column=0
    #if (output_end_row<0): output_end_row=0
    #if (output_end_column<0): output_end_column=0

    #--------------------------------------------------------------------------
    #Read data

    #Change folder to data folder
    Current_Folder=os.getcwd() #Current path
    os.chdir(file_directory)

    #Get file information
    file_path = os.path.join(file_directory, file_name)
    file_directory = os.path.dirname(file_path)
    file_name = os.path.basename(file_path)
    file_root, file_extension = os.path.splitext(file_path)

    #Read data
    if file_extension=='.xlsx' or file_extension=='.xls' or file_extension=='.ods':
        data = pd.read_excel(file_name, sheet_name=0, header=header_row, index_col=False, skiprows=skip_beginning_rows, \
            na_values=nan_values, \
            parse_dates={'Date_Column' : datetime_column}, keep_date_col=True)
    else:
        data = pd.read_csv(file_name, sep=column_delimiter, header=header_row, index_col=False, skiprows=skip_beginning_rows, \
            na_values=nan_values, \
            parse_dates={'Date_Column' : datetime_column}, keep_date_col=True)

    #Remove unit row
    if unit_row=='yes':
        data = data.drop(0, axis=0)

    #Return to current folder
    os.chdir(Current_Folder)

    #--------------------------------------------------------------------------
    #Set index
    
    if index_column is False:
        #Infer datetime format
        if datetime_format is None:
            InferDatetimeFormat = True #Try to infer the format of the datetime automatically
        else:
            InferDatetimeFormat = False

        #Set 'Date_Column' column as dataframe index
        data['Date_Column'] = pd.to_datetime(data['Date_Column'], format=datetime_format, infer_datetime_format=InferDatetimeFormat)
        data = data.set_index('Date_Column')

    else:
        #Set index_column column as dataframe index
        data = data.set_index(data.iloc[:,index_column])

    #--------------------------------------------------------------------------
    #Drop rows and columns

    if output_start_row<0 or output_start_row>data.shape[0]:
        output_start_row=0

    if output_end_row<=0 or output_end_row>data.shape[0]:
        output_end_row=data.shape[0]
    
    if output_start_column<0 or output_start_column>data.shape[1]:
        output_start_column=0

    if output_end_column<=0 or output_end_column>data.shape[1]:
        output_end_column=data.shape[1]

    #Remove rows
    if output_start_row!=0 or output_end_row!=data.shape[0]:
        if output_end_row==data.shape[0]:
            data = data.iloc[output_start_row:,:]
        elif output_end_row>output_start_row and output_end_row<data.shape[0]:
            data = data.iloc[output_start_row:output_end_row,:]

    #Remove columns
    if output_columns is None:
        if output_start_column!=0 or output_end_column!=data.shape[1]:
            if output_end_column==data.shape[1]:
                data = data.iloc[:,output_start_column:]
            elif output_end_column>output_start_column and output_end_column<data.shape[1]:
                data = data.iloc[:,output_start_column:output_end_column]
    elif output_columns is not None:
        data = data.iloc[:,output_columns]

    #--------------------------------------------------------------------------
    #Sort data by index (date)
    data = data.sort_index()

    #--------------------------------------------------------------------------
    #Drop duplicate index
    if drop_duplicate == 'yes':
        data = data[~data.index.duplicated(keep='first')]
        #Or data = data.drop_duplicates(subset=[index_column], keep='first')

    #--------------------------------------------------------------------------
    #Find missing dates
    if find_missing_rows == 'yes':
        start_date = data.index[0]
        end_date = data.index[-1]
        if date_interval is None:
            date_interval=pd.Timedelta(data.index[1]-data.index[0])
        date_range = pd.date_range(start_date, end_date, freq=date_interval)
        data = data.reindex(date_range)

    #--------------------------------------------------------------------------
    #Extract data for output dates
    
    if (output_start_date is not None) or (output_end_date is not None):

        if output_start_date is None:
            output_start_date=data.index[0]

        if output_end_date is None:
            output_end_date=data.index[-1]

        if output_date_interval is None:
            output_date_interval=pd.Timedelta(data.index[1]-data.index[0])
        
        output_date_range = pd.date_range(output_start_date, output_end_date, freq=output_date_interval)
        
        #data = data.loc[output_date_range,:] #Extract data for defined date range
        data = data.reindex(output_date_range) #Extract data for defined date range
    
    #--------------------------------------------------------------------------
    #Interpolate NaN data
    if interpolate_nan == 'yes':
        data = data.interpolate(method='time')

    #--------------------------------------------------------------------------
    #Prepare output

    data_DataFrame = data.copy()

    #Convert DataFrame to NumPy array
    data_array = data.values
    #Or data_array = data.to_numpy()

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
        
        data.plot()

    #--------------------------------------------------------------------------
    #Outputs
    return data_DataFrame, data_array

    #--------------------------------------------------------------------------
