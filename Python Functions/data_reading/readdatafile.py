def readdatafile(file_name, file_directory=None, column_delimiter='default', skip_beginning_rows=0, skip_beginning_columns=0, nan_values=None):
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

    scientimate.readdatafile
    ========================

    .. code:: python

        data = scientimate.readdatafile(file_name, file_directory=None, column_delimiter='default', skip_beginning_rows=0, skip_beginning_columns=0, nan_values=None)

    Description
    -----------

    Read and extract data from ASCII, text, Comma Separated Values (CSV), Matlab .mat file

    Inputs
    ------

    file_name
        | Name of data file between ' ' mark, example: 'data.csv'
        | Acceptable file types are:
        | ASCII file, example: 'data.xyz'
        |     Files without extension is considered as ASCII
        |     All values in ASCII file except for values in header lines should be number
        | Text file, example: 'data.txt'
        |     All values in text file except for values in header lines should be number
        | Comma Separated Values (CSV) file, example: 'data.csv'
        |     All values in CSV file except for values in header lines and skiped columns should be number
        | Matlab file with .mat extension, example: 'data.mat'
        |     All values in .mat file should be number
    file_directory=None
        Location of data file between ' ' mark, example: r'C:'
    column_delimiter='default'
        | Character (delimiter) that separates column from each other between ' ' mark, example: ','
        | For column_delimiter='default', it uses default delimiter
        | For ASCII or text file, if data are separated by a single space, use ' ' 
        | For ASCII or text file, if data are separated by tab, use '\t' 
        | For ASCII or text file, if data are separated by comma, use ',' 
        | For CSV (comma-separated values) file, use ','
        | For .mat file, it is not required
    skip_beginning_rows=0
        | Number of rows to be skipped at the beginning of the file
        | For skip_beginning_rows=0, no line will be skipped
        | For skip_beginning_rows=1, first line will be skipped
        | For skip_beginning_rows=2, first and second lines will be skipped
        | For skip_beginning_rows=n, first n lines will be skipped
        | Is not applied for .mat file
    skip_beginning_columns=0
        | Number of column to be skipped from the starting column in a file
        | For skip_beginning_columns=0, no column will be skipped
        | For skip_beginning_columns=1, first column will be skipped
        | For skip_beginning_columns=2, first and second columns will be skipped
        | For skip_beginning_columns=n, first n columns will be skipped
        | Is not applied for .mat file
    nan_values=None
        | Value to be considered as NaN
        | None: No value will be considede as NaN
        | 'empty': it replace empty value with NaN
        | Example: nan_values=-99

    Outputs
    -------

    data
        date extracted from input data file

    Examples
    --------

    .. code:: python

        import scientimate as sm

        data=readdatafile('data.csv')

        file_name='data.csv'
        file_directory=None
        column_delimiter='default'
        skip_beginning_rows=0
        skip_beginning_columns=0
        nan_values=None
        data=sm.readdatafile(file_name,file_directory,column_delimiter,skip_beginning_rows,skip_beginning_columns,nan_values)

    References
    ----------


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
    from scipy import io
    import csv 
    import os

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
    #Checking inputs

    if file_directory is None: file_directory=os.getcwd()
    if (skip_beginning_rows<0): skip_beginning_rows=0
    if (skip_beginning_columns<0): skip_beginning_columns=0

    #--------------------------------------------------------------------------
    #Defining file type

    #Extracting file path, file name and file extension
    fname,fextension=os.path.splitext(file_name)

    #Default file type is ASCII
    filetype='ascii'

    #Check for other file types
    #Check if file is Comma Separated Values (CSV) file
    if fextension=='.csv':
        filetype='csv'

    #Check if file is Matlab file with .mat extension
    elif fextension=='.mat':
        filetype='mat'

    #--------------------------------------------------------------------------
    #Reading file

    #Reading all data
    currentFolder=os.getcwd()
    os.chdir(file_directory)
    
    #Reading data from ASCII and text file
    if filetype=='ascii':

        #Importing data
        with open(file_name,'r') as datafile:
            
            #Reading data in datafile
            txtlines=[]
            #txtlinesNotSplit=[]
            for line in datafile:
                
                if column_delimiter=='default':
                    row=line.split()
                else:
                    row=line.split(sep=column_delimiter)

                #Replacing identified value with NaN
                if nan_values is not None:
                    if nan_values=='empty':
                        row=['nan' if not w else w for w in row]
                    else:
                        row=[w.replace(str(nan_values),'nan') for w in row]
                
                txtlines.append(row) #txtlines as a list
                #txtlinesNotSplit.append(line.strip()) #txtlinesNotSplit as a list

        #Removing header lines and converting to number
        txtlines=txtlines[skip_beginning_rows:] #Removing header lines
        data=np.array(txtlines).astype('float') #Converting a list to numpy array
        if skip_beginning_columns!=0:
            data=data[:,skip_beginning_columns:] #Removing skipped columns


    #Reading data from Comma Separated Values (CSV) file
    elif filetype=='csv':

        csvdata=[]
        with open(file_name,'r') as datafile:

            if column_delimiter=='default':
                datareader=csv.reader(datafile)
            else:
                datareader=csv.reader(datafile,delimiter=column_delimiter)

            #csvdata=list(datareader) #csvdata as a list
            for row in datareader:

                #Replacing identified value with NaN
                if nan_values is not None:
                    if nan_values=='empty':
                        row=['nan' if not w else w for w in row]
                    else:
                        row=[w.replace(str(nan_values),'nan') for w in row]

                csvdata.append(row) #csvdata as a list

        #Removing header lines and skipped columns
        csvdata=csvdata[skip_beginning_rows:] #Removing header lines
        if skip_beginning_columns!=0:
            csvdata=[csvdata[i][skip_beginning_columns:] for i in range(0,len(csvdata),1)] #Removing skipped columns

        data=np.array(csvdata).astype('float') #Converting a list to numpy array

    
    #Reading data from Matlab file with .mat extension
    elif filetype=='mat':

        matreader=sp.io.loadmat(file_name)
        matreaderkey=(list(matreader.keys()))[3:]

        #One .mat file saved
        if len(matreaderkey)==1:
            matdata=matreader[matreaderkey[0]]

        #More than one .mat files saved together, adds arrays in 3d axis and returns 3d array
        elif len(matreaderkey)>1:
            for i in range(0,len(matreaderkey),1):
                if i==0:
                    matdata=matreader[matreaderkey[0]]
                else:
                    matdata=np.dstack((matdata,matreader[matreaderkey[i]]))

        #Replacing identified value with NaN
        if nan_values is not None:
            if nan_values=='empty':
                #matdata=['nan' if not x else x for x in matdata]
                pass
            else:
                #matdata=['nan' if x==nan_values else x for x in matdata]
                matdata[matdata==nan_values]=np.nan

        data=np.array(matdata).astype('float') #Converting a list to numpy array
    
        #Converting output array with size of (:,1) to array with size of (:,)
        if np.ndim(data)==2:
            M,N=np.shape(data)
            if N==1:
                data=np.reshape(data,(M)) #Reshaping array

        #Converting output array with size of (:,:,1) to array with size of (:,:)
        elif np.ndim(data)==3:
            M,N,K=np.shape(data)
            if K==1:
                data=np.reshape(data,(M,N)) #Reshaping array

    #--------------------------------------------------------------------------
    #Changing directory to working directory

    os.chdir(currentFolder)

    #--------------------------------------------------------------------------
    #Outputs
    return data

    #--------------------------------------------------------------------------
