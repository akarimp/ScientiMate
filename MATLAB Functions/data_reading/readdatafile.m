function [data] = readdatafile(file_name, file_directory, column_delimiter, skip_beginning_rows, skip_beginning_columns, nan_values)
%{
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

readdatafile
============

.. code:: MATLAB

    [data] = readdatafile(file_name, column_delimiter, skip_beginning_rows, skip_beginning_columns, nan_values, file_directory)

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
file_directory=pwd;
    Location of data file between ' ' mark, example: 'C:\'
column_delimiter='default';
    | Character (delimiter) that separates column from each other between ' ' mark, example: ','
    | For column_delimiter='default', it uses default delimiter
    | For ASCII or text file, if data are separated by a single space, use ' ' 
    | For ASCII or text file, if data are separated by tab, use '\t' 
    | For ASCII or text file, if data are separated by comma, use ',' 
    | For CSV (comma-separated values) file, use ','
    | For .mat file, it is not required
skip_beginning_rows=0;
    | Number of rows to be skipped at the beginning of the file
    | For skip_beginning_rows=0, no line will be skipped
    | For skip_beginning_rows=1, first line will be skipped
    | For skip_beginning_rows=2, first and second lines will be skipped
    | For skip_beginning_rows=n, first n lines will be skipped
    | Is not applied for .mat file
skip_beginning_columns=0;
    | Number of column to be skipped from the starting column in a file
    | For skip_beginning_columns=0, no column will be skipped
    | For skip_beginning_columns=1, first column will be skipped
    | For skip_beginning_columns=2, first and second columns will be skipped
    | For skip_beginning_columns=n, first n columns will be skipped
    | Is not applied for .mat file
nan_values='none';
    | Value to be considered as NaN
    | 'none': No value will be considede as NaN
    | 'empty': it replace empty value with NaN
    | Example: nan_values=-99

Outputs
-------

data
    date extracted from input data file

Examples
--------

.. code:: MATLAB

    [data]=readdatafile('data.csv');

    file_name='data.csv';
    file_directory=pwd;
    column_delimiter='default';
    skip_beginning_rows=0;
    skip_beginning_columns=0;
    nan_values='none';
    [data]=readdatafile(file_name,file_directory,column_delimiter,skip_beginning_rows,skip_beginning_columns,nan_values);

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
%}

%--------------------------------------------------------------------------
%CODE
%--------------------------------------------------------------------------
%Assign default values

switch nargin
    case 1
        file_directory=pwd; column_delimiter='default'; skip_beginning_rows=0; skip_beginning_columns=0; nan_values='none';
    case 2
        column_delimiter='default'; skip_beginning_rows=0; skip_beginning_columns=0; nan_values='none';
    case 3
        skip_beginning_rows=0; skip_beginning_columns=0; nan_values='none';
    case 4
        skip_beginning_columns=0; nan_values='none';
    case 5
        nan_values='none';
end

%--------------------------------------------------------------------------
%Checking inputs

skip_beginning_rows(skip_beginning_rows<0)=0;
skip_beginning_columns(skip_beginning_columns<0)=0;

%--------------------------------------------------------------------------
%Defining file type

%Extracting file path, file name and file extension
[fpath,fname,fextension]=fileparts(file_name);

%Default file type is ASCII
filetype='ascii';

%Check for other file types
%Check if file is Comma Separated Values (CSV) file
if strcmp(fextension,'.csv')==1
    filetype='csv';

%Check if file is Matlab file with .mat extension
elseif strcmp(fextension,'.mat')==1
    filetype='mat';

end

%--------------------------------------------------------------------------
%Reading file

%Reading all data
currentFolder=pwd;
cd(file_directory)

%Reading data from ASCII and text file
if strcmp(filetype,'ascii')==1

    fid=fopen(file_name);
    txtline=fgetl(fid); %Text line
    txtlines=cell(0,1); %All text lines
    while ischar(txtline)
        txtlines{end+1,1}=txtline;
        txtline=fgetl(fid);
    end
    fclose(fid);

    %Removing header lines and converting to number
    data=[];
    for i=skip_beginning_rows+1:length(txtlines) %Removing header lines
        
        %Reading each line of data
        if strcmp(column_delimiter,'default')==1
            C=textscan(txtlines{i},'');
        else
            C=textscan(txtlines{i},'','delimiter',column_delimiter);
        end
        
        %Replacing identified value with NaN
        if strcmp(nan_values,'none')~=1 
            if strcmp(nan_values,'empty')==1 
                for j=1:length(C)
                    if isempty(C{j})==1
                        C{j}=NaN;
                    end
                end
            else 
                for j=1:length(C)
                    if C{j}==nan_values
                        C{j}=NaN;
                    end
                end
            end    
        end

        %Converting string to number
        data(end+1,:)=cell2mat(C);
    
    end

    if skip_beginning_columns~=0
        data=data(:,skip_beginning_columns+1:end); %Removing skipped columns
    end

%Reading data from Comma Separated Values (CSV) file
elseif strcmp(filetype,'csv')==1

    startrowoffset=skip_beginning_rows; %Starting row offset
    startcolumnoffset=skip_beginning_columns; %Starting column offset
    data=csvread(file_name,startrowoffset,startcolumnoffset);

    %Replacing identified value with NaN
    if strcmp(nan_values,'none')~=1 
        if strcmp(nan_values,'empty')==1
            data(isempty(data)==1)=NaN;
        else
            data(data==nan_values)=NaN;
        end
    end

%Reading data from Matlab file with .mat extension
elseif strcmp(filetype,'mat')==1

    data=importdata(file_name);

    %Replacing identified value with NaN
    if strcmp(nan_values,'none')~=1 
        if strcmp(nan_values,'empty')==1
            data(isempty(data)==1)=NaN;
        else
            data(data==nan_values)=NaN;
        end
    end

end

%--------------------------------------------------------------------------
%Changing directory to working directory

cd(currentFolder)

%--------------------------------------------------------------------------
