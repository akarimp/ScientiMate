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
