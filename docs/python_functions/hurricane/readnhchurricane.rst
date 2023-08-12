.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-10-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

scientimate.readnhchurricane
============================

.. code:: python

    tsdate, tstime, lattrack, longtrack, MaxSustWindvelSI, MinPressureSI, rMeter, MaxSustWindvelKt, MinPressureMb, rMile, recordid, systemstatus \
        = scientimate.readnhchurricane(filename='hurdat2.txt', filelocation=None, hurricanename='KATRINA', hurricaneyear=2005, hearderlinelen=37, dispout='no')

Description
-----------

Read and extracts hurricane data from National Hurricane Center (NHC) HURDAT2 file

Inputs
------

filename='hurdat2.txt'
    | Name of HURDAT2 file between ' ' mark, example: 'hurdat2.txt'
    | HURDAT2 file can be obtained from www.nhc.noaa.gov/data
    | Columns in HURDAT2 file should be as follow:
    |     column 1: Time series date in YYYYMMDD format
    |     column 2: Time series time in HHMM format
    |     column 3: Record identifier
    |     column 4: Status of system
    |     column 5: Latitude of best track as string in degree
    |     column 6: Longitude of best track as string in degree
    |     column 7: Maximum sustained wind velocity (in knots) 
    |     column 8: Minimum Pressure (in millibars)
    |     column 9: 34 kt (17.49 m/s) wind radii maximum extent in northeastern quadrant (in nautical miles) 
    |     column 10: 34 kt (17.49 m/s) wind radii maximum extent in southeastern quadrant (in nautical miles) 
    |     column 11: 34 kt (17.49 m/s) wind radii maximum extent in southwestern quadrant (in nautical miles) 
    |     column 12: 34 kt (17.49 m/s) wind radii maximum extent in northwestern quadrant (in nautical miles) 
    |     column 13: 50 kt (25.72 m/s) wind radii maximum extent in northeastern quadrant (in nautical miles) 
    |     column 14: 50 kt (25.72 m/s) wind radii maximum extent in southeastern quadrant (in nautical miles) 
    |     column 15: 50 kt (25.72 m/s) wind radii maximum extent in southwestern quadrant (in nautical miles) 
    |     column 16: 50 kt (25.72 m/s) wind radii maximum extent in northwestern quadrant (in nautical miles) 
    |     column 17: 64 kt (32.92 m/s) wind radii maximum extent in northeastern quadrant (in nautical miles) 
    |     column 18: 64 kt (32.92 m/s) wind radii maximum extent in southeastern quadrant (in nautical miles) 
    |     column 19: 64 kt (32.92 m/s) wind radii maximum extent in southwestern quadrant (in nautical miles) 
    |     column 20: 64 kt (32.92 m/s) wind radii maximum extent in northwestern quadrant (in nautical miles) 
filelocation=pwd
    Location of HURDAT2 file between ' ' mark, example: 'C:\'
hurricanename='KATRINA'
    | Hurricane name between ' ' mark, example: 'KATRINA'
    | 'all' will read all data in the file
hurricaneyear=2005
    | Year that hurricane occurs
    | Hurricane Katrina occured in 2005
    | if hurricanename='all'; then hurricaneyear is ignored
hearderlinelen=37
    | Number of charachters in header line
    | In HURDAT2 file, hurricane header line has 37 charachters (including spaces)
    | In HURDAT2 file, hurricane data line has 120 charachters (including spaces)
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

tsdate
    Time series date in YYYYMMDD format
tstime
    Time series time in HHMM format
lattrack
    | Latitude of best track in (Degree)
    | from -90 degree to +90 degree
longtrack
    | Longitude of best track (Degree)
    | from -180 degree to +180 degree
MaxSustWindvelSI
    | Maximum sustained wind velocity in (m/s) 
    | It is defined as the maximum 1-min average wind at an elevation of 10 m
MinPressureSI
    Minimum (central) pressure in (Pa)
rMeter
    | Wind radii maximum in (m)
    | 1st column: 34 kt (17.49 m/s) wind radii maximum extent in northeastern quadrant in (m) 
    | 2nd column: 34 kt (17.49 m/s) wind radii maximum extent in southeastern quadrant in (m) 
    | 3rd column: 34 kt (17.49 m/s) wind radii maximum extent in southwestern quadrant in (m) 
    | 4th column: 34 kt (17.49 m/s) wind radii maximum extent in northwestern quadrant in (m) 
    | 5th column: 50 kt (25.72 m/s) wind radii maximum extent in northeastern quadrant in (m) 
    | 6th column: 50 kt (25.72 m/s) wind radii maximum extent in southeastern quadrant in (m) 
    | 7th column: 50 kt (25.72 m/s) wind radii maximum extent in southwestern quadrant in (m) 
    | 8th column: 50 kt (25.72 m/s) wind radii maximum extent in northwestern quadrant in (m) 
    | 9th column: 64 kt (32.92 m/s) wind radii maximum extent in northeastern quadrant in (m) 
    | 10th column: 64 kt (32.92 m/s) wind radii maximum extent in southeastern quadrant in (m) 
    | 11th column: 64 kt (32.92 m/s) wind radii maximum extent in southwestern quadrant in (m) 
    | 12th column: 64 kt (32.92 m/s) wind radii maximum extent in northwestern quadrant in (m) 
MaxSustWindvelKt
    | Maximum sustained wind velocity in (knots)
    | It is defined as the maximum 1-min average wind at an elevation of 10 m
MinPressureMb
    Minimum (central) pressure in (millibars)
rMile
    | Wind radii maximum in (nautical miles)
    | 1st column: 34 kt (17.49 m/s) wind radii maximum extent in northeastern quadrant in (nautical miles) 
    | 2nd column: 34 kt (17.49 m/s) wind radii maximum extent in southeastern quadrant in (nautical miles) 
    | 3rd column: 34 kt (17.49 m/s) wind radii maximum extent in southwestern quadrant in (nautical miles) 
    | 4th column: 34 kt (17.49 m/s) wind radii maximum extent in northwestern quadrant in (nautical miles) 
    | 5th column: 50 kt (25.72 m/s) wind radii maximum extent in northeastern quadrant in (nautical miles) 
    | 6th column: 50 kt (25.72 m/s) wind radii maximum extent in southeastern quadrant in (nautical miles) 
    | 7th column: 50 kt (25.72 m/s) wind radii maximum extent in southwestern quadrant in (nautical miles) 
    | 8th column: 50 kt (25.72 m/s) wind radii maximum extent in northwestern quadrant in (nautical miles) 
    | 9th column: 64 kt (32.92 m/s) wind radii maximum extent in northeastern quadrant in (nautical miles) 
    | 10th column: 64 kt (32.92 m/s) wind radii maximum extent in southeastern quadrant in (nautical miles) 
    | 11th column: 64 kt (32.92 m/s) wind radii maximum extent in southwestern quadrant in (nautical miles) 
    | 12th column: 64 kt (32.92 m/s) wind radii maximum extent in northwestern quadrant in (nautical miles) 
recordid
    | Record identifier
    |    L: Landfall (center of system crossing a coastline)
    |    P: Minimum in central pressure
    |    I: An intensity peak in terms of both pressure and maximum wind
    |    S: Change of status of the system
    |    T: Provides additional detail on the track (position) of the cyclone
systemstatus
    | Status of system
    |     TD: Tropical cyclone of tropical depression intensity (< 34 knots)
    |     TS: Tropical cyclone of tropical storm intensity (34-63 knots)
    |     HU: Tropical cyclone of hurricane intensity (> 64 knots)
    |     EX: Extratropical cyclone (of any intensity)
    |     SD: Subtropical cyclone of subtropical depression intensity (< 34 knots)
    |     SS: Subtropical cyclone of subtropical storm intensity (> 34 knots)
    |     LO: A low that is neither a tropical cyclone, a subtropical cyclone, nor an extratropical cyclone (of any intensity)
    |     DB: Disturbance (of any intensity) 
    | Note: In all original outputs, missing values are noted as '-999'
    |     In all SI outputs, missing values are noted as 'NaN'

Examples
--------

.. code:: python

    import scientimate as sm
    import os

    filename='hurdat2.txt'
    filelocation=os.getcwd()
    filelocation='C:'
    tsdate,tstime,lattrack,longtrack,MaxSustWindvelSI,MinPressureSI,rMeter,MaxSustWindvelKt,MinPressureMb,rMile,recordid,systemstatus\
        =sm.readnhchurricane(filename,filelocation,'all',2005,37,'yes')

References
----------

Data
* www.nhc.noaa.gov/data/
* www.nhc.noaa.gov/data/hurdat/hurdat2-format-nencpac.pdf
* coast.noaa.gov/hurricanes
* www.aoml.noaa.gov/hrd/data_sub/re_anal.html

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
