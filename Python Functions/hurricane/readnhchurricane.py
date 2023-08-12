def readnhchurricane(filename='hurdat2.txt', filelocation=None, hurricanename='KATRINA', hurricaneyear=2005, hearderlinelen=37, dispout='no'):
    """
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
    """

    #--------------------------------------------------------------------------
    #CODE
    #--------------------------------------------------------------------------
    #Import required packages

    import os
    import re
    import numpy as np
    import scipy as sp
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
    #Assigning default values

    if filelocation is None: filelocation=os.getcwd()

    #--------------------------------------------------------------------------
    #Reading xyz file

    #Reading all data
    currentFolder=os.getcwd()
    os.chdir(filelocation)

    #Reading HURDAT2 data into cell array
    with open(filename,'r') as datafile:
        
        #Reading data in datafile, first method
        #datalist=list(datafile)
        ##datalist=datafile.readlines()
        #txtlines=[]
        #txtlinesNotSplit=[]
        #for i in range(0,len(datalist),1):
        #    #row=datalist[i].strip()
        #    row=datalist[i].split()
        #    txtlines.append(row)
        #    txtlinesNotSplit.append(datalist[i].strip())
           
        #Reading data in datafile, second method
        txtlines=[]
        txtlinesNotSplit=[]
        for line in datafile:
            row=line.split()
            row=[w.replace(',','') for w in row] #Removing comma from each element in the list
            txtlines.append(row)
            txtlinesNotSplit.append(line.strip())
        
        #Converting a list to numpy array
        #txtlines=np.array(txtlines).astype('float')

    #--------------------------------------------------------------------------
    #Defining hearder lines and data lines for each hurricane

    #Reading data of all hurricanes in the file
    if hurricanename=='all':
    
        #Reading data 
        data=[]
        for i in range (0,len(txtlines),1):
            if len(txtlinesNotSplit[i])!=hearderlinelen: #Checking that txtlines{i} is not header line
                data.append(txtlines[i])
    
    #Reading data for hurricane with specific name in hurricanename
    else:
    
        #Obtaining length of each line in text cell array
        lentxtlines=np.zeros(len(txtlines)) #Pre-assigning array
        for i in range (0,len(txtlines),1):
            lentxtlines[i]=len(txtlinesNotSplit[i])
    

        #Defining index of header lines
        #Line associated to hurricane headline has 37 charachters
        #Line associated to hurricane data has 120 charachters
        IndxHurrHeader=np.int_((np.nonzero(lentxtlines==hearderlinelen))[0])
    
        #Defining header lines contains a desired hurricane name, which might be more than one
        #For example hurricane with name Katrina occured 7 times between 1967 to 2005
        IndxHurrName=[]
        for i in range(0,len(txtlines),1):
            if hurricanename in txtlinesNotSplit[i]:
                IndxHurrName.append(i)

        #Defining lines contains a desired hurricane year
        IndxHurrYear=[]
        for i in range(0,len(txtlines),1):
            if str(hurricaneyear) in txtlinesNotSplit[i]:
                IndxHurrYear.append(i)
    
        #Defining a header line associated to a desired hurricane name in that specific year
        #IndxHurrHeaderReq=list(set(IndxHurrYear).intersection(set(IndxHurrName)))
        #IndxHurrHeaderReq=list(set(IndxHurrYear) & set(IndxHurrName))
        IndxHurrHeaderReq=list(set.intersection(set(IndxHurrYear) and set(IndxHurrName)))
        IndxHurrHeaderReq=int(IndxHurrHeaderReq[0]) #Converting to int
    
        #Defining index of a desired hurricane header line in IndxHurrHeader
        IndxinIndxHurrHeader=np.int_((np.nonzero(IndxHurrHeader==IndxHurrHeaderReq))[0])
    
        #Defining first line and last line containing data of a desired hurricane
        #DataStartIndx=[x+1 for x in IndxHurrHeaderReq] #Index of first line containing data of a desired hurricane
        DataStartIndx=IndxHurrHeaderReq+1 #Index of first line containing data of a desired hurricane
        if IndxHurrHeader[IndxinIndxHurrHeader]<IndxHurrHeader[-1]:
            DataEndIndx=int(IndxHurrHeader[IndxinIndxHurrHeader+1]-1) #Index of last line containing data of a desired hurricane
        else: #If hurricane is the last hurricane in HURDAT2 file
            DataEndIndx=len(txtlines) #Index of last line containing data of a desired hurricane
    
        #Reading data 
        data=[]
        for i in range (DataStartIndx,DataEndIndx+1,1):
            data.append(txtlines[i])
    
    
    #Extracting data
    tsdate=[]
    tstime=[]
    recordid=[]
    systemstatus=[]
    latbtStr=[]
    longbtStr=[]
    MaxSustWindvelKt=[]
    MinPressureMb=[]
    rNE34ktMile=[]
    rSE34ktMile=[]
    rSW34ktMile=[]
    rNW34ktMile=[]
    rNE50ktMile=[]
    rSE50ktMile=[]
    rSW50ktMile=[]
    rNW50ktMile=[]
    rNE64ktMile=[]
    rSE64ktMile=[]
    rSW64ktMile=[]
    rNW64ktMile=[]
    for i in range(0,len(data),1):
        tsdate.append(data[i][0]) #Time series date in YYYYMMDD format
        tstime.append(data[i][1]) #Time series time in HHMM format
        recordid.append(data[i][2]) #Record identifier
        systemstatus.append(data[i][3]) #Status of system
        latbtStr.append(data[i][4]) #Latitude of best track as string in degree
        longbtStr.append(data[i][5]) #Longitude of best track as string in degree
        MaxSustWindvelKt.append(data[i][6]) #Maximum sustained wind velocity (in knots) 
        MinPressureMb.append(data[i][7]) #Minimum Pressure (in millibars)
        rNE34ktMile.append(data[i][8]) #34 kt (17.49 m/s) wind radii maximum extent in northeastern quadrant (in nautical miles) 
        rSE34ktMile.append(data[i][9]) #34 kt (17.49 m/s) wind radii maximum extent in southeastern quadrant (in nautical miles) 
        rSW34ktMile.append(data[i][10]) #34 kt (17.49 m/s) wind radii maximum extent in southwestern quadrant (in nautical miles) 
        rNW34ktMile.append(data[i][11]) #34 kt (17.49 m/s) wind radii maximum extent in northwestern quadrant (in nautical miles) 
        rNE50ktMile.append(data[i][12]) #50 kt (25.72 m/s) wind radii maximum extent in northeastern quadrant (in nautical miles) 
        rSE50ktMile.append(data[i][13]) #50 kt (25.72 m/s) wind radii maximum extent in southeastern quadrant (in nautical miles) 
        rSW50ktMile.append(data[i][14]) #50 kt (25.72 m/s) wind radii maximum extent in southwestern quadrant (in nautical miles) 
        rNW50ktMile.append(data[i][15]) #50 kt (25.72 m/s) wind radii maximum extent in northwestern quadrant (in nautical miles) 
        rNE64ktMile.append(data[i][16]) #64 kt (32.92 m/s) wind radii maximum extent in northeastern quadrant (in nautical miles) 
        rSE64ktMile.append(data[i][17]) #64 kt (32.92 m/s) wind radii maximum extent in southeastern quadrant (in nautical miles) 
        rSW64ktMile.append(data[i][18]) #64 kt (32.92 m/s) wind radii maximum extent in southwestern quadrant (in nautical miles) 
        rNW64ktMile.append(data[i][19]) #64 kt (32.92 m/s) wind radii maximum extent in northwestern quadrant (in nautical miles) 
    
    #Separating best track latitude and longitude values from N, S, E and W text
    latbtNumStr=[]
    longbtNumStr=[]
    for i in range(0,len(data),1):
        
        #Split and only keeps number and ingores strings
        #latbtNumStr.append(re.findall(r'[-+]?\d*\.\d+|\d+',latbtStr[i]))
        #longbtNumStr.append(re.findall(r'[-+]?\d*\.\d+|\d+',longbtStr[i]))
        
        #Split and keeps both number and strings
        latbtNumStr.append(re.findall(r'([-+]?\d*\.\d+|\d+)([a-zA-Z]+)',latbtStr[i]))
        longbtNumStr.append(re.findall(r'([-+]?\d*\.\d+|\d+)([a-zA-Z]+)',longbtStr[i]))

        #Note:Parentheses around each group makes to include both ones that match and the ones that do not match
        #Two other examples that split and keep both number and strings
        #latbtNumStr.append(re.split('([a-zA-Z]+)',latbtStr[i]))
        #latbtNumStr.append(re.split(r'([-+]?\d*\.\d+|\d+)',latbtStr[i]))
    
    #Converting best track latitude and longitude from N, S, E and W format to -90 to 90 and -180 to 180
    lattrack=[]
    longtrack=[]
    for i in range(0,len(latbtNumStr),1):
    
        if list(latbtNumStr[i])[0][1]=='N':
            lattrack.append([float(x) for x in [(latbtNumStr[i])[0][0]]][0])
        elif list(latbtNumStr[i])[0][1]=='S':
            lattrack.append([-float(x) for x in [(latbtNumStr[i])[0][0]]][0])

    
        if list(longbtNumStr[i])[0][1]=='E':
            longtrack.append([float(x) for x in [(longbtNumStr[i])[0][0]]][0])
        elif list(longbtNumStr[i])[0][1]=='W':
            longtrack.append([-float(x) for x in [(longbtNumStr[i])[0][0]]][0])
    
    
    #Converting string to float
    MaxSustWindvelKt=[float(x) for x in MaxSustWindvelKt] #Maximum sustained wind velocity (in knots) 
    MinPressureMb=[float(x) for x in MinPressureMb] #Minimum Pressure (in millibars)
    rNE34ktMile=[float(x) for x in rNE34ktMile] #34 kt (17.49 m/s) wind radii maximum extent in northeastern quadrant (in nautical miles) 
    rSE34ktMile=[float(x) for x in rSE34ktMile] #34 kt (17.49 m/s) wind radii maximum extent in southeastern quadrant (in nautical miles) 
    rSW34ktMile=[float(x) for x in rSW34ktMile] #34 kt (17.49 m/s) wind radii maximum extent in southwestern quadrant (in nautical miles) 
    rNW34ktMile=[float(x) for x in rNW34ktMile] #34 kt (17.49 m/s) wind radii maximum extent in northwestern quadrant (in nautical miles) 
    rNE50ktMile=[float(x) for x in rNE50ktMile] #50 kt (25.72 m/s) wind radii maximum extent in northeastern quadrant (in nautical miles) 
    rSE50ktMile=[float(x) for x in rSE50ktMile] #50 kt (25.72 m/s) wind radii maximum extent in southeastern quadrant (in nautical miles) 
    rSW50ktMile=[float(x) for x in rSW50ktMile] #50 kt (25.72 m/s) wind radii maximum extent in southwestern quadrant (in nautical miles) 
    rNW50ktMile=[float(x) for x in rNW50ktMile] #50 kt (25.72 m/s) wind radii maximum extent in northwestern quadrant (in nautical miles) 
    rNE64ktMile=[float(x) for x in rNE64ktMile] #64 kt (32.92 m/s) wind radii maximum extent in northeastern quadrant (in nautical miles) 
    rSE64ktMile=[float(x) for x in rSE64ktMile] #64 kt (32.92 m/s) wind radii maximum extent in southeastern quadrant (in nautical miles) 
    rSW64ktMile=[float(x) for x in rSW64ktMile] #64 kt (32.92 m/s) wind radii maximum extent in southwestern quadrant (in nautical miles) 
    rNW64ktMile=[float(x) for x in rNW64ktMile] #64 kt (32.92 m/s) wind radii maximum extent in northwestern quadrant (in nautical miles) 

    #Converting values to SI Units
    MaxSustWindvelSI=[float(x)*0.514444 for x in MaxSustWindvelKt] #Maximum sustained wind velocity in (m/s) 
    MinPressureSI=[float(x)*100 for x in MinPressureMb] #Minimum Pressure in (Pa)
    rNE34ktSI=[float(x)*1852 for x in rNE34ktMile] #34 kt (17.49 m/s) wind radii maximum extent in northeastern quadrant in (m) 
    rSE34ktSI=[float(x)*1852 for x in rSE34ktMile] #34 kt (17.49 m/s) wind radii maximum extent in southeastern quadrant in (m) 
    rSW34ktSI=[float(x)*1852 for x in rSW34ktMile] #34 kt (17.49 m/s) wind radii maximum extent in southwestern quadrant in (m) 
    rNW34ktSI=[float(x)*1852 for x in rNW34ktMile] #34 kt (17.49 m/s) wind radii maximum extent in northwestern quadrant in (m) 
    rNE50ktSI=[float(x)*1852 for x in rNE50ktMile] #50 kt (25.72 m/s) wind radii maximum extent in northeastern quadrant in (m) 
    rSE50ktSI=[float(x)*1852 for x in rSE50ktMile] #50 kt (25.72 m/s) wind radii maximum extent in southeastern quadrant in (m) 
    rSW50ktSI=[float(x)*1852 for x in rSW50ktMile] #50 kt (25.72 m/s) wind radii maximum extent in southwestern quadrant in (m) 
    rNW50ktSI=[float(x)*1852 for x in rNW50ktMile] #50 kt (25.72 m/s) wind radii maximum extent in northwestern quadrant in (m) 
    rNE64ktSI=[float(x)*1852 for x in rNE64ktMile] #64 kt (32.92 m/s) wind radii maximum extent in northeastern quadrant in (m) 
    rSE64ktSI=[float(x)*1852 for x in rSE64ktMile] #64 kt (32.92 m/s) wind radii maximum extent in southeastern quadrant in (m) 
    rSW64ktSI=[float(x)*1852 for x in rSW64ktMile] #64 kt (32.92 m/s) wind radii maximum extent in southwestern quadrant in (m) 
    rNW64ktSI=[float(x)*1852 for x in rNW64ktMile] #64 kt (32.92 m/s) wind radii maximum extent in northwestern quadrant in (m) 
    
    #Replacing Missing values with 'NaN'
    for i in range(0,len(MaxSustWindvelSI),1):
        if MaxSustWindvelKt[i]==-999: MaxSustWindvelSI[i]=np.nan #Maximum sustained wind velocity in (m/s) 
        if MaxSustWindvelKt[i]==-99: MaxSustWindvelSI[i]=np.nan #Maximum sustained wind velocity in (m/s) 
        if MinPressureMb[i]==-999: MinPressureSI[i]=np.nan #Minimum Pressure in (Pa)
        if rNE34ktMile[i]==-999: rNE34ktSI[i]=np.nan #34 kt (17.49 m/s) wind radii maximum extent in northeastern quadrant in (m) 
        if rSE34ktMile[i]==-999: rSE34ktSI[i]=np.nan #34 kt (17.49 m/s) wind radii maximum extent in southeastern quadrant in (m) 
        if rSW34ktMile[i]==-999: rSW34ktSI[i]=np.nan #34 kt (17.49 m/s) wind radii maximum extent in southwestern quadrant in (m) 
        if rNW34ktMile[i]==-999: rNW34ktSI[i]=np.nan #34 kt (17.49 m/s) wind radii maximum extent in northwestern quadrant in (m) 
        if rNE50ktMile[i]==-999: rNE50ktSI[i]=np.nan #50 kt (25.72 m/s) wind radii maximum extent in northeastern quadrant in (m) 
        if rSE50ktMile[i]==-999: rSE50ktSI[i]=np.nan #50 kt (25.72 m/s) wind radii maximum extent in southeastern quadrant in (m) 
        if rSW50ktMile[i]==-999: rSW50ktSI[i]=np.nan #50 kt (25.72 m/s) wind radii maximum extent in southwestern quadrant in (m) 
        if rNW50ktMile[i]==-999: rNW50ktSI[i]=np.nan #50 kt (25.72 m/s) wind radii maximum extent in northwestern quadrant in (m) 
        if rNE64ktMile[i]==-999: rNE64ktSI[i]=np.nan #64 kt (32.92 m/s) wind radii maximum extent in northeastern quadrant in (m) 
        if rSE64ktMile[i]==-999: rSE64ktSI[i]=np.nan #64 kt (32.92 m/s) wind radii maximum extent in southeastern quadrant in (m) 
        if rSW64ktMile[i]==-999: rSW64ktSI[i]=np.nan #64 kt (32.92 m/s) wind radii maximum extent in southwestern quadrant in (m) 
        if rNW64ktMile[i]==-999: rNW64ktSI[i]=np.nan #64 kt (32.92 m/s) wind radii maximum extent in northwestern quadrant in (m) 
    
    #Two examples of replaving elements in list
    #MaxSustWindvelSI=[x if x!=-999 else np.nan for x in MaxSustWindvelSI]
    #MaxSustWindvelSI=[np.nan if x==-999 else x for x in MaxSustWindvelSI]
    

    #Combining all radii values in one array
    rMile=np.column_stack((rNE34ktMile,rSE34ktMile,rSW34ktMile,rNW34ktMile,rNE50ktMile,rSE50ktMile,rSW50ktMile,rNW50ktMile,rNE64ktMile,rSE64ktMile,rSW64ktMile,rNW64ktMile))
    rMeter=np.column_stack((rNE34ktSI,rSE34ktSI,rSW34ktSI,rNW34ktSI,rNE50ktSI,rSE50ktSI,rSW50ktSI,rNW50ktSI,rNE64ktSI,rSE64ktSI,rSW64ktSI,rNW64ktSI))

    #--------------------------------------------------------------------------
    #Changing directory to working directory

    os.chdir(currentFolder)

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
    
        #Plotting
        plt.scatter(longtrack,lattrack,label='Best Track')
    
        plt.xlabel('Longitude (Degree)')
        plt.ylabel('Latitude (Degree)')
        plt.legend()
    
    #--------------------------------------------------------------------------
    #Output
    return tsdate, tstime, lattrack, longtrack, MaxSustWindvelSI, MinPressureSI, rMeter, MaxSustWindvelKt, MinPressureMb, rMile, recordid, systemstatus

    #--------------------------------------------------------------------------
