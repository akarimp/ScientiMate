function [tsdate, tstime, lattrack, longtrack, MaxSustWindvelSI, MinPressureSI, rMeter, MaxSustWindvelKt, MinPressureMb, rMile,recordid, systemstatus] ...
    = readnhchurricane(filename, filelocation, hurricanename, hurricaneyear, hearderlinelen, dispout)
%{
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

readnhchurricane
================

.. code:: MATLAB

    [tsdate, tstime, lattrack, longtrack, MaxSustWindvelSI, MinPressureSI, rMeter, MaxSustWindvelKt, MinPressureMb, rMile,recordid, systemstatus] ...
    = readnhchurricane(filename, filelocation, hurricanename, hurricaneyear, hearderlinelen, dispout)

Description
-----------

Read and extracts hurricane data from National Hurricane Center (NHC) HURDAT2 file

Inputs
------

filename='hurdat2.txt';
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
filelocation=pwd;
    Location of HURDAT2 file between ' ' mark, example: 'C:\'
hurricanename='KATRINA';
    | Hurricane name between ' ' mark, example: 'KATRINA'
    | 'all' will read all data in the file
hurricaneyear=2005;
    | Year that hurricane occurs
    | Hurricane Katrina occured in 2005
    | if hurricanename='all'; then hurricaneyear is ignored
hearderlinelen=37;
    | Number of charachters in header line
    | In HURDAT2 file, hurricane header line has 37 charachters (including spaces)
    | In HURDAT2 file, hurricane data line has 120 charachters (including spaces)
dispout='no';
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
    | Minimum (central) pressure in (Pa)
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

.. code:: MATLAB

    filename='hurdat2.txt';
    filelocation=pwd;
    filelocation='C:';
    [tsdate,tstime,lattrack,longtrack,MaxSustWindvelSI,MinPressureSI,rMeter,MaxSustWindvelKt,MinPressureMb,rMile,recordid,systemstatus]...
        =readnhchurricane(filename,filelocation,'KATRINA',2005,37,'yes');

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
%}

%--------------------------------------------------------------------------
%CODE
%--------------------------------------------------------------------------
%Assign default values

switch nargin
    case 0
        filename='hurdat2.txt'; filelocation=pwd; hurricanename='KATRINA'; hurricaneyear=2005; hearderlinelen=37; dispout='no';
    case 1
        filelocation=pwd; hurricanename='KATRINA'; hurricaneyear=2005; hearderlinelen=37; dispout='no';
    case 2
        hurricanename='KATRINA'; hurricaneyear=2005; hearderlinelen=37; dispout='no';
    case 3
        hurricaneyear=2005; hearderlinelen=37; dispout='no';
    case 4
        hearderlinelen=37; dispout='no';
    case 5
        dispout='no';
end

%--------------------------------------------------------------------------
%Reading xyz file

%Reading all data
currentFolder=pwd;
cd(filelocation)

%Reading HURDAT2 data into cell array
fid=fopen(filename);
txtline=fgetl(fid); %Text line
txtlines=cell(0,1); %All text lines
while ischar(txtline)
    txtlines{end+1,1}=txtline;
    txtline=fgetl(fid);
end
fclose(fid);

%--------------------------------------------------------------------------
%Defining hearder lines and data lines for each hurricane

%Reading data of all hurricanes in the file
if strcmp(hurricanename,'all')==1

    %Reading data 
    data=cell(0,1);
    for i=1:length(txtlines)
        if length(txtlines{i})~=hearderlinelen %Checking that txtlines{i} is not header line
            data{end+1,1}=textscan(txtlines{i},'%f %f %s %s %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f','delimiter',',');
        end
    end

%Reading data for hurricane with specific name in hurricanename
else

    %Obtaining length of each line in text cell array
    for i=1:length(txtlines)
        lentxtlines(i,1)=length(txtlines{i});
    end

    %Defining index of header lines
    %Line associated to hurricane headline has 37 charachters
    %Line associated to hurricane data has 120 charachters
    IndxHurrHeader=find(lentxtlines==hearderlinelen);

    %Defining header lines contains a desired hurricane name, which might be more than one
    %For example hurricane with name Katrina occured 7 times between 1967 to 2005
    Indx1=strfind(txtlines,hurricanename);
    IndxHurrName=find(not(cellfun('isempty',Indx1)));

    %Defining lines contains a desired hurricane year
    Indx1=strfind(txtlines,num2str(hurricaneyear));
    IndxHurrYear=find(not(cellfun('isempty',Indx1)));

    %Defining a header line associated to a desired hurricane name in that specific year
    IndxHurrHeaderReq=intersect(IndxHurrYear,IndxHurrName);

    %Defining index of a desired hurricane header line in IndxHurrHeader
    IndxinIndxHurrHeader=find(IndxHurrHeader==IndxHurrHeaderReq);

    %Defining first line and last line containing data of a desired hurricane
    DataStartIndx=IndxHurrHeaderReq+1; %Index of first line containing data of a desired hurricane
    if IndxHurrHeader(IndxinIndxHurrHeader,1)<IndxHurrHeader(end,1)
        DataEndIndx=IndxHurrHeader(IndxinIndxHurrHeader+1,1)-1; %Index of last line containing data of a desired hurricane
    else %If hurricane is the last hurricane in HURDAT2 file
        DataEndIndx=length(txtlines); %Index of last line containing data of a desired hurricane
    end  

    %Reading data 
    data=cell(0,1);
    for i=DataStartIndx:DataEndIndx
        data{end+1,1}=textscan(txtlines{i},'%f %f %s %s %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f','delimiter',',');
    end

end

%Extracting data
recordid=cell(0,1);
systemstatus=cell(0,1);
latbtStr=cell(0,1);
longbtStr=cell(0,1);
for i=1:length(data)
    tsdate(i,1)=data{i}{1}; %Time series date in YYYYMMDD format
    tstime(i,1)=data{i}{2}; %Time series time in HHMM format
    recordid{end+1,1}=data{i}{3}{1}; %Record identifier
    systemstatus{end+1,1}=data{i}{4}{1}; %Status of system
    latbtStr{end+1,1}=data{i}{5}{1}; %Latitude of best track as string in degree
    longbtStr{end+1,1}=data{i}{6}{1}; %Longitude of best track as string in degree
    MaxSustWindvelKt(i,1)=data{i}{7}; %Maximum sustained wind velocity (in knots) 
    MinPressureMb(i,1)=data{i}{8}; %Minimum Pressure (in millibars)
    rNE34ktMile(i,1)=data{i}{9}; %34 kt (17.49 m/s) wind radii maximum extent in northeastern quadrant (in nautical miles) 
    rSE34ktMile(i,1)=data{i}{10}; %34 kt (17.49 m/s) wind radii maximum extent in southeastern quadrant (in nautical miles) 
    rSW34ktMile(i,1)=data{i}{11}; %34 kt (17.49 m/s) wind radii maximum extent in southwestern quadrant (in nautical miles) 
    rNW34ktMile(i,1)=data{i}{12}; %34 kt (17.49 m/s) wind radii maximum extent in northwestern quadrant (in nautical miles) 
    rNE50ktMile(i,1)=data{i}{13}; %50 kt (25.72 m/s) wind radii maximum extent in northeastern quadrant (in nautical miles) 
    rSE50ktMile(i,1)=data{i}{14}; %50 kt (25.72 m/s) wind radii maximum extent in southeastern quadrant (in nautical miles) 
    rSW50ktMile(i,1)=data{i}{15}; %50 kt (25.72 m/s) wind radii maximum extent in southwestern quadrant (in nautical miles) 
    rNW50ktMile(i,1)=data{i}{16}; %50 kt (25.72 m/s) wind radii maximum extent in northwestern quadrant (in nautical miles) 
    rNE64ktMile(i,1)=data{i}{17}; %64 kt (32.92 m/s) wind radii maximum extent in northeastern quadrant (in nautical miles) 
    rSE64ktMile(i,1)=data{i}{18}; %64 kt (32.92 m/s) wind radii maximum extent in southeastern quadrant (in nautical miles) 
    rSW64ktMile(i,1)=data{i}{19}; %64 kt (32.92 m/s) wind radii maximum extent in southwestern quadrant (in nautical miles) 
    rNW64ktMile(i,1)=data{i}{20}; %64 kt (32.92 m/s) wind radii maximum extent in northwestern quadrant (in nautical miles) 
end    

%Separating best track latitude and longitude values from N, S, E and W text
latbtNumStr=cell(0,1);
longbtNumStr=cell(0,1);
for i=1:length(data)
    latbtNumStr{end+1,1}=textscan(latbtStr{i},'%f%s');
    longbtNumStr{end+1,1}=textscan(longbtStr{i},'%f%s');
end

%Converting best track latitude and longitude from N, S, E and W format to -90 to 90 and -180 to 180
for i=1:length(latbtNumStr)

    if latbtNumStr{i}{2}{1}=='N'
        lattrack(i,1)=latbtNumStr{i}{1};
    elseif latbtNumStr{i}{2}{1}=='S'
        lattrack(i,1)=-latbtNumStr{i}{1};
    end        

    if longbtNumStr{i}{2}{1}=='E'
        longtrack(i,1)=longbtNumStr{i}{1};
    elseif longbtNumStr{i}{2}{1}=='W'
        longtrack(i,1)=-longbtNumStr{i}{1};
    end        

end    

%Converting values to SI Units
MaxSustWindvelSI=MaxSustWindvelKt.*0.514444; %Maximum sustained wind velocity in (m/s) 
MinPressureSI=MinPressureMb.*100; %Minimum Pressure in (Pa)
rNE34ktSI=rNE34ktMile.*1852; %34 kt (17.49 m/s) wind radii maximum extent in northeastern quadrant in (m) 
rSE34ktSI=rSE34ktMile.*1852; %34 kt (17.49 m/s) wind radii maximum extent in southeastern quadrant in (m) 
rSW34ktSI=rSW34ktMile.*1852; %34 kt (17.49 m/s) wind radii maximum extent in southwestern quadrant in (m) 
rNW34ktSI=rNW34ktMile.*1852; %34 kt (17.49 m/s) wind radii maximum extent in northwestern quadrant in (m) 
rNE50ktSI=rNE50ktMile.*1852; %50 kt (25.72 m/s) wind radii maximum extent in northeastern quadrant in (m) 
rSE50ktSI=rSE50ktMile.*1852; %50 kt (25.72 m/s) wind radii maximum extent in southeastern quadrant in (m) 
rSW50ktSI=rSW50ktMile.*1852; %50 kt (25.72 m/s) wind radii maximum extent in southwestern quadrant in (m) 
rNW50ktSI=rNW50ktMile.*1852; %50 kt (25.72 m/s) wind radii maximum extent in northwestern quadrant in (m) 
rNE64ktSI=rNE64ktMile.*1852; %64 kt (32.92 m/s) wind radii maximum extent in northeastern quadrant in (m) 
rSE64ktSI=rSE64ktMile.*1852; %64 kt (32.92 m/s) wind radii maximum extent in southeastern quadrant in (m) 
rSW64ktSI=rSW64ktMile.*1852; %64 kt (32.92 m/s) wind radii maximum extent in southwestern quadrant in (m) 
rNW64ktSI=rNW64ktMile.*1852; %64 kt (32.92 m/s) wind radii maximum extent in northwestern quadrant in (m) 

%Replacing Missing values with 'NaN'
MaxSustWindvelSI(MaxSustWindvelKt==-999)=NaN; %Maximum sustained wind velocity in (m/s) 
MaxSustWindvelSI(MaxSustWindvelKt==-99)=NaN; %Maximum sustained wind velocity in (m/s) 
MinPressureSI(MinPressureMb==-999)=NaN; %Minimum Pressure in (Pa)
rNE34ktSI(rNE34ktMile==-999)=NaN; %34 kt (17.49 m/s) wind radii maximum extent in northeastern quadrant in (m) 
rSE34ktSI(rSE34ktMile==-999)=NaN; %34 kt (17.49 m/s) wind radii maximum extent in southeastern quadrant in (m) 
rSW34ktSI(rSW34ktMile==-999)=NaN; %34 kt (17.49 m/s) wind radii maximum extent in southwestern quadrant in (m) 
rNW34ktSI(rNW34ktMile==-999)=NaN; %34 kt (17.49 m/s) wind radii maximum extent in northwestern quadrant in (m) 
rNE50ktSI(rNE50ktMile==-999)=NaN; %50 kt (25.72 m/s) wind radii maximum extent in northeastern quadrant in (m) 
rSE50ktSI(rSE50ktMile==-999)=NaN; %50 kt (25.72 m/s) wind radii maximum extent in southeastern quadrant in (m) 
rSW50ktSI(rSW50ktMile==-999)=NaN; %50 kt (25.72 m/s) wind radii maximum extent in southwestern quadrant in (m) 
rNW50ktSI(rNW50ktMile==-999)=NaN; %50 kt (25.72 m/s) wind radii maximum extent in northwestern quadrant in (m) 
rNE64ktSI(rNE64ktMile==-999)=NaN; %64 kt (32.92 m/s) wind radii maximum extent in northeastern quadrant in (m) 
rSE64ktSI(rSE64ktMile==-999)=NaN; %64 kt (32.92 m/s) wind radii maximum extent in southeastern quadrant in (m) 
rSW64ktSI(rSW64ktMile==-999)=NaN; %64 kt (32.92 m/s) wind radii maximum extent in southwestern quadrant in (m) 
rNW64ktSI(rNW64ktMile==-999)=NaN; %64 kt (32.92 m/s) wind radii maximum extent in northwestern quadrant in (m) 

%Combining all radii values in one array
rMile=[rNE34ktMile,rSE34ktMile,rSW34ktMile,rNW34ktMile,rNE50ktMile,rSE50ktMile,rSW50ktMile,rNW50ktMile,rNE64ktMile,rSE64ktMile,rSW64ktMile,rNW64ktMile];
rMeter=[rNE34ktSI,rSE34ktSI,rSW34ktSI,rNW34ktSI,rNE50ktSI,rSE50ktSI,rSW50ktSI,rNW50ktSI,rNE64ktSI,rSE64ktSI,rSW64ktSI,rNW64ktSI];

%--------------------------------------------------------------------------
%Changing directory to working directory

cd(currentFolder)

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1

    %Plotting
    scatter(longtrack,lattrack)

    xlabel('Longitude (Degree)')
    ylabel('Latitude (Degree)')
    legend('Best Track')

end

%--------------------------------------------------------------------------
