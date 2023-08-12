function [x, y, z] = readxyzfile(xyzfilename, xyzfilelocation, zscale, domain, xmin, xmax, ymin, ymax, savedata, outfilename, outfilelocation)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-11-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

readxyzfile
===========

.. code:: MATLAB

    [x, y, z] = readxyzfile(xyzfilename, xyzfilelocation, zscale, domain, xmin, xmax, ymin, ymax, savedata, outfilename, outfilelocation)

Description
-----------

| Read and extract x (longitude), y (latitude) and z (elevation) data from ASCII gridded (tabular) xyz file
| Use readdatafile function for more options

Inputs
------

xyzfilename
    | Name of xyz file between ' ' mark, example: 'xyzfile.xyz'
    | xyz file should be in form of 3 coloumn format
xyzfilelocation=pwd;
    Location of xyz file between ' ' mark, example: 'C:\'
zscale=1;
    Scale z (elevation) data by factor of zscale
domain='all';
    | Define a domain to be extracted from data
    | 'all': all xyz data in input file are extracted 
    | 'domain': only data within a defined domain are extracted
xmin=-180;
    Minimum x (longitude) of domain to be extracted
xmax=180;
    Maximum x (longitude) of domain to be extracted
ymin=-90;
    Minimum y (latitude) of domain to be extracted
ymax=90;
    Maximum y (latitude) of domain to be extracted
savedata='no';
    | Define if save xyz data in a file or not in outfilelocation folder
    | 'no': does not save
    | 'yes': save xyz data as csv file
outfilename='xyzdata.csv';
    | Name of output file between ' ' mark, example: 'xyzdata.csv'
    | outfilename should have '.csv' extension
outfilelocation=pwd;
    Location of output file between ' ' mark, example: 'C:\'

Outputs
-------

x
    x (longitude) data extracted from xyz file
y
    y (latitude) data extracted from xyz file
z
    z (elevation) data extracted from xyz file

Examples
--------

.. code:: MATLAB

    xyzfilename='xyzfile.xyz'; %e.g. xyzfilename='PersianGulf_ETOPO1.xyz'
    xyzfilelocation='C:\'; %e.g. xyzfilelocation='C:\datafolder'
    [x,y,z]=readxyzfile(xyzfilename,xyzfilelocation);

    xyzfilename='xyzfile.xyz'; %e.g. xyzfilename='PersianGulf_ETOPO1.xyz'
    xyzfilelocation='C:\'; %e.g. xyzfilelocation='C:\datafolder'
    [x,y,z]=readxyzfile(xyzfilename,xyzfilelocation,1,'all',-180,180,-90,90,'no');

References
----------

Geospatial data

* https://www.mathworks.com/help/map/finding-geospatial-data.html
* https://maps.ngdc.noaa.gov/viewers/wcs-client/
* https://www.ngdc.noaa.gov/mgg/global/global.html
* https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/
* https://www.ngdc.noaa.gov/mgg/image/2minrelief.html
* https://www.ngdc.noaa.gov/mgg/coastal/crm.html
* https://viewer.nationalmap.gov/launch/
* https://earthexplorer.usgs.gov
* http://www.shadedrelief.com/cleantopo2/index.html

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
        xyzfilelocation=pwd; zscale=1; domain='all'; xmin=-180; xmax=180; ymin=-90; ymax=90; savedata='no'; outfilename='xyzdata.csv'; outfilelocation=pwd;
    case 2
        zscale=1; domain='all'; xmin=-180; xmax=180; ymin=-90; ymax=90; savedata='no'; outfilename='xyzdata.csv'; outfilelocation=pwd;
    case 3
        domain='all'; xmin=-180; xmax=180; ymin=-90; ymax=90; savedata='no'; outfilename='xyzdata.csv'; outfilelocation=pwd;
    case 4
        xmin=-180; xmax=180; ymin=-90; ymax=90; savedata='no'; outfilename='xyzdata.csv'; outfilelocation=pwd;
    case 5
        xmax=180; ymin=-90; ymax=90; savedata='no'; outfilename='xyzdata.csv'; outfilelocation=pwd;
    case 6
        ymin=-90; ymax=90; savedata='no'; outfilename='xyzdata.csv'; outfilelocation=pwd;
    case 7
        ymax=90; savedata='no'; outfilename='xyzdata.csv'; outfilelocation=pwd;
    case 8
        savedata='no'; outfilename='xyzdata.csv'; outfilelocation=pwd;
    case 9
        outfilename='xyzdata.csv'; outfilelocation=pwd;
    case 10
        outfilelocation=pwd;
end

%--------------------------------------------------------------------------
%Reading xyz file

%Reading all data
currentFolder=pwd;
cd(xyzfilelocation)
xyzdata=importdata(xyzfilename);

%Reading all xyz data in input file
if strcmp(domain,'all')==1

    x=xyzdata(:,1); %x (longitude) data in (Degree)
    y=xyzdata(:,2); %y (latitude) data in (Degree)
    z=zscale.*xyzdata(:,3); %z (elevation) data in (m)

    xyzdata1(:,1)=x;  %x (longitude) data in (Degree)
    xyzdata1(:,2)=y;  %y (latitude) data in (Degree)
    xyzdata1(:,3)=z; %z (elevation) data in (m)

%Reading data within a defined domain
elseif strcmp(domain,'domain')==1
    
    %All data
    x1=xyzdata(:,1); %x (longitude) data in (Degree)
    y1=xyzdata(:,2); %y (latitude) data in (Degree)
    z1=zscale.*xyzdata(:,3); %z (elevation) data in (m)
    
    %Finding a location data that are within [xmin:xmax] range
    xIndx=find(x1>=xmin & x1<=xmax);
    
    %Retaining data that are within [xmin:xmax] range
    x2=x1(xIndx,1);   %xdata
    y2=y1(xIndx,1);   %ydata
    z2=z1(xIndx,1); %zdata
    
    %Finding a location data that are within [ymin:ymax] range
    yIndx=find(y2>=ymin & y2<=ymax);
    
    %Retaining data that are within [ymin:ymax] range
    x=x2(yIndx,1);   %x data
    y=y2(yIndx,1);   %y data
    z=z2(yIndx,1); %z data
    
    xyzdata1(:,1)=x;  %x (longitude) data in (Degree)
    xyzdata1(:,2)=y;  %y (latitude) data in (Degree)
    xyzdata1(:,3)=z; %z (elevation) data in (m)
    
end

%--------------------------------------------------------------------------
%Saving xyz data

if strcmp(savedata,'yes')==1
    cd(outfilelocation)
    %csvwrite('xyzdata.csv',xyzdata1)
    %dlmwrite('xyzdata.txt',xyzdata1)
    dlmwrite(outfilename,xyzdata1,'delimiter',',')
end

%--------------------------------------------------------------------------
%Changing directory to working directory

cd(currentFolder)

%--------------------------------------------------------------------------
