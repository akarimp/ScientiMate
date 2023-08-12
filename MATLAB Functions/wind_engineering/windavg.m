function [windvelavg, winddiravg] = windavg(windvel, winddir, NPointsAvg, NPointsInterval, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-07-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

windavg
=======

.. code:: MATLAB

    [windvelavg, winddiravg] = windavg(windvel, winddir, NPointsAvg, NPointsInterval, dispout)

Description
-----------

Average wind velocity and wind direction

Inputs
------

windvel
    Wind velocity time series data
winddir
    Wind direction time series data in (degree)
NPointsAvg=length(windvel(:,1));
    Number of data points from start of each section (interval) to be averaged
NPointsInterval=length(windvel(:,1));
    Number of points that each section (interval) has
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

windvelavg
    Averaged wind velocity data
winddiravg
    Averaged wind direction data in (degree)

Examples
--------

.. code:: MATLAB

    windvel(:,1)=10.*rand(10,1);
    winddir(:,1)=45.*rand(10,1);
    [windvelavg,winddiravg]=windavg(windvel,winddir);

    windvel(:,1)=10.*rand(5*60,1); %One data point every minute for 5 hours
    winddir(:,1)=225.*rand(5*60,1); %One data point every minute for 5 hours
    [windvelavg,winddiravg]=windavg(windvel,winddir,10,60,'yes');

References
----------

Yamartino, R. J. (1984). 
A comparison of several "single-pass" estimators of the standard deviation of wind direction. 
Journal of Climate and Applied Meteorology, 23(9), 1362-1366.

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
    case 2
        NPointsAvg=length(windvel(:,1)); NPointsInterval=length(windvel(:,1)); dispout='no';
    case 3
        NPointsInterval=length(windvel(:,1)); dispout='no';
    case 4
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(windvel)==1
    windvel=windvel';
end

if isrow(winddir)==1
    winddir=winddir';
end

%--------------------------------------------------------------------------
%Calculating number of bursts in the input file

burst=fix(length(windvel(:,1))/NPointsInterval); %Number of burst in input file
sample=NPointsInterval; %Number of data points in 1 burst
NPointsAvg(NPointsAvg>NPointsInterval)=NPointsInterval; %Checking number of points to be averaged

%--------------------------------------------------------------------------
%Wind velocity averaging 

windvelavg=zeros(burst,1); #Pre-assigning array
for i=1:burst
    j1=(i-1)*sample+1;
    j2=j1+(NPointsAvg-1);
    windvelavg(i,1)=nanmean(windvel(j1:j2,1));
end

%--------------------------------------------------------------------------
%Wind direction averaging using Yamartino (1984) method

%Notes:
%Use arctan2 instead of arctan to get the averaged wind direction between -180 and 180
%Add 360 to all numbers to have them all positive
%Use mod(360) to take care of the ones larger than 360

winddiravg=zeros(burst,1); %Pre-assigning array
DirSimpleAvg=zeros(burst,1); %Pre-assigning array
for i=1:burst
    j1=(i-1)*sample+1;
    j2=j1+(NPointsAvg-1);
    DirSimpleAvg(i,1)=nanmean(winddir(j1:j2,1));
    x=cos(deg2rad(winddir(j1:j2,1)));
    y=sin(deg2rad(winddir(j1:j2,1)));
    xmean=nanmean(x(:,1));
    ymean=nanmean(y(:,1));
    Dir=rad2deg(atan2(ymean,xmean));
    winddiravg(i,1)=mod((Dir+360),360);
end

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1

    dt1(:,1)=[0:1:length(windvel(:,1))-1]; %Delat_t for all data points 
    dt2(:,1)=[0:NPointsInterval:(burst-1)*NPointsInterval]; %Delat_t for averaged data points 

    subplot(2,1,1)
    plot(dt1,windvel)
    hold on
    plot(dt2,windvelavg)
    xlabel('dt (s)');
    ylabel('Velocity')
    title('Velocity')
    legend('Input Data','Avg Data')

    subplot(2,1,2)
    plot(dt1,winddir)
    hold on
    plot(dt2,winddiravg)
    xlabel('dt (s)');
    ylabel('Direction (Degree)')
    title('Direction')
    legend('Input Data','Avg Data')

end 

%--------------------------------------------------------------------------
 