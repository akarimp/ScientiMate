function [UpCrossIndx, UpCrossTime, UpCrossValue, CrestIndx, CrestTime, CrestValue, TroughIndx, TroughTime, TroughValue, t] = levelcrossing(x, CrossingLevel, fs, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-04-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

levelcrossing
=============

.. code:: MATLAB

    [UpCrossIndx, UpCrossTime, UpCrossValue, CrestIndx, CrestTime, CrestValue, TroughIndx, TroughTime, TroughValue, t] = levelcrossing(x, CrossingLevel, fs, dispout)

Description
-----------

Calculate crossing point for a given level by using an upward zero crossing method

Inputs
------

x
    Input time series (oscillatory) data
CrossingLevel=mean(x);
    | Level that crossing reported for
    | If x is deterended wave (mean=0), then set CrossingLevel=0
fs=1;
    | Sampling frequency that data collected at in (Hz)
    | If fs is not given, then default fs is fs=1;
    | If fs=1, then index of data points represents time as well
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

UpCrossIndx
    Index of up-crossing location
UpCrossTime
    Time of up-crossing location (s)
UpCrossValue
    Value of up-crossing
CrestIndx
    Index of crest location
CrestTime
    Time of wave crest location (s)
CrestValue
    Value of wave crest
TroughIndx
    Index of trough location
TroughTime
    Time of wave trough location (s)
TroughValue
    Value of wave trough
t
    Time (s)

Examples
--------

.. code:: MATLAB

    fs=2; %Sampling frequency
    duration=1024; %Duration of the data
    N=fs*duration; %Total number of points
    df=fs/N; %Frequency difference 
    dt=1/fs; %Time difference, dt=1/fs
    t(:,1)=linspace(0,duration-dt,N); %Time
    x(:,1)=detrend(0.5.*cos(2*pi*0.2*t)+(-0.1+(0.1-(-0.1))).*rand(N,1));
    [UpCrossIndx,UpCrossTime,UpCrossValue,CrestIndx,CrestTime,CrestValue,TroughIndx,TroughTime,TroughValue,t]=levelcrossing(x,mean(x),fs,'yes');

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
        CrossingLevel=mean(x); fs=1; dispout='no';
    case 2
        fs=1; dispout='no';
    case 3
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(x)==1
    x=x';
end

%--------------------------------------------------------------------------
%Checking if a crossing level has any intersection with input data

if (CrossingLevel>=max(x(:,1))) | (CrossingLevel<=min(x(:,1)))
    error('No intersection between input data and given level, try another CrossingLevel value')
end

%--------------------------------------------------------------------------

%Deterending input data
xDetrended=x-CrossingLevel;

%--------------------------------------------------------------------------

%Generating time vector
N=length(x(:,1)); %Length of a time series, total number of points in input file, N=fs*duration
duration=N/fs; %Duration time that data are collected (second), duration=N/fs
dt=1/fs; %Time difference between consecutive samples, dt=1/fs=duration/N
%t(:,1)=linspace(0,duration-dt,N); %Time from 0 to T-dt, equally spaced at dt
t(:,1)=[0:dt:duration-dt]; %Time from 0 to T-dt, equally spaced at dt

%--------------------------------------------------------------------------
%Detecting first and last wave (first and last complete crest-trough cycle)

%Detecting the start point of the first wave (first complete crest-trough cycle)
if xDetrended(1,1)==0 & xDetrended(2,1)>0 
    len3=1;
end

if (xDetrended(1,1)==0 & xDetrended(2,1)<0) | (xDetrended(1,1)~=0)
    for i=1:N-1
        if xDetrended(i,1)<0 & xDetrended(i+1,1)>0 
            len3=i;
            break
        end
    end
end

%Detecting the end point of last wave (last complete crest-trough cycle)
if xDetrended(end,1)==0 & xDetrended(end-1,1)<0 
    len4=N;
end

if (xDetrended(end,1)==0 & xDetrended(end-1,1)>0) | (xDetrended(end,1)~=0)
    for i=N-1:-1:1
        if xDetrended(i,1)<0 & xDetrended(i+1,1)>0 
            len4=i;
            break
        end
    end
end

%--------------------------------------------------------------------------
%Detecting zero crossing points

m=0;
n=0;

for i=len3:len4 
    
    %Detecting up-crossing zero-crossing points
    if i==1
        m=m+1;
        xupcross(m,1)=dt;
        yupcross(m,1)=0;
        xupcrossIndx(m,1)=i;
    end
    
    if i>1 & i<N
        if xDetrended(i,1)<0 & xDetrended(i+1,1)>0
            m=m+1;
            xupcross(m,1)=t(i,1)-(t(i+1,1)-t(i,1))/(xDetrended(i+1,1)-xDetrended(i,1))*xDetrended(i,1);
            yupcross(m,1)=0;
            xupcrossIndx(m,1)=i;
        end
    end
    
    if i==N
        m=m+1;
        xupcross(m,1)=t(end,1);
        yupcross(m,1)=0;
        xupcrossIndx(m,1)=i;
    end
    
    %Detecting down-crossing zero-crossing points
    if i>1 & i<N
        if xDetrended(i,1)>0 & xDetrended(i+1,1)<0
            n=n+1;
            xdowncross(n,1)=t(i,1)-(t(i+1,1)-t(i,1))/(xDetrended(i+1,1)-xDetrended(i,1))*xDetrended(i,1);
            ydowncross(n,1)=0;
            xdowncrossIndx(n,1)=i;
            if xdowncrossIndx(n,1)>=N
                xdowncrossIndx(n,1)=N-1;
            end
        end
    end
    
end

%--------------------------------------------------------------------------
%Detecting crest and trough

m=0;
n=0;

for i=xupcrossIndx(1,1):xupcrossIndx(end,1)
    
    %Detecting crest
    if i>1 & i<N
        
        if xDetrended(i,1)>xDetrended(i-1,1) & xDetrended(i,1)>xDetrended(i+1,1)
            
            %Check if after last crest, a trough is detected or not 
            %m==n means it is detected and the new y(i,1) is the crest of the next wave
            if xDetrended(i,1)>0
                
                if m==n
                    
                    m=m+1;
                    xmax(m,1)=t(i,1);
                    ymax(m,1)=xDetrended(i,1);
                    xmaxIndx(m,1)=i;
                    
                elseif m~=0
                        
                    %Replacingthe old crest location with new one if the new one is larger
                    if xDetrended(i,1)>ymax(m,1) & m~=0
                        
                        xmax(m,1)=t(i,1);
                        ymax(m,1)=xDetrended(i,1);
                        xmaxIndx(m,1)=i;
                            
                    end
                    
                end
                
            end
            
        end
        
    end
    
    %Detecting trough
    if i>1 & i<N
        
        if xDetrended(i,1)<xDetrended(i-1,1) & xDetrended(i,1)<xDetrended(i+1,1)
            
            if xDetrended(i,1)<0
                
                %Check if after last trough, a crest is detected or not 
                %n==m-1 means it is detected and the new y(i,1) is the trough of next wave
                if n==m-1
                    
                    n=n+1;
                    xmin(n,1)=t(i,1);
                    ymin(n,1)=xDetrended(i,1);
                    xminIndx(n,1)=i;
                    
                elseif n~=0
                        
                    %Replacingthe old crest location with new one if the new one is smaller
                    if xDetrended(i,1)<ymin(n,1)
                        
                        xmin(n,1)=t(i,1);
                        ymin(n,1)=xDetrended(i,1);
                        xminIndx(n,1)=i;
                            
                    end
                    
                end
                
            end
            
        end
        
    end
    
    
end

%x=xDetrended; %Water surface elevation time series

%--------------------------------------------------------------------------
%Assigning output

UpCrossIndx=xupcrossIndx; %Up crossing Index
UpCrossTime=xupcross; %Up crossing time
UpCrossValue=yupcross+CrossingLevel; %Up crossing value

CrestIndx=xmaxIndx; %Crest index
CrestTime=xmax; %Crest time
CrestValue=ymax+CrossingLevel; %Crest value

TroughIndx=xminIndx; %Trough index
TroughTime=xmin; %Trough time
TroughValue=ymin+CrossingLevel; %Trough value

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1
    
    %Plotting data time series
    plot(t,x)
    hold on
    plot([t(1,1),t(end,1)],[CrossingLevel,CrossingLevel])
    
    scatter(xmax,ymax+CrossingLevel,'*')
    scatter(xmin,ymin+CrossingLevel,'o')
    
    scatter(xupcross,yupcross+CrossingLevel,'^')
    scatter(xdowncross,ydowncross+CrossingLevel,'v')

    xlabel('Time (s)')
    ylabel('Wave')
    legend('Time Series','Crossing Level','Crest','Trough','Up Crossing','Down Crossing')
    
end

%--------------------------------------------------------------------------
