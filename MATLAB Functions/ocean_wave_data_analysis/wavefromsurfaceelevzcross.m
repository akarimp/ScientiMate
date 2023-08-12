function [Hs, Ts, Hz, Tz, Hrms, H, T] = wavefromsurfaceelevzcross(Eta, fs, dispout)
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

wavefromsurfaceelevzcross
=========================

.. code:: MATLAB

    [Hs, Ts, Hz, Tz, Hrms, H, T] = wavefromsurfaceelevzcross(Eta, fs, dispout)

Description
-----------

Calculate wave properties from water surface elevation by using an upward zero crossing method

Inputs
------

Eta
    Water surface elevation time series data in (m)
fs
    Sampling frequency that data collected at in (Hz)
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

Hs
    Significant Wave Height (m)
Ts
    Significant Wave Period (second)
Hz
    Zero Crossing Mean Wave Height (m)
Tz
    Zero Crossing Mean Wave Period (second)
Hrms
    Root Mean Square Wave Height (m)
H
    Wave Height Data Series array (m)
T
    Wave Period Data Series array (second)

Examples
--------

.. code:: MATLAB

    fs=2; %Sampling frequency
    duration=1024; %Duration of the data
    N=fs*duration; %Total number of points
    df=fs/N; %Frequency difference 
    dt=1/fs; %Time difference, dt=1/fs
    t(:,1)=linspace(0,duration-dt,N); %Time
    Eta(:,1)=detrend(0.5.*cos(2*pi*0.2*t)+(-0.1+(0.1-(-0.1))).*rand(N,1));
    [Hs,Ts,Hz,Tz,Hrms,H,T]=wavefromsurfaceelevzcross(Eta,fs,'yes');

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
    case 2
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(Eta)==1
    Eta=Eta';
end

%--------------------------------------------------------------------------

%Deterending input data
EtaDetrended=detrend(Eta,'linear');

%--------------------------------------------------------------------------

%Generating time vector
N=length(Eta(:,1)); %Length of a time series, total number of points in input file, N=fs*duration
duration=N/fs; %Duration time that data are collected (second), duration=N/fs
dt=1/fs; %Time difference between consecutive samples, dt=1/fs=duration/N
%t(:,1)=linspace(0,duration-dt,N); %Time from 0 to T-dt, equally spaced at dt
t(:,1)=[0:dt:duration-dt]; %Time from 0 to T-dt, equally spaced at dt

%--------------------------------------------------------------------------
%Detecting first and last wave (first and last complete crest-trough cycle)

%Detecting the start point of the first wave (first complete crest-trough cycle)
if EtaDetrended(1,1)==0 & EtaDetrended(2,1)>0 
    len3=1;
end

if (EtaDetrended(1,1)==0 & EtaDetrended(2,1)<0) | (EtaDetrended(1,1)~=0)
    for i=1:N-1
        if EtaDetrended(i,1)<0 & EtaDetrended(i+1,1)>0 
            len3=i;
            break
        end
    end
end

%Detecting the end point of last wave (last complete crest-trough cycle)
if EtaDetrended(end,1)==0 & EtaDetrended(end-1,1)<0 
    len4=N;
end

if (EtaDetrended(end,1)==0 & EtaDetrended(end-1,1)>0) | (EtaDetrended(end,1)~=0)
    for i=N-1:-1:1
        if EtaDetrended(i,1)<0 & EtaDetrended(i+1,1)>0 
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
        if EtaDetrended(i,1)<0 & EtaDetrended(i+1,1)>0
            m=m+1;
            xupcross(m,1)=t(i,1)-(t(i+1,1)-t(i,1))/(EtaDetrended(i+1,1)-EtaDetrended(i,1))*EtaDetrended(i,1);
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
        if EtaDetrended(i,1)>0 & EtaDetrended(i+1,1)<0
            n=n+1;
            xdowncross(n,1)=t(i,1)-(t(i+1,1)-t(i,1))/(EtaDetrended(i+1,1)-EtaDetrended(i,1))*EtaDetrended(i,1);
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
        
        if EtaDetrended(i,1)>EtaDetrended(i-1,1) & EtaDetrended(i,1)>EtaDetrended(i+1,1)
            
            %Check if after last crest, a trough is detected or not 
            %m==n means it is detected and the new y(i,1) is the crest of the next wave
            if EtaDetrended(i,1)>0
                
                if m==n
                    
                    m=m+1;
                    xmax(m,1)=t(i,1);
                    ymax(m,1)=EtaDetrended(i,1);
                    xmaxIndx(m,1)=i;
                    
                elseif m~=0
                        
                    %Replacingthe old crest location with new one if the new one is larger
                    if EtaDetrended(i,1)>ymax(m,1) & m~=0
                        
                        xmax(m,1)=t(i,1);
                        ymax(m,1)=EtaDetrended(i,1);
                        xmaxIndx(m,1)=i;
                            
                    end
                    
                end
                
            end
            
        end
        
    end
    
    %Detecting trough
    if i>1 & i<N
        
        if EtaDetrended(i,1)<EtaDetrended(i-1,1) & EtaDetrended(i,1)<EtaDetrended(i+1,1)
            
            if EtaDetrended(i,1)<0
                
                %Check if after last trough, a crest is detected or not 
                %n==m-1 means it is detected and the new y(i,1) is the trough of next wave
                if n==m-1
                    
                    n=n+1;
                    xmin(n,1)=t(i,1);
                    ymin(n,1)=EtaDetrended(i,1);
                    xminIndx(n,1)=i;
                    
                elseif n~=0
                        
                    %Replacingthe old crest location with new one if the new one is smaller
                    if EtaDetrended(i,1)<ymin(n,1)
                        
                        xmin(n,1)=t(i,1);
                        ymin(n,1)=EtaDetrended(i,1);
                        xminIndx(n,1)=i;
                            
                    end
                    
                end
                
            end
            
        end
        
    end
    
    
end

%Eta=EtaDetrended; %Water surface elevation time series

%--------------------------------------------------------------------------
%Calculating wave height and surface elevation of the wave crest and trough

len1=length(xmax(:,1));
len2=length(xmin(:,1));

if len1~=len2
    error('Numbers of crests and troughs are not the same')
end

for i=1:len1
    
    H(i,1)=ymax(i,1)-ymin(i,1); %Wave height
    xmean(i,1)=(xmax(i,1)+xmin(i,1))/2;
    Etac(i,1)=ymax(i,1); %Water surface elevation of the wave crest
    Etat(i,1)=ymin(i,1); %Water surface elevation of the wave trough
    
end
        
%--------------------------------------------------------------------------
%Calculating wave period 

for i=1:length(H(:,1))
    for j=1:length(xupcross(:,1))-1
        
        if xupcross(j,1)<xmean(i,1) & xupcross(j+1,1)>xmean(i,1)
            T(i,1)=xupcross(j+1,1)-xupcross(j,1);
        end
        
    end
end

%--------------------------------------------------------------------------
%Calculating wave properties

numberofwaves=length(H(:,1));
Etarms=sqrt(sum(EtaDetrended.^2)/N);

[Hsort,HsortIndex]=sort(H,'descend');
HTop3rdIndx=round(1/3*length(Hsort(:,1))); %Index where top 3rd wave height starts
Hs=mean(Hsort(1:HTop3rdIndx,1)); %Zero-crossing significant wave height

Tsort=T(HsortIndex);
TTop3rdIndx=round(1/3*length(Tsort(:,1))); %Index where top 3rd wave period starts
Ts=mean(Tsort(1:TTop3rdIndx,1)); %Zero-crossing significant wave period

Hrms=sqrt(sum(H.^2)/length(H(:,1))); %Root Mean Square Wave Height
Hz=mean(H(:,1)); %Zero-crossing mean wave height
Tz=mean(T(:,1)); %Zero-crossing mean wave period

%Hs=sqrt(2)*sqrt(sum(H.^2)/length(H(:,1))); %Significant wave height, Hs=(2^0.5)*Hrms
%Hz=sum(H)/numberofwaves; %Zero-crossing mean wave height
%Tz=duration/numberofwaves; %Zero-crossing mean wave period

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1
    
    val=[Hs Ts Hz Tz Hrms];
    name={'Hs','Ts','Hz','Tz','Hrms'};
    for i=1:length(val)
        fprintf('%14s   %g\n',name{i},val(i));
    end
    
    %Plotting data time series
    plot(t,EtaDetrended)
    hold on
    
    scatter(xmax,ymax,'*')
    scatter(xmin,ymin,'o')
    
    scatter(xupcross,yupcross,'^')
    scatter(xdowncross,ydowncross,'v')

    for i=1:length(H(:,1))
        
        plot([xmean(i,1),xmean(i,1)],[(ymax(i,1)+ymin(i,1))/2-H(i,1)/2,(ymax(i,1)+ymin(i,1))/2+H(i,1)/2],':m')
        
    end
    
    xlabel('Time (s)')
    ylabel('Water Surface Elevatio (m)')
    legend('Time Series','Crest','Trough','Up Crossing','Down Crossing','Wave Height')
    
end

%--------------------------------------------------------------------------
