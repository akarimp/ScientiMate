function [SustWindDur] = sustainedwindduration(windvel, winddir, MaxwindvelVar, MaxwinddirVar, WindTimeInterval, window_length, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2021-01-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sustainedwindduration
=====================

.. code:: MATLAB

    [SustWindDur] = sustainedwindduration(windvel, winddir, MaxwindvelVar, MaxwinddirVar, WindTimeInterval, window_length, dispout)

Description
-----------

Calculate the sustained wind duration

Inputs
------

windvel
    Wind velocity time series data in (m/s)
winddir
    Wind direction time series data in (Degree)
MaxwindvelVar=2.5;
    | Maximum allowed wind velocity variation around a mean to allow wind to be considered sustained in (m/s)
    | Coastal Engineering Manual (2015) suggests 2.5 m/s
MaxwinddirVar=15;
    | Maximum allowed wind direction variation around a mean to allow wind to be considered sustained in (Degree)
    | Coastal Engineering Manual (2015) suggests 15 degree
WindTimeInterval=3600;
    | Time interval (time step) between two consecutive wind measurements (data points) in (second)
    | Example: for 60-min measurement interval: WindTimeInterval=3500 (s)
window_length=7;
    Number of data points in sliding window for calculating moving average values, must be odd
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

SustWindDur
    Sustained Wind Duration in (second)

Examples
--------

.. code:: MATLAB

    %Data from https://tidesandcurrents.noaa.gov for Grand Isle, LA, USA (8761724), for June 1st, 2017, reported hourly
    windvel=[3;4.7;4.9;5.3;3.3;3.4;3.3;3.8;3.1;2;1.3;1.2;1.5;3.2;2.9;3;2.9;3.7;3.7;3.1;3.4;2.6;2.5;2.5]; %24 Hour wind velocity
    winddir=[78;86;88;107;131;151;163;163;153;150;148;105;105;75;95;103;97;103;108;111;124;183;171;113]; %24 Hour wind direction
    [SustWindDur]=sustainedwindduration(windvel,winddir,2.5,15,3600,7,'yes');

References
----------

U.S. Army Corps of Engineers (2015). 
Coastal Engineering Manual. 
Engineer Manual 1110-2-1100, Washington, D.C.: U.S. Army Corps of Engineers.

Yamartino, R. J. (1984). 
A comparison of several "single-pass" estimators of the standard deviation of wind direction. 
Journal of Climate and Applied Meteorology, 23(9), 1362-1366.

Change point detection:
* https://en.wikipedia.org/wiki/Change_detection
* https://en.wikipedia.org/wiki/Step_detection
* https://www.mathworks.com/help/signal/ref/findchangepts.html
* https://www.mathworks.com/help/matlab/ref/ischange.html
* https://github.com/deepcharles/ruptures

.. License & Disclaimer
.. --------------------
..
.. Copyright (c) 2021 Arash Karimpour
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
        MaxwindvelVar=2.5; MaxwinddirVar=15; WindTimeInterval=3600; window_length=7; dispout='no';
    case 3
        MaxwinddirVar=15; WindTimeInterval=3600; window_length=7; dispout='no';
    case 4
        WindTimeInterval=3600; window_length=7; dispout='no';
    case 4
        window_length=7; dispout='no';
    case 5
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
%Calculate moving average wind

%Input size
[M,N]=size(windvel);

%Recatngular normal window
window_length_half=(window_length-1)/2;
rect_window_norm=ones(window_length,1)./window_length;

%Moving average wind velocity
windvel_moveavg = conv(windvel, rect_window_norm, 'same');
windvel_moveavg(1:window_length_half,1)=nanmean(windvel(1:window_length_half,1));
windvel_moveavg(end-window_length_half+1:end,1)=nanmean(windvel(end-window_length_half+1:end,1));

%Moving average wind velocity difference
windvel_diff = diff(windvel);
%windvel_grad = gradient(windvel,2);

%Find large variations and replace them with original data
index_windvel_diff_large=find(abs(windvel_diff)>=30);
if isempty(index_windvel_diff_large)==0
    for i=1:length(index_windvel_diff_large(:,1));
        index_to_recover=index_windvel_diff_large(i,1)-window_length_half:index_windvel_diff_large(i,1)+window_length_half;
        index_to_recover(index_to_recover<1)=1;
        index_to_recover(index_to_recover>M)=M;
        windvel_moveavg(index_to_recover)=windvel(index_to_recover);
    end
end

%Moving average wind direction
winddir_moveavg = conv(winddir, rect_window_norm, 'same');
winddir_moveavg(1:window_length_half,1)=nanmean(winddir(1:window_length_half,1));
winddir_moveavg(end-window_length_half+1:end,1)=nanmean(winddir(end-window_length_half+1:end,1));

%Moving average wind direction difference
winddir_diff = diff(winddir);
%winddir_grad = gradient(winddir,2);

%Find large variations and replace them with original data
index_winddir_diff_large=find(abs(winddir_diff)>=30);
if isempty(index_winddir_diff_large)==0
    for i=1:length(index_winddir_diff_large(:,1));
        index_to_recover=index_winddir_diff_large(i,1)-window_length_half:index_winddir_diff_large(i,1)+window_length_half;
        index_to_recover(index_to_recover<1)=1;
        index_to_recover(index_to_recover>M)=M;
        winddir_moveavg(index_to_recover)=winddir(index_to_recover);
    end
end

%--------------------------------------------------------------------------
%Calculate sustained wind duration

%Pre-assigning array
delta_windvel=zeros(M,1); %Pre-assigning array
delta_winddir=zeros(M,1); %Pre-assigning array
SustWind_Identifier=zeros(M,1); %Pre-assigning array
SustWindDur_Counter=zeros(M,1); %Pre-assigning array
SustWindDur=zeros(M,1); %Pre-assigning array

%Assigning first element
SustWind_Identifier(1,1)=1; %Sustained wind index (1: Sustained, 0: Un-sustained)
SustWindDur_Counter(1,1)=SustWind_Identifier(1,1); %Sustained wind duration index
SustWindDur(1,1)=SustWindDur_Counter(1,1)*WindTimeInterval; % Sustained Wind Duration

%Calculating the sustained wind duration
k=1; %Sustained wind index counter
for i=2:M
    
    delta_windvel(i,1)=windvel(i,1)-windvel_moveavg(i-1,1);
    delta_winddir(i,1)=winddir(i,1)-winddir_moveavg(i-1,1);
    
    %Check if winddir_moveavg is right before 360 degree and winddir is right after 360 (0) degree
    if winddir_moveavg(i-1,1)>(360-MaxwinddirVar) & winddir_moveavg(i-1,1)<=360 & winddir_moveavg(i,1)>=0 & winddir_moveavg(i,1)<MaxwinddirVar
        delta_winddir(i,1)=(winddir(i,1)+360)-winddir_moveavg(i-1,1);
    end
    
    %Check if winddir_moveavg is right after 360 (0) degree and winddir is right before 360 degree
    if winddir_moveavg(i-1,1)>=0 & winddir_moveavg(i-1,1)<MaxwinddirVar & winddir_moveavg(i,1)>(360-MaxwinddirVar) & winddir_moveavg(i,1)<=360
        delta_winddir(i,1)=(360-winddir(i,1))-winddir_moveavg(i-1,1);
    end
    
    %Calculate sustained wind duration counter, each counter means duration equal to 1 WindTimeInterval
    if abs(delta_windvel(i,1))<=MaxwindvelVar & abs(delta_winddir(i,1))<=MaxwinddirVar %Check for Sustained Wind (1: Sustained, 0: Un-sustained)
        SustWind_Identifier(i,1)=1; %Check for Sustained Wind (1: Sustained, 0: Un-sustained)
        SustWindDur_Counter(i,1)=SustWindDur_Counter(i-1,1)+SustWind_Identifier(i,1); %Sustained Wind Duration Unit
    else
        SustWind_Identifier(i,1)=0; %Check for Sustained Wind (1: Sustained, 0: Un-sustained)
        SustWindDur_Counter(i,1)=1;
        k=i-1;
    end
    
end

SustWindDur=SustWindDur_Counter.*WindTimeInterval; %Sustained Wind Duration
SustWindDur(SustWindDur_Counter==0)=1*WindTimeInterval; %Shortest Sustained Wind Duration is WindTimeInterval seceond
SustWindDur(SustWindDur<WindTimeInterval)=WindTimeInterval; %Shortest Sustained Wind Duration is WindTimeInterval seceond

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1

    Time(:,1)=[1:1:M]*WindTimeInterval;

    subplot(3,1,1)
    plot(Time./3600,windvel)
    hold on
    plot(Time./3600,windvel_moveavg)
    xlabel('Time (Hour)');
    ylabel('Velocity (m/s)')
    title('Velocity')
    legend('Wind Velcity','Moving Average')

    subplot(3,1,2)
    plot(Time./3600,winddir)
    hold on
    plot(Time./3600,winddir_moveavg)
    xlabel('Time (Hour)');
    ylabel('Direction (Degree)')
    title('Direction')
    legend('Wind Direction','Moving Average')

    subplot(3,1,3)
    plot(Time./3600,SustWindDur./3600)
    xlabel('Time (Hour)');
    ylabel('Sustained Wind Duration (Hour)')
    title('Sustained Wind Duration')

end 

%--------------------------------------------------------------------------
