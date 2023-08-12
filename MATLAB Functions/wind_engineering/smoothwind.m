function [windvel_smooth, winddir_smooth] = smoothwind(windvel, winddir, window_length, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2021-02-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

smoothwind
==========

.. code:: MATLAB

    [windvel_smooth, winddir_smooth] = smoothwind(windvel, winddir, window_length, dispout)

Description
-----------

Smooth wind data using moving average window

Inputs
------

windvel=[];
    | Wind velocity time series data
    | Leave wind velocity empty if only winddir is available
winddir=[];
    | Wind direction time series data
    | Leave wind direction empty if only windvel is available
window_length=5;
    Number of data points in sliding window for calculating moving average values, must be odd
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

windvel_smooth
    Smoothed wind velocity
winddir_smooth
    | Smoothed wind direction
    | Wind directions with large variations are not smoothed

Examples
--------

.. code:: MATLAB

    %Data from https://tidesandcurrents.noaa.gov for Grand Isle, LA, USA (8761724), for June 1st, 2017, reported hourly
    windvel=[3;4.7;4.9;5.3;3.3;3.4;3.3;3.8;3.1;2;1.3;1.2;1.5;3.2;2.9;3;2.9;3.7;3.7;3.1;3.4;2.6;2.5;2.5]; %24 Hour wind velocity
    winddir=[78;86;88;107;131;151;163;163;153;150;148;105;105;75;95;103;97;103;108;111;124;183;171;113]; %24 Hour wind direction
    [windvel_smooth,winddir_smooth]=smoothwind(windvel,winddir,5,'yes');

References
----------

U.S. Army Corps of Engineers (2015). 
Coastal Engineering Manual. 
Engineer Manual 1110-2-1100, Washington, D.C.: U.S. Army Corps of Engineers.

Yamartino, R. J. (1984). 
A comparison of several "single-pass" estimators of the standard deviation of wind direction. 
Journal of Climate and Applied Meteorology, 23(9), 1362-1366.

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
    case 0
        windvel=[]; winddir=[]; window_length=5; dispout='no';
    case 1
        winddir=[]; window_length=5; dispout='no';
    case 2
        window_length=5; dispout='no';
    case 3
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(windvel)==1
    windvel=windvel.';
end

if isrow(winddir)==1
    winddir=winddir.';
end

%--------------------------------------------------------------------------
%Define rectangular moving window

%Input size
if isempty(windvel)==0
    [M,N]=size(windvel);
elseif isempty(winddir)==0
    [M,N]=size(winddir);
end

%Recatngular normal window
window_length_half=(window_length-1)/2;
rect_window_norm=ones(window_length,1)./window_length;

%--------------------------------------------------------------------------
%Calculate moving average velocity

%Moving average wind velocity
if isempty(windvel)==0

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

else
    windvel_moveavg=[];

end

%--------------------------------------------------------------------------
%Calculate moving average direction

%Moving average wind direction
if isempty(winddir)==0

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

else
    winddir_moveavg=[];

end

%--------------------------------------------------------------------------
%Calculate sustained wind duration

windvel_smooth=windvel_moveavg; %Smoothed wind velocity
winddir_smooth=winddir_moveavg; %Smoothed wind direction

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1

    subplot(2,1,1)
    if isempty(windvel)==0
        plot(windvel)
        hold on
        plot(windvel_smooth)
        xlabel('Time')
        ylabel('Velocity')
        title('Velocity')
        legend('Wind Velcity','Smoothed Wind Velcity')
    end

    subplot(2,1,2)
    if isempty(winddir)==0
        plot(winddir)
        hold on
        plot(winddir_smooth)
        xlabel('Time')
        ylabel('Direction')
        title('Direction')
        legend('Wind Direction','Smoothed Wind Direction')
    end

end 

%--------------------------------------------------------------------------
