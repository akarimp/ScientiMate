function plot2dtimeseries(x, starttime, endtime, plottype, cmapcolors, sizestyle)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2019-02-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

plot2dtimeseries
================

.. code:: MATLAB

    plot2dtimeseries(x, starttime, endtime, plottype, cmapcolors, sizestyle)

Description
-----------

Plot x data in 2-d timeseries

Inputs
------

x
    x data as a 2-d array of (M,N)
starttime='2000-10-20 00:00:00';
    | Start date and time for time series generation
    | Format: 'yyyy-mm-dd HH:MM:SS'
endtime='2100-10-20 23:00:00';
    | End date and time for time series generation
    | Format: 'yyyy-mm-dd HH:MM:SS'
plottype='bar';
    | Plot type
    | 'line': line plot
    | 'line_grid': line plot with grid lines
    | 'line_subplot': line plot with subplots
    | 'line_subplot_grid': line plot with subplots and grid lines
    | 'bar': bar plot
    | 'bar_grid': bar plot with grid lines
    | 'barh': horizontal bar plot
    | 'barh_grid': horizontal bar plot with grid lines
cmapcolors='seq';
    | Colormap style
    | 'blue': blue colormap
    | 'red': red colormap
    | 'green': green colormap
    | 'yellow': yellow colormap
    | 'purple': purple colormap
    | 'brown': brown colormap
    | 'gray': gray colormap
    | 'blue_red': blue-red colormap
    | 'red_blue': red-blue colormap
    | 'blue_green': blue-green colormap
    | 'green_blue': green-blue colormap
    | 'green_yellow': green-yellow colormap
    | 'yellow_green': yellow-green colormap
    | 'red_yellow': red-yellow colormap
    | 'yellow_red': yellow-red colormap
    | 'cyclic': cyclic/oscillation colormap 
    | 'seq': sequential colormap
    | User-defined colors may be used to generate colormap
    | User-defined colors should be defined as (M,3) array in RGB color format
    | At least two colors should be defined, i.e. M should be equal or larger than 2
    | User-defined colors values should be between 0 and 255
    | Any available colormap name such as 'cool', 'winter', etc also can be used
sizestyle='medium';
    | Plot drawing size style
    | 'small': small plot size
    | 'medium': medium plot size
    | 'large': large plot size

Outputs
-------


Examples
--------

.. code:: MATLAB

    x=rand(366,1);
    starttime='2020-01-01 00:00:00';
    endtime='2020-12-31 00:00:00';
    plot2dtimeseries(x,starttime,endtime,'line','purple','medium')

    x=rand(31,3);
    starttime='2020-01-01 00:00:00';
    endtime='2020-01-31 00:00:00';
    plot2dtimeseries(x,starttime,endtime,'line_subplot','seq','medium')

References
----------

Colormap

* https://matplotlib.org/tutorials/colors/colormaps.html
* https://www.mathworks.com/help/matlab/ref/colormap.html
* http://colorbrewer2.org
* http://matplotlib.org/cmocean/
* http://jdherman.github.io/colormap/

Color

* http://htmlcolorcodes.com

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
    starttime='2000-10-20 00:00:00'; endtime='2100-10-20 23:00:00'; plottype='line'; cmapcolors='blue'; sizestyle='medium';
case 2
    endtime='2100-10-20 23:00:00'; plottype='line'; cmapcolors='blue'; sizestyle='medium';
case 3
    plottype='line'; cmapcolors='blue'; sizestyle='medium';
case 4
    cmapcolors='blue'; sizestyle='medium';
case 5
    sizestyle='medium';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(x)==1
    x=x';
end


%--------------------------------------------------------------------------
%Set plot drawing size

if strcmp(plottype,'line')==1 | strcmp(plottype,'line_grid')==1 ...
    | strcmp(plottype,'line_subplot')==1 | strcmp(plottype,'line_subplot_grid')==1

    if strcmp(sizestyle,'small')==1
        drawsize=1; %Size of plot
    elseif strcmp(sizestyle,'medium')==1
        drawsize=2; %Size of plot
    elseif strcmp(sizestyle,'large')==1
        drawsize=4; %Size of plot
    end

elseif strcmp(plottype,'bar')==1 | strcmp(plottype,'bar_grid')==1 ...
    | strcmp(plottype,'barh')==1 | strcmp(plottype,'barh_grid')==1

    if strcmp(sizestyle,'small')==1
        drawsize=0.5; %Size of plot
    elseif strcmp(sizestyle,'medium')==1
        drawsize=0.8; %Size of plot
    elseif strcmp(sizestyle,'large')==1
        drawsize=1; %Size of plot
    end

end

plotfontsize=18; %Size of plot fonts
axisfontsize=18; %Size of axis fonts


%--------------------------------------------------------------------------
%Input variable names

%if isempty(y)==1
%    xaxislabel='x'; %x axis label
%    yaxislabel=inputname(1); %y axis label
%else
%    xaxislabel=inputname(1); %x axis label
%    yaxislabel=inputname(2); %y axis label
%end 

xaxislabel='time'; %x axis label
yaxislabel='y'; %y axis label


%--------------------------------------------------------------------------
%Check input data

y=[];
%Create x and y if y=[]
if isempty(y)==1
    [Mx,Nx]=size(x);
    y=x;
    x=[1:1:Mx]';
end

%Input data shape
[Mx,Nx]=size(x);
[My,Ny]=size(y);


%--------------------------------------------------------------------------
%Reshape input data to (M,N)

%x and y have the same shape of (M,N)
if Mx==My & Nx==Ny
    xdata=x;
    ydata=y;

%x has a shape of (M,1) and y has a shape of (M,N)
elseif Mx==My & Nx==1
    xdata=repmat(x,1,Ny);
    ydata=y;

%x has a shape of (M,1) and y has a shape of (N,M)
elseif Mx==Ny & Nx==1
    xdata=repmat(x,1,My);
    ydata=y';

%x has a shape of (1,N) and y has a shape of (M,N)
elseif Mx==1 & Nx==Ny
    xdata=x';
    xdata=repmat(xdata,1,My);
    ydata=y';

%x has a shape of (1,N) and y has a shape of (N,M)
elseif Mx==1 & Nx==My
    xdata=x';
    xdata=repmat(xdata,1,Ny);
    ydata=y;

end

%Shape of reshaped data
[Mxdata,Nxdata]=size(xdata);
[Mydata,Nydata]=size(ydata);


%--------------------------------------------------------------------------
%Colors

light_blue=[212, 230, 241];
mid_blue=[41, 128, 185];
dark_blue=[21, 67, 96];

light_red=[242, 215, 213];
mid_red=[192, 57, 43];
dark_red=[100, 30, 22];

light_green=[200, 230, 201];
mid_green=[76, 175, 80];
dark_green=[27, 94, 32];

light_yellow=[255, 249, 196];
mid_yellow=[255, 235, 59];
dark_yellow=[245, 127, 23];

light_purple=[235, 222, 240];
mid_purple=[155, 89, 182];
dark_purple=[81, 46, 95];

light_brown=[246, 221, 204];
mid_brown=[211, 84, 0];
dark_brown=[110, 44, 0];

light_gray=[250, 250, 250];
dark_gray=[1, 1, 1];


%--------------------------------------------------------------------------
%Colormaps

%zclmap
%'blue': blue colormap
if strcmp(cmapcolors,'blue')==1
    zclmap=[light_blue;mid_blue;dark_blue]./255;

%'red': red colormap
elseif strcmp(cmapcolors,'red')==1
    zclmap=[light_red;mid_red;dark_red]./255;

%'green': green colormap
elseif strcmp(cmapcolors,'green')==1
    zclmap=[light_green;mid_green;dark_green]./255;

%'yellow': yellow colormap
elseif strcmp(cmapcolors,'yellow')==1
    zclmap=[light_yellow;mid_yellow;dark_yellow]./255;

%'purple': purple colormap
elseif strcmp(cmapcolors,'purple')==1
    zclmap=[light_purple;mid_purple;dark_purple]./255;

%'brown': brown colormap, 'cooper'
elseif strcmp(cmapcolors,'brown')==1
    zclmap=[light_brown;mid_brown;dark_brown]./255;

%'gray': gray colormap, 'gray'
elseif strcmp(cmapcolors,'gray')==1
    zclmap=[light_gray;dark_gray]./255;

%'blue_red': blue-red colormap, 'cool'
elseif strcmp(cmapcolors,'blue_red')==1
    zclmap=[mid_blue;mid_red]./255;

%'red_blue': red-blue colormap, 'cool-1'
elseif strcmp(cmapcolors,'red_blue')==1
    zclmap=[mid_red;mid_blue]./255;

%'blue_green': blue-green colormap, 'winter'
elseif strcmp(cmapcolors,'blue_green')==1
    zclmap=[mid_blue;mid_green]./255;

%'green_blue': green-blue colormap, 'winter-1'
elseif strcmp(cmapcolors,'green_blue')==1
    zclmap=[mid_green;mid_blue]./255;

%'green_yellow': green-yellow colormap, 'summer'
elseif strcmp(cmapcolors,'green_yellow')==1
    zclmap=[mid_green;mid_yellow]./255;

%'yellow_green': yellow-green colormap, 'summer-1'
elseif strcmp(cmapcolors,'yellow_green')==1
    zclmap=[mid_yellow;mid_green]./255;

%'red_yellow': red-yellow colormap, 'autumn'
elseif strcmp(cmapcolors,'red_yellow')==1
    zclmap=[mid_red;mid_yellow]./255;

%'yellow_red': yellow-red colormap, 'autumn-1'
elseif strcmp(cmapcolors,'yellow_red')==1
    zclmap=[mid_yellow;mid_red]./255;

%'cyclic': cyclic/oscillation colormap, 'hsv'
elseif strcmp(cmapcolors,'cyclic')==1
    %colormap('hsv')
    %zclmap=[mid_red;light_red;light_blue;mid_blue;mid_blue;light_blue;light_red;mid_red]./255;
    zclmap=[light_red;light_blue;mid_blue;mid_red]./255;

%'seq': sequential colormap, 'line'
elseif strcmp(cmapcolors,'seq')==1
    zclmap=lines(Nxdata);

%User input colormap
elseif isnumeric(cmapcolors)==1
    zclmap=cmapcolors./255; %Converting RGB to values between 0 and 1
    zclmap(zclmap<0 )=0;
    zclmap(zclmap>1 )=1;

%Pre-defined colormap such as 'cool', 'winter', etc
else
    zclmap=colormap(cmapcolors);

end


%--------------------------------------------------------------------------
%Interpolating a colormap to include ncolor zclmap

if Nxdata==1
    ncolor=3; %Number of zclmap to be used in colormap
else
    ncolor=Nxdata; %Number of zclmap to be used in colormap
end

if strcmp(cmapcolors,'seq')==1
    cmap_ncolor=lines(ncolor);
else
    cmap_ncolor=zeros(ncolor,3);
    cmap_ncolor(:,1:3)=interp1(linspace(1,ncolor,length(zclmap(:,1))),zclmap(:,1:3),linspace(1,ncolor,ncolor));
end

if Nxdata==1
    cmap_ncolor=cmap_ncolor(2:end,:); %Remove first (light) color
end


%--------------------------------------------------------------------------
%Generate time vector

timeformat='yyyy-mm-dd HH:MM:SS'; %Time format

startdate=datenum(starttime,timeformat); %Convert start date and time to numerical date
enddate=datenum(endtime,timeformat); %Convert end date and time to numerical date

timedata(:,1)=linspace(startdate,enddate,Mxdata); %Time seris in numerical date format
dt=timedata(2,1)-timedata(1,1); %Delta t

datestring=datestr(timedata,timeformat); %Time seris in string format


%--------------------------------------------------------------------------
%The line plot

if strcmp(plottype,'line')==1 | strcmp(plottype,'line_grid')==1

    hold on
    for i=1:Nxdata
        plot(timedata(:,1),ydata(:,i),'Color',cmap_ncolor(i,:),'LineWidth',drawsize)
    end

end


%--------------------------------------------------------------------------
%The line plot with subplots

if strcmp(plottype,'line_subplot')==1 | strcmp(plottype,'line_subplot_grid')==1

    plotfontsize=18; %Size of plot fonts
    axisfontsize=18; %Size of axis fonts
    hold on
    for i=1:Nxdata
        subplot(Nxdata,1,i)
        plot(timedata(:,1),ydata(:,i),'Color',cmap_ncolor(i,:),'LineWidth',drawsize)
        box on
        title(num2str(i))
        datetick('x','yyyy-mm-dd','keepticks')
        xlim([timedata(1,1),timedata(end,1)])
        ylim([nanmin(y(:)),nanmax(y(:))])
        set(gca,'fontsize',plotfontsize)
        xlabel(xaxislabel,'fontsize',axisfontsize)
        ylabel(yaxislabel,'fontsize',axisfontsize)
    end

end


%--------------------------------------------------------------------------
%The line plot with subplots

if strcmp(plottype,'line_subplot_grid')==1

    plotfontsize=18; %Size of plot fonts
    axisfontsize=18; %Size of axis fonts
    hold on
    for i=1:Nxdata
        subplot(Nxdata,1,i)
        plot(timedata(:,1),ydata(:,i),'Color',cmap_ncolor(i,:),'LineWidth',drawsize)
        box on
        grid on
        title(num2str(i))
        datetick('x','yyyy-mm-dd','keepticks')
        xlim([timedata(1,1),timedata(end,1)])
        ylim([nanmin(y(:)),nanmax(y(:))])
        set(gca,'fontsize',plotfontsize)
        xlabel(xaxislabel,'fontsize',axisfontsize)
        ylabel(yaxislabel,'fontsize',axisfontsize)
    end

end


%--------------------------------------------------------------------------
%Plot bar plot

if strcmp(plottype,'bar')==1 | strcmp(plottype,'bar_grid')==1

    if Nxdata==1
        bar(timedata,ydata,0.8,'FaceColor',cmap_ncolor(1,:),'EdgeColor','none')
    else
        bar(timedata,ydata(:,1),0.8,'FaceColor',cmap_ncolor(1,:),'EdgeColor','none')
    end

end


%--------------------------------------------------------------------------
%Plot bar plot

%The bar plot
if strcmp(plottype,'barh')==1 | strcmp(plottype,'barh_grid')==1

    if Nxdata==1
        barh(timedata,ydata,0.8,'FaceColor',cmap_ncolor(1,:),'EdgeColor','none')
    else
        barh(timedata,ydata(:,1),0.8,'FaceColor',cmap_ncolor(1,:),'EdgeColor','none')
    end

end


%--------------------------------------------------------------------------
%Set plot properties

if strcmp(plottype,'line_subplot')==0 & strcmp(plottype,'line_subplot_grid')==0

    %Set time axis ticks
    %set(gca,'XTick',ceil(startdate):dt:floor(enddate))
    if strcmp(plottype,'barh')==1 | strcmp(plottype,'barh_grid')==1
        datetick('y','yyyy-mm-dd','keepticks')
    else
        datetick('x','yyyy-mm-dd','keepticks')
    end

    %Draw outline around a plot
    box on

    %Set axis limits
    if strcmp(plottype,'barh')==1 | strcmp(plottype,'barh_grid')==1
        xlim([nanmin(y(:)),nanmax(y(:))])
        ylim([timedata(1,1),timedata(end,1)])
    else
        xlim([timedata(1,1),timedata(end,1)])
        ylim([nanmin(y(:)),nanmax(y(:))])
    end

    %Plot grid lines
    if strcmp(plottype,'line_grid')==1 | strcmp(plottype,'line_subplot_grid')==1 ...
        | strcmp(plottype,'bar_grid')==1 | strcmp(plottype,'barh_grid')==1
        grid on
    end

    %Set figure font size
    set(gca,'fontsize',plotfontsize)

    %Plot axis label
    if strcmp(plottype,'barh')==1 | strcmp(plottype,'barh_grid')==1
        xlabel(yaxislabel,'fontsize',axisfontsize)
        ylabel(xaxislabel,'fontsize',axisfontsize)
    else
        xlabel(xaxislabel,'fontsize',axisfontsize)
        ylabel(yaxislabel,'fontsize',axisfontsize)
    end

    %Plot legends
    legend_all_names={'y_1','y_2','y_3','y_4','y_5','y_6','y_7','y_8','y_9'};
    if Nxdata>=2 & Nxdata<=3
        legend_names=legend_all_names(1:Nxdata);
        legend(legend_names)
    end

end

%--------------------------------------------------------------------------
