function [fdensityxy, fdensityx, fdensityy, fdensitycumulativex, fdensitycumulativey, bincenterx, bincentery, xmean, ymean, xstd, ystd] = probability2d(x, y, binedgex, binedgey, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-06-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

probability2d
=============

.. code:: MATLAB

    [fdensityxy, fdensityx, fdensityy, fdensitycumulativex, fdensitycumulativey, bincenterx, bincentery, xmean, ymean, xstd, ystd] = probability2d(x, y, binedgex, binedgey, dispout)

Description
-----------

Calculate 2D (joint) probability density distribution for two given datasets

Inputs
------

x
    First input dataset 
y
    Second input dataset 
binedgex
    | Bin edges for x data  
    | length(binedgex)=number of bin +1   
    | If there are N bins in histogram/distribution, then values in binedgex are as:   
    | 1st value: left edge of 1st bin, 2nd value: left edge of 2nd bin, ...   
    | (N)th value: left edge of last bin, (N+1)th value: right edge of last bin   
binedgey
    | Bin edge for y data  
    | length(binedgey)=number of bin +1   
    | If there are N bins in histogram/distribution, then values in binedgey are as:   
    | 1st value: left edge of 1st bin, 2nd value: left edge of 2nd bin, ...   
    | (N)th value: left edge of last bin, (N+1)th value: right edge of last bin   
dispout='no';
    | Define to display outputs or not
    | 'no': not display 
    | 'bar': bar plot
    | 'imagesc': 2 dimensional plot using imagesc or imshow
    | 'pcolor': 2 dimensional plot using pcolor
    | 'contour': 2 dimensional contour plot, number of contour=32
    | 'contourf': 2 dimensional filled contour plot, number of contour=32
    | 'surface': 3 dimensional surface plot 
    | 'binscatter': 2 dimensional histogram plot using imagesc or imshow

Outputs
-------

fdensityxy
    2D (joint) probability density distribution for x-y data
fdensityx
    Probability density distribution for x data
fdensityy
    Probability density distribution for y data
fdensitycumulativex
    Cumulative probability density distribution for x data
fdensitycumulativey
    Cumulative probability density distribution for y data
bincenterx
    Bin center for x data
bincentery
    Bin center for y data
xmean
    Mean value of x data
ymean
    Mean value of y data
xstd
    Standard deviation of x data
ystd
    Standard deviation of y data

Examples
--------

.. code:: MATLAB

    x(:,1)=(-0.1+(0.1-(-0.1))).*randn(1024*2,1);
    y(:,1)=(-0.2+(0.2-(-0.2))).*randn(1024*2,1);
    binedgex(:,1)=linspace(min(x),max(x),11);
    binedgey(:,1)=linspace(min(y),max(y),11);
    [fdensityxy,fdensityx,fdensityy,fdensitycumulativex,fdensitycumulativey,bincenterx,bincentery,xmean,ymean,xstd,ystd]=probability2d(x,y,binedgex,binedgey,'surface');

    x(:,1)=randn(100000,1);
    y(:,1)=1.5.*x+randn(100000,1);
    binedgex(:,1)=linspace(min(x),max(x),101);
    binedgey(:,1)=linspace(min(y),max(y),101);
    probability2d(x,y,binedgex,binedgey,'binscatter');

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
        binedgex(:,1)=linspace(nanmin(x(:,1)),nanmax(x(:,1)),11); binedgey(:,1)=linspace(nanmin(y(:,1)),nanmax(y(:,1)),11); dispout='no';
    case 3
        binedgey(:,1)=linspace(nanmin(y(:,1)),nanmax(y(:,1)),11); dispout='no';
    case 4
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(x)==1
    x=x';
end

if isrow(y)==1
    y=y';
end

if isrow(binedgex)==1
    binedgex=binedgex';
end

if isrow(binedgey)==1
    binedgey=binedgey';
end

%--------------------------------------------------------------------------
%Removing NaN and Inf values

%Total number of elements in dataset
NTotal=length(x(:,1));

%Finding NaN and Inf in x dataset
Indx1=find(isnan(x(:,1))~=1 & isinf(x(:,1))~=1);
x1=x(Indx1,1);
y1=y(Indx1,1);

%Finding NaN and Inf in y dataset
Indx2=find(isnan(y1(:,1))~=1 & isinf(y1(:,1))~=1);
x2=x1(Indx2,1);
y2=y1(Indx2,1);

%Datasets without NaN and Inf
x=x2;
y=y2;

%Total number of elements in dataset after removal of NaN and Inf
N=length(x(:,1));

%Total number of NaN and Inf elements in dataset
NNaNInf=NTotal-N;

%--------------------------------------------------------------------------
%Calculating probability density distribution for x data

binwidthx(:,1)=diff(binedgex(:,1)); %Bin width for x data, length(binwidth) is one element less than length(binedge) 
bincenterx(:,1)=binedgex(1:end-1,1)+binwidthx./2; %Bin center for x data

[NPointsx(:,1),bincenteroutx(:,1)]=hist(x,bincenterx); %Calculating number of data points in each bin for x data
histareax=NPointsx.*binwidthx; %Area under each element of histogram for x data
fdensityx=NPointsx./(sum(histareax)); %Calculating probability density distribution for x data
%Note: sum(fdensityx.*binwidthx)=1

%--------------------------------------------------------------------------
%Calculating cumulative probability density distribution for x data

fdensitycumulativex=zeros(length(bincenterx(:,1)),1); %Pre-assigning vector
for i=1:length(bincenterx(:,1))
    fdensitycumulativex(i,1)=sum(fdensityx(1:i,1).*binwidthx(1:i,1));
end

%--------------------------------------------------------------------------
%Calculating statistical properties for x data

xmean=mean(x(:,1)); %Mean value of x data
xstd=std(x(:,1)); %Standard deviation of x data

%--------------------------------------------------------------------------
%Calculating probability density distribution for y data

binwidthy(:,1)=diff(binedgey(:,1)); %Bin width for y data, length(binwidth) is one element less than length(binedge) 
bincentery(:,1)=binedgey(1:end-1,1)+binwidthy./2; %Bin center for y data

[NPointsy(:,1),bincenterouty(:,1)]=hist(y,bincentery); %Calculating number of data points in each bin for y data
histareay=NPointsy.*binwidthy; %Area under each element of histogram for y data
fdensityy=NPointsy./(sum(histareay)); %Calculating probability density distribution for y data
%Note: sum(fdensityy.*binwidthy)=1

%--------------------------------------------------------------------------
%Calculating cumulative probability density distribution for y data

fdensitycumulativey=zeros(length(bincentery(:,1)),1); %Pre-assigning vector
for i=1:length(bincentery(:,1))
    fdensitycumulativey(i,1)=sum(fdensityy(1:i,1).*binwidthy(1:i,1));
end

%--------------------------------------------------------------------------
%Calculating statistical properties for y data

ymean=mean(y(:,1)); %Mean value of y data
ystd=std(y(:,1)); %Standard deviation of y data

%--------------------------------------------------------------------------
%Calculating 2D (joint) probability density distribution for x and y

[NPointsxy,bincenteroutxy]=hist3([x y],{bincenterx bincentery});

[binwidthxGrid,binwidthyGrid]=meshgrid(binwidthx(:,1),binwidthy(:,1));
histxyvolume=NPointsxy.*binwidthxGrid'.*binwidthyGrid'; %Volume under each element of histogram
fdensityxy=NPointsxy./(sum(histxyvolume(:))); %Calculating 2D (joint) probability density distribution for x-y data
%Note: sum(sum(fdensityxy.*binwidthxGrid'.*binwidthyGrid'))=1

%fdensityxy=zeros(length(bincenterx(:,1)),length(bincentery(:,1))); %Pre-assigning vector
%for i=1:length(bincenterx(:,1))
%    for j=1:length(bincentery(:,1))
%        fdensityxy(i,j)=NPointsxy(i,j)/(sum(histxyvolume(:))); %Calculating 2D (joint) probability density distribution for x-y data
%        %Note: sum(sum(fdensityxy.*binwidthxGrid'.*binwidthyGrid'))=1
%    end
%end

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'bar')==1 %Bar plot
    
    val=[xmean xstd ymean ystd NTotal N NNaNInf];
    name={'x mean','x st. deviation','y mean','y st. dev.','N Input','N Used','N NaN Inf'};
    for i=1:length(val)
        %fprintf('%14s   %g\n',name{i},val(i));
        fprintf('     %-15s     %-+g\n',name{i},val(i));
    end
    
    %Plotting x Data
    subplot(2,2,1)
    bar(bincenterx,fdensityx)
    title('Probability Density Distribution')
    ylabel('Probability Density, fx')
    xlabel('x Data')
    
    subplot(2,2,3)
    plot(bincenterx,fdensitycumulativex)
    title('Cumulative Probability Density Distribution')
    ylabel('Cumulative Probability Density, fx')
    xlabel('x Data')
    
    %Plotting y Data
    subplot(2,2,2)
    bar(bincentery,fdensityy)
    title('Probability Density Distribution') %Title of plot
    ylabel('Probability Density, fy')
    xlabel('y Data')

    subplot(2,2,4)
    plot(bincentery,fdensitycumulativey)
    title('Cumulative Probability Density Distribution')
    ylabel('Cumulative Probability Density, fy')
    xlabel('y Data')
    
elseif strcmp(dispout,'imagesc')==1 %Imagesc plot
    
    val=[xmean xstd ymean ystd NTotal N NNaNInf];
    name={'x mean','x st. deviation','y mean','y st. dev.','N Input','N Used','N NaN Inf'};
    for i=1:length(val)
        fprintf('     %-15s     %-+g\n',name{i},val(i));
    end

    [bincenterxGrid,bincenteryGrid] = meshgrid(bincenterx(:,1),bincentery(:,1));
    imagesc([min(bincenterxGrid(:)),max(bincenterxGrid(:))],[min(bincenteryGrid(:)),max(bincenteryGrid(:))],fdensityxy')
    set(gca,'YDir','normal'); %Correct y axis direction
    colorbar
    title('2D Probability Density Distribution')
    xlabel('x Data')
    ylabel('y Data')
    
elseif strcmp(dispout,'pcolor')==1 %Pcolor plot

    val=[xmean xstd ymean ystd NTotal N NNaNInf];
    name={'x mean','x st. deviation','y mean','y st. dev.','N Input','N Used','N NaN Inf'};
    for i=1:length(val)
        fprintf('     %-15s     %-+g\n',name{i},val(i));
    end

    [bincenterxGrid,bincenteryGrid] = meshgrid(bincenterx(:,1),bincentery(:,1));
    pcol=pcolor(bincenterxGrid,bincenteryGrid,fdensityxy');
    set(pcol,'edgecolor','none')
    colorbar
    title('2D Probability Density Distribution')
    xlabel('x Data')
    ylabel('y Data')
    
elseif strcmp(dispout,'contour')==1 %Contour plot
    
    val=[xmean xstd ymean ystd NTotal N NNaNInf];
    name={'x mean','x st. deviation','y mean','y st. dev.','N Input','N Used','N NaN Inf'};
    for i=1:length(val)
        fprintf('     %-15s     %-+g\n',name{i},val(i));
    end

    [bincenterxGrid,bincenteryGrid] = meshgrid(bincenterx(:,1),bincentery(:,1));
    contour(bincenterxGrid,bincenteryGrid,fdensityxy');
    colorbar
    title('2D Probability Density Distribution')
    xlabel('x Data')
    ylabel('y Data')
    
elseif strcmp(dispout,'contourf')==1 %Contourf plot
    
    val=[xmean xstd ymean ystd NTotal N NNaNInf];
    name={'x mean','x st. deviation','y mean','y st. dev.','N Input','N Used','N NaN Inf'};
    for i=1:length(val)
        fprintf('     %-15s     %-+g\n',name{i},val(i));
    end

    [bincenterxGrid,bincenteryGrid] = meshgrid(bincenterx(:,1),bincentery(:,1));
    contourf(bincenterxGrid,bincenteryGrid,fdensityxy');
    colorbar
    title('2D Probability Density Distribution')
    xlabel('x Data')
    ylabel('y Data')
    
elseif strcmp(dispout,'surface')==1 %Surface plot
    
    val=[xmean xstd ymean ystd NTotal N NNaNInf];
    name={'x mean','x st. deviation','y mean','y st. dev.','N Input','N Used','N NaN Inf'};
    for i=1:length(val)
        fprintf('     %-15s     %-+g\n',name{i},val(i));
    end

    [bincenterxGrid,bincenteryGrid] = meshgrid(bincenterx(:,1),bincentery(:,1));
    surf1=surf(bincenterxGrid,bincenteryGrid,fdensityxy');
    %set(surf1,'xscale','log','zscale','log')
    set(surf1,'FaceColor','interp','EdgeColor','none')
    colorbar
    title('2D Probability Density Distribution')
    xlabel('x Data')
    ylabel('y Data')
    zlabel('2D Probability Density, f_{x,y}')
    
elseif strcmp(dispout,'binscatter')==1 %Binscatter plot
    
    val=[xmean xstd ymean ystd NTotal N NNaNInf];
    name={'x mean','x st. deviation','y mean','y st. dev.','N Input','N Used','N NaN Inf'};
    for i=1:length(val)
        fprintf('     %-15s     %-+g\n',name{i},val(i));
    end

    imagesc([min(bincenteroutxy{1}),max(bincenteroutxy{1})],[min(bincenteroutxy{2}),max(bincenteroutxy{2})],NPointsxy')
    set(gca,'YDir','normal'); %Correct y axis direction
    colorbar
    title('2D Histogram')
    xlabel('x Data')
    ylabel('y Data')

end

%--------------------------------------------------------------------------
