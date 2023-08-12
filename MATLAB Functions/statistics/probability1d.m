function [fdensity, fdensitycumulative, bincenter, xmean, xstd] = probability1d(x, binedge, dispout)
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

probability1d
=============

.. code:: MATLAB

    [fdensity, fdensitycumulative, bincenter, xmean, xstd] = probability1d(x, binedge, dispout)

Description
-----------

Calculate 1D probability density distribution for a given dataset

Inputs
------

x
    Input data 
binedge
    | Bin edges  
    | length(binedge)=number of bin +1   
    | If there are N bins in histogram/distribution, then values in binedge are as:   
    | 1st value: left edge of 1st bin, 2nd value: left edge of 2nd bin, ...   
    | (N)th value: left edge of last bin, (N+1)th value: right edge of last bin   
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

fdensity
    Probability density distribution
fdensitycumulative
    Cumulative probability density distribution
bincenter
    Bin center
xmean
    Mean value of input data
xstd
    Standard deviation of input data

Examples
--------

.. code:: MATLAB

    x(:,1)=(-0.1+(0.1-(-0.1))).*randn(1024*2,1);
    binedge(:,1)=linspace(min(x),max(x),11);
    [fdensity,fdensitycumulative,bincenter,xmean,xstd]=probability1d(x,binedge,'yes');

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
        binedge(:,1)=linspace(nanmin(x(:,1)),nanmax(x(:,1)),11); dispout='no';
    case 2
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(x)==1
    x=x';
end

if isrow(binedge)==1
    binedge=binedge';
end

%--------------------------------------------------------------------------
%Removing NaN and Inf values

%Total number of elements in dataset
NTotal=length(x(:,1));

%Finding NaN and Inf in x dataset
Indx1=find(isnan(x(:,1))~=1 & isinf(x(:,1))~=1);
x1=x(Indx1,1);

%Datasets without NaN and Inf
x=x1;

%Total number of elements in dataset after removal of NaN and Inf
N=length(x(:,1));

%Total number of NaN and Inf elements in dataset
NNaNInf=NTotal-N;

%--------------------------------------------------------------------------
%Calculating probability density distribution

binwidth(:,1)=diff(binedge(:,1)); %Bin width, length(binwidth) is one element less than length(binedge) 
bincenter(:,1)=binedge(1:end-1,1)+binwidth./2; %Bin center

[NPoints(:,1),bincenterout(:,1)]=hist(x,bincenter); %Calculating number of data points in each bin
histarea=NPoints.*binwidth; %Area under each element of histogram
fdensity=NPoints./(sum(histarea)); %Calculating probability density distribution
%Note: sum(fdensity.*binwidth)=1

%--------------------------------------------------------------------------
%Calculating cumulative probability density distribution

fdensitycumulative=zeros(length(bincenter(:,1)),1); %Pre-assigning vector
for i=1:length(bincenter(:,1))
    fdensitycumulative(i,1)=sum(fdensity(1:i,1).*binwidth(1:i,1));
end

%--------------------------------------------------------------------------
%Calculating statistical properties

xmean=mean(x(:,1)); %Mean value of the input data
xstd=std(x(:,1)); %Standard deviation of the input data

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1 
    
    val=[xmean xstd NTotal N NNaNInf];
    name={'mean','st. dev.','N Input','N Used','N NaN Inf'};
    for i=1:length(val)
        %fprintf('%14s   %g\n',name{i},val(i));
        fprintf('     %-15s     %-+g\n',name{i},val(i));
    end

    %Plotting
    subplot(2,1,1)
    bar(bincenter,fdensity)
    title('Probability Density Distribution')
    ylabel('Probability Density, f')
    xlabel('Data')
    
    subplot(2,1,2)
    plot(bincenter,fdensitycumulative)
    title('Cumulative Probability Density Distribution')
    ylabel('Cumulative Probability Density, f')
    xlabel('Data')
    
end

%--------------------------------------------------------------------------
