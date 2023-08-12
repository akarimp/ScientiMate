function [xReplaced, outlier_Indx] = replaceoutlier(x, WindowSize, zscore_threshold, interpMethod, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2020-05-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

replaceoutlier
==============

.. code:: MATLAB

    [xReplaced, outlier_Indx] = replaceoutlier(x, WindowSize, zscore_threshold, interpMethod, dispout)

Description
-----------

Remove outliers in the time series using moving z-score window

Inputs
------

x
    Input data
WindowSize=15;
    Window size (number of adjacent elements) that is used for moving window, should be equal or larger than 3
zscore_threshold=2;
    | z-score threshold to define outliers
    | data in range of x < (xmean-zscore_threshold*std) or x > (xmean+zscore_threshold*std) considered outliers
interpMethod='linear';
    | Interpolation method for replacing spike points:
    | Matlab/Octave: 'linear', 'nearest', 'next', 'previous', 'pchip', 'cubic', 'spline'
    | Python/Scipy : 'linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic'
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

xReplaced
    Replaced data
outlier_Indx
    Logical index of replaced points

Examples
--------

.. code:: MATLAB

    fs=128;
    t(:,1)=linspace(0,9.5,10*fs);
    x(:,1)=sin(2*pi*0.3*t)+0.1*sin(2*pi*4*t);
    spikeloc(:,1)=[10:100:length(t(:,1))];
    x(spikeloc+round(2*randn(length(spikeloc(:,1)),1)),1)=sign(randn(length(spikeloc(:,1)),1));
    x(220:225,1)=1.5;
    x=x+5;
    [xReplaced,outlier_Indx]=replaceoutlier(x,37,2,'linear','yes');

    fs=2;
    t(:,1)=linspace(0,1023.5,1024*fs);
    x(:,1)=detrend(0.5.*cos(2*pi*0.2*t)+(-0.1+(0.1-(-0.1))).*rand(1024*fs,1));
    spikeloc(:,1)=[10:100:length(t(:,1))];
    x(spikeloc+round(2*randn(length(spikeloc(:,1)),1)),1)=sign(randn(length(spikeloc(:,1)),1));
    [xReplaced,outlier_Indx]=replaceoutlier(x,21,2,'linear','yes');

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
        WindowSize=15; zscore_threshold=2; interpMethod='linear'; dispout='no';
    case 2
        zscore_threshold=2; interpMethod='linear'; dispout='no';
    case 3
        interpMethod='linear'; dispout='no';
    case 4
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(x)==1
    x=x';
end

%--------------------------------------------------------------------------

%Preserving an input data
xInput=x;

%--------------------------------------------------------------------------

N=length(x(:,1)); %Length of the input data
Nw=WindowSize; %Window size

%Checking size of the window
Nw(Nw>N)=N;
Nw(Nw<=3)=3;

%--------------------------------------------------------------------------
%Removing outliers

%Assigning empty to NaN_Indx
NaN_Indx=[];

%Moving mean
rolled_mean=movmean(x,Nw,'Endpoints','shrink');

%Moving std
rolled_std=movstd(x,Nw,'Endpoints','shrink');

%Calculate moving z-score
rolled_z_score=(x-rolled_mean)./rolled_std;

%Replace outliers with NaN
xInput(abs(rolled_z_score)>zscore_threshold)=NaN;

%Find NaN
NaN_Indx=(isnan(xInput)==1);
Valid_Indx=(isnan(xInput)~=1);

%Assigning x to xReplaced
xReplaced=x;

%Replacing outliers
samples(:,1)=[1:1:length(x(:,1))];
if sum(NaN_Indx(:))>0
    xReplaced(NaN_Indx==1)=interp1(samples(Valid_Indx==1),xReplaced(Valid_Indx==1),samples(NaN_Indx==1),interpMethod,'extrap');
end

%Assigning NaN_Indx to outlier_Indx
outlier_Indx=NaN_Indx;

%--------------------------------------------------------------------------
%Defining replaced points

if sum(NaN_Indx(:))>0
    NumberReplacedPoint=sum(NaN_Indx(:));
    PercentReplacedPoint=sum(NaN_Indx(:))/numel(x)*100;
else
    NumberReplacedPoint=0;
    PercentReplacedPoint=0;
end

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1

    disp(['Total number of pints replaced = ',num2str(NumberReplacedPoint)])
    disp(['Percent of replaced points     = ',num2str(PercentReplacedPoint),' %'])
    
    subplot(2,1,1)
    plot(samples,x)
    xlabel('Sample')
    ylabel('Data')
    title('Input Data')

    subplot(2,1,2)
    plot(samples,xReplaced)
    hold on
    scatter(samples(NaN_Indx),xReplaced(NaN_Indx))
    xlabel('Sample')
    ylabel('Data')
    title('Replaced Data')
    
end

%--------------------------------------------------------------------------
