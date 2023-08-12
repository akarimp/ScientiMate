function [xReplaced, NaN_Indx] = replacemissing1d(x, what2replace, interpMethod, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-02-01/2020-02-01                  +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

replacemissing1d
================

.. code:: MATLAB

    [xReplaced, NaN_Indx] = replacemissing1d(x, what2replace, interpMethod, dispout)

Description
-----------

Replace missing data points in 1d data such as time series

Inputs
------

x
    Input data
what2replace='both';
    | What needs to be replaced
    | 'NaN': replacing NaN data points
    | 'Inf': replacing Inf data points
    | 'both': replacing NaN and Inf data points
    | Number: replacing data points equal to Number
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
NaN_Indx
    Logical index of replaced points

Examples
--------

.. code:: MATLAB

    fs=128;
    t(:,1)=linspace(0,9.5,10*fs);
    x(:,1)=sin(2*pi*0.3*t)+0.1*sin(2*pi*4*t);
    spikeloc(:,1)=[10:100:length(t(:,1))];
    x(spikeloc+round(2*randn(length(spikeloc(:,1)),1)),1)=NaN;
    x=x+5;
    x(220:225,1)=NaN;
    [xReplaced,NaN_Indx]=replacemissing1d(x,'NaN','linear','yes');

    fs=2;
    t(:,1)=linspace(0,1023.5,1024*fs);
    x(:,1)=detrend(0.5.*cos(2*pi*0.2*t)+(-0.1+(0.1-(-0.1))).*rand(1024*fs,1));
    spikeloc(:,1)=[10:100:length(t(:,1))];
    x=x+1;
    x(spikeloc+round(2*randn(length(spikeloc(:,1)),1)),1)=0;
    [xReplaced,NaN_Indx]=replacemissing1d(x,0,'linear','yes');

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
        what2replace='both'; interpMethod='linear'; dispout='no';
    case 2
        interpMethod='linear'; dispout='no';
    case 3
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
%Replacing missing data points

%Assigning empty to NaN_Indx
NaN_Indx=[];

%Defining missing data points
if strcmp(what2replace,'NaN')==1
    %Define NaN data points
    NaN_Indx=(isnan(x)==1);
    Valid_Indx=(isnan(x)~=1);
elseif strcmp(what2replace,'Inf')==1
    %Define Inf data points
    NaN_Indx=(isinf(x)==1);
    Valid_Indx=(isinf(x)~=1);
elseif strcmp(what2replace,'both')==1
    %Define NaN and Inf data points
    NaN_Indx=(isnan(x)==1 | isinf(x)==1);
    Valid_Indx=(isnan(x)~=1 | isinf(x)~=1);
elseif isnumeric(what2replace)==1
    %Define NaN and Inf data points
    NaN_Indx=(x==what2replace);
    Valid_Indx=(x~=what2replace);
end

%Assigning x to xReplaced
xReplaced=x;

%Replacing missing points
samples(:,1)=[1:1:length(x(:,1))];
if sum(NaN_Indx(:))>0
    xReplaced(NaN_Indx==1)=interp1(samples(Valid_Indx==1),xReplaced(Valid_Indx==1),samples(NaN_Indx==1),interpMethod,'extrap');
end
    
%Assigning a replaced data to x for next repeat
%x=xReplaced;


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
