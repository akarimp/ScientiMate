function [xDespiked, Indx] = replacespikediff(x, WindowSize, spikedefineMethod, nrpeat, DespikeScaleFactor, interpMethod, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-02-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

replacespikediff
================

.. code:: MATLAB

    [xDespiked, Indx] = replacespikediff(x, WindowSize, spikedefineMethod, nrpeat, DespikeScaleFactor, interpMethod, dispout)

Description
-----------

Remove spikes in the time series using a local difference of data respect to a moving average window

Inputs
------

x
    Input data
WindowSize=5;
    Window size (number of adjacent elements) that is used for smoothing, should be equal or larger than 3
spikedefineMethod='ellipse';
    | Method to define spike points
    | 'ellipse': use both local difference and its gradient
    | 'circle': use local difference only
nrpeat=1;
    Number of time despiking procedure is repeating
DespikeScaleFactor=1;
    Scaling a despiking threshold
interpMethod='linear';
    | Interpolation method for replacing spike points:
    | Matlab/Octave: 'linear', 'nearest', 'next', 'previous', 'pchip', 'cubic', 'spline'
    | Python/Scipy : 'linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic'
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

xDespiked
    Dispiked data
Indx
    Index of despiked points

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
    [xDespiked,Indx]=replacespikediff(x,5,'ellipse',2,1,'linear','yes');

    fs=2;
    t(:,1)=linspace(0,1023.5,1024*fs);
    x(:,1)=detrend(0.5.*cos(2*pi*0.2*t)+(-0.1+(0.1-(-0.1))).*rand(1024*fs,1));
    spikeloc(:,1)=[10:100:length(t(:,1))];
    x(spikeloc+round(2*randn(length(spikeloc(:,1)),1)),1)=sign(randn(length(spikeloc(:,1)),1));
    [xDespiked,Indx]=replacespikediff(x,5,'ellipse',1,1,'linear','yes');

References
----------

Goring, D. G., & Nikora, V. I. (2002). 
Despiking acoustic Doppler velocimeter data. 
Journal of Hydraulic Engineering, 128(1), 117-126.

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
        WindowSize=5; spikedefineMethod='ellipse'; nrpeat=1; DespikeScaleFactor=1; interpMethod='linear'; dispout='no';
    case 2
        spikedefineMethod='ellipse'; nrpeat=1; DespikeScaleFactor=1; interpMethod='linear'; dispout='no';
    case 3
        nrpeat=1; DespikeScaleFactor=1; interpMethod='linear'; dispout='no';
    case 4
        DespikeScaleFactor=1; interpMethod='linear'; dispout='no';
    case 5
         interpMethod='linear'; dispout='no';
    case 6
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

%Creating window function in time domain

nwin(:,1)=[0:1:Nw-1]; %Counter from 0 to Nw-1

%Calculating window function
WnRec=ones(Nw,1); %Rectangular window function (moving average filter)
Wn=WnRec;
WnNorm=Wn./sum(Wn(:,1)); %Normalizing a window function in time domain

%--------------------------------------------------------------------------
%Removing spikes

for i=1:nrpeat
    
    %Smoothing input using convolution
    xSmoothed=conv(x,WnNorm,'same'); %Smoothing input using convolution

    %Removing alias from edge of smooted data
    for n=1:floor(Nw/2)
        WnEdge=WnNorm(ceil(Nw/2)-(n-1):ceil(Nw/2)+(n-1),1);
        WnEdgeNorm=WnEdge./sum(WnEdge(:,1)); %Normalizing a window function in time domain
        xSmoothed(n,1)=sum(x(1:n+(n-1),1).*WnEdgeNorm);
        xSmoothed(N-n+1,1)=sum(x((N-n+1)-(n-1):N,1).*WnEdgeNorm);
    end
    
    %Local difference
    xDiff=abs(x-xSmoothed); %Local difference of x respecte to moving average window

    %Standard deviation
    xStd=std(x); %Standard deviation of x
    xDiffStd=std(xDiff); %Standard deviation of xDiff

    %Expected absolute maximum based on normal distributaion
    Lambda=sqrt(2*log(N));

    %Axis of ellipse
    a1=Lambda*xStd; %Major axis of ellipse for xDiff vs x
    b1=Lambda*xDiffStd; %Minor axis of ellipse for xDiff vs x

    %Calculating mean values
    xmean=mean(x(:,1));
    xDiffmean=mean(xDiff(:,1));

    %Calculating distances based on ellipse equation
    distellipse=((x-xmean).^2./a1.^2+(xDiff-xDiffmean).^2./b1.^2);

    %Defining spike points 
    if strcmp(spikedefineMethod,'ellipse')==1
        %If a point is locatated outside an ellipse, it is considered as spike
        Indx1=find(distellipse>1*DespikeScaleFactor);
    elseif strcmp(spikedefineMethod,'circle')==1
        %If a point is locatated outside a circle, it is considered as spike
        %Considering a Normal distributaion, all points are smaller than 3 times of standard deviation
        Indx1=find(xDiff>(xDiffmean+3*xDiffStd)*DespikeScaleFactor);
    end
    
    %Setting spike points as NaN
    xDespiked=x;
    xDespiked(Indx1,1)=NaN;

    %Locating all spike pints
    Indx2=find(isnan(xDespiked(:,1))==1);

    %Locating all non-spike pints
    Indx3=find(isnan(xDespiked(:,1))~=1);

    %Replacing spike points
    samples(:,1)=[1:1:length(x(:,1))];
    if length(Indx2)~=0
        xDespiked(Indx2,1)=interp1(samples(Indx3,1),xDespiked(Indx3,1),samples(Indx2,1),interpMethod,'extrap');
    end

    %Assigning a despiked data to x for next repeat
    x=xDespiked;

end

%--------------------------------------------------------------------------
%Defining despiked points

Indx(:,1)=find(xInput~=xDespiked);
NumberDespikedPoint=length(Indx(:,1));
PercentDespikedPoint=length(Indx(:,1))/length(x(:,1))*100;

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1

    disp(['Total number of pints despiked = ',num2str(NumberDespikedPoint)])
    disp(['Percent of despiked points     = ',num2str(PercentDespikedPoint),' %'])
    
    subplot(2,1,1)
    plot(samples,xInput)
    hold on
    scatter(samples(Indx,1),xInput(Indx,1))
    xlabel('Sample')
    ylabel('Data')
    title('Input Data')

    subplot(2,1,2)
    plot(samples,xDespiked)
    hold on
    scatter(samples(Indx,1),xDespiked(Indx,1))
    xlabel('Sample')
    ylabel('Data')
    title('Despiked Data')
    
end

%--------------------------------------------------------------------------
