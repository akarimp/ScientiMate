function [xDespiked, Indx] = replacespikeenvelope(x, lowbound, upbound, nrpeat, interpMethod, dispout)
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

replacespikeenvelope
====================

.. code:: MATLAB

    [xDespiked, Indx] = replacespikeenvelope(x, lowbound, upbound, nrpeat, interpMethod, dispout)

Description
-----------

Remove spikes in the time series that are outside a defined envelope

Inputs
------

x
    Input data
lowbound
    Lower boundary of the data, all data points should be larger than that
upbound
    Upper boundary of the data, all data points should be smaller than that
nrpeat=1;
    Number of time despiking procedure is repeating
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
    x(220:225,1)=1.5;
    x=x+5;
    [xDespiked,Indx]=replacespikeenvelope(x,3.9,6.1,1,'linear','yes');

    fs=2;
    t(:,1)=linspace(0,1023.5,1024*fs);
    x(:,1)=detrend(0.5.*cos(2*pi*0.2*t)+(-0.1+(0.1-(-0.1))).*rand(1024*fs,1));
    spikeloc(:,1)=[10:100:length(t(:,1))];
    x(spikeloc+round(2*randn(length(spikeloc(:,1)),1)),1)=sign(randn(length(spikeloc(:,1)),1));
    [xDespiked,Indx]=replacespikeenvelope(x,-0.6,0.6,1,'linear','yes');

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
    case 3
        nrpeat=1; interpMethod='linear'; dispout='no';
    case 4
        interpMethod='linear'; dispout='no';
    case 5
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
%Removing spikes

for i=1:nrpeat
    
    %Checking the defined boundary
    Indx1=find(x>lowbound & x<upbound);
    if length(Indx1)==0
        disp('Warning-----------------------------------------')
        disp('No data points are between defined boundaries')
        disp('Data are not de-spiked')
        disp('------------------------------------------------')
        lowbound=min(x(:,1));
        upbound=max(x(:,1));
    end
    
    %Defining spike points, if a point locatated outside boundary limits is considered as spike
    Indx2=find(x(:,1)<lowbound);
    Indx3=find(x(:,1)>upbound);

    %Setting spike points as NaN
    xDespiked=x;
    xDespiked(Indx2,1)=NaN;
    xDespiked(Indx3,1)=NaN;

    %Locating all spike pints
    Indx4=find(isnan(xDespiked(:,1))==1);

    %Locating all non-spike pints
    Indx5=find(isnan(xDespiked(:,1))~=1);

    %Replacing spike points
    samples(:,1)=[1:1:length(x(:,1))];
    if length(Indx4)~=0
        xDespiked(Indx4,1)=interp1(samples(Indx5,1),xDespiked(Indx5,1),samples(Indx4,1),interpMethod,'extrap');
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
    plot([1;max(samples(:,1))],[lowbound;lowbound])
    plot([1;max(samples(:,1))],[upbound;upbound])
    xlabel('Sample')
    ylabel('Data')
    title('Input Data')

    subplot(2,1,2)
    plot(samples,xDespiked)
    hold on
    scatter(samples(Indx,1),xDespiked(Indx,1))
    plot([1;max(samples(:,1))],[lowbound;lowbound])
    plot([1;max(samples(:,1))],[upbound;upbound])
    xlabel('Sample')
    ylabel('Data')
    title('Despiked Data')
    
end

%--------------------------------------------------------------------------
    