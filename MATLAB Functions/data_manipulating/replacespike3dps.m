function [xDespiked, Indx] = replacespike3dps(x, nrpeat, DespikeScaleFactor, interpMethod, dispout)
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

replacespike3dps
================

.. code:: MATLAB

    [xDespiked, Indx] = replacespike3dps(x, nrpeat, DespikeScaleFactor, interpMethod, dispout)

Description
-----------

Remove spikes in the time series based on 3D phase space method by Goring and Nikora (2002)

Inputs
------

x
    Input data
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
    [xDespiked,Indx]=replacespike3dps(x,2,1,'linear','yes');

    fs=2;
    t(:,1)=linspace(0,1023.5,1024*fs);
    x(:,1)=detrend(0.5.*cos(2*pi*0.2*t)+(-0.1+(0.1-(-0.1))).*rand(1024*fs,1));
    spikeloc(:,1)=[10:100:length(t(:,1))];
    x(spikeloc+round(2*randn(length(spikeloc(:,1)),1)),1)=sign(randn(length(spikeloc(:,1)),1));
    [xDespiked,Indx]=replacespike3dps(x,1,1,'linear','yes');

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
        nrpeat=1; DespikeScaleFactor=1; interpMethod='linear'; dispout='no';
    case 2
        DespikeScaleFactor=1; interpMethod='linear'; dispout='no';
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

N=length(x(:,1)); %Length of the input data

%Preserving an input data
xInput=x;
xInputmean=mean(xInput(:,1)); %Mean of input data

%Removing a mean from input data
x=x-xInputmean;

%--------------------------------------------------------------------------
%Removing spikes

for i=1:nrpeat

    %Gradient
    x1stGrad=gradient(x); %First gradient of x
    x2ndGrad=gradient(x1stGrad); %Second gradient of x

    %Standard deviation
    xStd=std(x); %Standard deviation of x
    x1stGradStd=std(x1stGrad); %Standard deviation of x1stGrad
    x2ndGradStd=std(x2ndGrad); %Standard deviation of x2ndGrad

    %Rotation  angle
    %Theta=atan2(sum(x.*x2ndGrad),sum(x.^2)); %Rotation  angle  of  the  principal  axis of x2ndGrad vs x
    Theta=atan(sum(x.*x2ndGrad)./sum(x.^2)); %Rotation  angle  of  the  principal  axis of x2ndGrad vs x
     
    %Expected absolute maximum based on normal distributaion
    Lambda=sqrt(2*log(N));

    %Axis of ellipse
    a1=Lambda*xStd; %Major axis of ellipse for x1stGrad vs x
    b1=Lambda*x1stGradStd; %Minor axis of ellipse for x1stGrad vs x

    a3=Lambda*x1stGradStd; %Major axis of ellipse for x2ndGradStd vs x1stGrad
    b3=Lambda*x2ndGradStd; %Minor axis of ellipse for x2ndGradStd vs x1stGrad

    %Initial values for a2 and b2
    a2=a1; %Major axis of ellipse for x2ndGradStd vs x
    b2=b1; %Minor axis of ellipse for x2ndGradStd vs x
    a2prevstep=0;
    b2prevstep=0;

    %Loop with 5 iteration to define a2 and b2 
    for n=1:5
        a2=sqrt(((Lambda*xStd)^2-(b2*sin(Theta))^2)/((cos(Theta))^2)); %Major axis of ellipse for x2ndGradStd vs x
        b2=sqrt(((Lambda*x2ndGradStd)^2-(a2*sin(Theta))^2)/((cos(Theta))^2)); %Minor axis of ellipse for x2ndGradStd vs x
        a2diff=a2-a2prevstep;
        b2diff=b2-b2prevstep;
        a2prevstep=a2;
        b2prevstep=b2;
    end

    %Additional loop (if needed) with (100-5) iterations to define a2 and b2
    while n<=100 & (abs(a2diff)>0.0001 | abs(b2diff)>0.0001)
        a2=sqrt(((Lambda*xStd)^2-(b2*sin(Theta))^2)/((cos(Theta))^2)); %Major axis of ellipse for x2ndGradStd vs x
        b2=sqrt(((Lambda*x2ndGradStd)^2-(a2*sin(Theta))^2)/((cos(Theta))^2)); %Minor axis of ellipse for x2ndGradStd vs x
        a2diff=a2-a2prevstep;
        b2diff=b2-b2prevstep;
        a2prevstep=a2;
        b2prevstep=b2;
        n=n+1;
    end

    %Calculating mean values
    xmean=mean(x(:,1));
    x1stGradmean=mean(x1stGrad(:,1));
    x2ndGradmean=mean(x2ndGrad(:,1));

    %Calculating distances based on ellipse equation
    distellipse12=((x-xmean).^2./a1.^2+(x1stGrad-x1stGradmean).^2./b1.^2);
    distellipse13=((x-xmean).^2./a2.^2+(x2ndGrad-x2ndGradmean).^2./b2.^2);
    distellipse23=((x1stGrad-x1stGradmean).^2./a3.^2+(x2ndGrad-x2ndGradmean).^2./b3.^2);

    %Defining spike points, if a point locatated outside ellipse it is considered as spike
    Indx1=find(distellipse12>1*DespikeScaleFactor);
    Indx2=find(distellipse13>1*DespikeScaleFactor);
    Indx3=find(distellipse23>1*DespikeScaleFactor);

    %Fixing an indexing, method may pick 2 points right before and after the spike
    for m=2:length(Indx1)
        if Indx1(m,1)-Indx1(m-1,1)==2
            Indxmean=(Indx1(m,1)+Indx1(m-1,1))/2;
            Indx1(m,1)=Indxmean;
            Indx1(m-1,1)=Indxmean;
        end
    end

    for m=2:length(Indx2)
        if Indx2(m,1)-Indx2(m-1,1)==2
            Indxmean=(Indx2(m,1)+Indx2(m-1,1))/2;
            Indx2(m,1)=Indxmean;
            Indx2(m-1,1)=Indxmean;
        end
    end

    for m=2:length(Indx3)
        if Indx3(m,1)-Indx3(m-1,1)==2
            Indxmean=(Indx3(m,1)+Indx3(m-1,1))/2;
            Indx3(m,1)=Indxmean;
            Indx3(m-1,1)=Indxmean;
        end
    end

    %Setting spike points as NaN
    xDespiked=x;
    xDespiked(Indx1,1)=NaN;
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

%Adding mean back to data
xDespiked=xDespiked+xInputmean;

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
    