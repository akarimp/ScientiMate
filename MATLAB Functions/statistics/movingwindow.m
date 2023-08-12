function [moving_statistics] = movingwindow(x, WindowSize, StatisticalMethod, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2020-08-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

movingwindow
============

.. code:: MATLAB

    [moving_statistics] = movingwindow(x, WindowSize, StatisticalMethod, dispout)

Description
-----------

Calculate statistics of moving window through 1-d x data

Inputs
------

x
    Input data
WindowSize=3;
    | Window size (number of adjacent elements) that is used for moving window, should be equal or larger than 3
    | Window size should be an odd integer
StatisticalMethod='mean';
    | Statistical value of moving window to be reported:
    | 'mean': Moving mean
    | 'std': Moving standard deviation
    | 'min': Moving minimum
    | 'max': Moving maximum
    | 'sum': Moving sum
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

moving_statistics
    Statistical value of moving window

Examples
--------

.. code:: MATLAB

    fs=128;
    t(:,1)=linspace(0,9.5,10*fs);
    x(:,1)=sin(2*pi*0.3*t)+0.1*sin(2*pi*4*t);
    [moving_statistics]=movingwindow(x,37,'mean','yes');

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
        WindowSize=3; StatisticalMethod='mean'; dispout='no';
    case 2
        StatisticalMethod='mean'; dispout='no';
    case 3
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(x)==1
    x=x';
end

%--------------------------------------------------------------------------

N=length(x(:,1)); %Length of the input data
Nw=WindowSize; %Window size

%Checking size of the window
Nw(Nw>N)=N;
Nw(Nw<=3)=3;

%--------------------------------------------------------------------------

Nw_half=fix(Nw/2); %Half of the range width

moving_statistics=zeros(N,1);

%Calculate moving mean
if strcmp(StatisticalMethod,'mean')==1

    for i=1:N
        if i>=1 & i<=Nw_half
            moving_statistics(i,1)=mean(x(1:i+Nw_half,1));
            
        elseif i>=Nw_half+1 & i<=N-Nw_half-1
            moving_statistics(i,1)=mean(x(i-Nw_half:i+Nw_half,1));
            
        elseif i>=N-Nw_half & i<=N
            moving_statistics(i,1)=mean(x(i-Nw_half:end,1));
            
        end
    end

    %Or moving_statistics=movmean(x,Nw,'Endpoints','shrink');

%Calculate moving standard deviation
elseif strcmp(StatisticalMethod,'std')==1

    for i=1:N
        if i>=1 & i<=Nw_half
            moving_statistics(i,1)=std(x(1:i+Nw_half,1));
            
        elseif i>=Nw_half+1 & i<=N-Nw_half-1
            moving_statistics(i,1)=std(x(i-Nw_half:i+Nw_half,1));
            
        elseif i>=N-Nw_half & i<=N
            moving_statistics(i,1)=std(x(i-Nw_half:end,1));
            
        end
    end
    
    %Or moving_statistics=movstd(x,Nw,'Endpoints','shrink');

%Calculate moving minimum
elseif strcmp(StatisticalMethod,'min')==1

    for i=1:N
        if i>=1 & i<=Nw_half
            moving_statistics(i,1)=min(x(1:i+Nw_half,1));
            
        elseif i>=Nw_half+1 & i<=N-Nw_half-1
            moving_statistics(i,1)=min(x(i-Nw_half:i+Nw_half,1));
            
        elseif i>=N-Nw_half & i<=N
            moving_statistics(i,1)=min(x(i-Nw_half:end,1));
            
        end
    end
    
    %Or moving_statistics=movmin(x,Nw,'Endpoints','shrink');

%Calculate moving maximum
elseif strcmp(StatisticalMethod,'max')==1

    for i=1:N
        if i>=1 & i<=Nw_half
            moving_statistics(i,1)=max(x(1:i+Nw_half,1));
            
        elseif i>=Nw_half+1 & i<=N-Nw_half-1
            moving_statistics(i,1)=max(x(i-Nw_half:i+Nw_half,1));
            
        elseif i>=N-Nw_half & i<=N
            moving_statistics(i,1)=max(x(i-Nw_half:end,1));
            
        end
    end
    
    %Or moving_statistics=movmax(x,Nw,'Endpoints','shrink');

%Calculate moving sum
elseif strcmp(StatisticalMethod,'sum')==1

    for i=1:N
        if i>=1 & i<=Nw_half
            moving_statistics(i,1)=sum(x(1:i+Nw_half,1));
            
        elseif i>=Nw_half+1 & i<=N-Nw_half-1
            moving_statistics(i,1)=sum(x(i-Nw_half:i+Nw_half,1));
            
        elseif i>=N-Nw_half & i<=N
            moving_statistics(i,1)=sum(x(i-Nw_half:end,1));
            
        end
    end
    
    %Or moving_statistics=movsum(x,Nw,'Endpoints','shrink');

end

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1
    
    plot(x)
    hold on
    plot(moving_statistics)
    
    xlabel('Sample')
    ylabel('Data')
    legend('Input','Moving Statistics')
    
end

%--------------------------------------------------------------------------
