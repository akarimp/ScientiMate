function [r, R2, RMSE, MAE, SI, NSE, d, Bias, NMBias, RE] = fitgoodness(x, y, dispout)
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

fitgoodness
===========

.. code:: MATLAB

    [r, R2, RMSE, MAE, SI, NSE, d, Bias, NMBias, RE] = fitgoodness(x, y, dispout)

Description
-----------

Calculate goodness of fit parameters

Inputs
------

x
    Dataset with true (exact or expected) values, such as theoretical values
y
    | Dataset that needs to be evaluated, such as model results or estimated values
    | Accuracy of y dataset is evaluated against x dataset
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

r
    Pearson correlation coefficient
R2
    Coefficient of determination
RMSE
    Root mean square error
MAE
    Mean absolute error
SI
    Scatter index
NSE
    Nash Sutcliffe efficiency coefficient
d
    Index of agreement
Bias
    Bias
NMBias
    Normalized mean bias
RE
    Relative error

Examples
--------

.. code:: MATLAB

    x(:,1)=(-0.1+(0.1-(-0.1))).*randn(1024*2,1);
    y(:,1)=x+(-0.01+(0.01-(-0.01))).*randn(1024*2,1);
    [r,R2,RMSE,MAE,SI,NSE,d,Bias,NMBias,RE]=fitgoodness(x,y,'yes');

    x=[1;2;3;4;5;6;7;8;9;10];
    y=[1.1;1.98;3.3;4.2;4.8;5.95;7.5;7.7;8.99;10.5];
    [r,R2,RMSE,MAE,SI,NSE,d,Bias,NMBias,RE]=fitgoodness(x,y,'yes');

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
%Calculating statistical properties for x and y data

xmean=mean(x(:,1)); %Mean value of x data
xstd=std(x(:,1)); %Standard deviation of x data

ymean=mean(y(:,1)); %Mean value of y data
ystd=std(y(:,1)); %Standard deviation of y data

%--------------------------------------------------------------------------
%Calculating goodness of fit parameters (quantitavie evaluation of model) 

%Pearson's correlation coefficient
r=(sum((x-xmean).*(y-ymean)))/(sqrt((sum((x-xmean).^2))*(sum((y-ymean).^2)))); %Pearson's correlation coefficient 

%Coefficient of determination
%R2=1-(sum((x-y).^2))/(sum((x-xmean).^2)); %Coefficient of determination
R2=r^2; %Coefficient of determination

%Root-mean-square error
RMSE=sqrt((sum((y-x).^2))/N); %Root-mean-square error

%Mean absolute error 
MAE=(sum(abs(y-x)))/N; %Mean absolute error

%Scatter index
SI=RMSE/((sum(x(:,1)))/N); %Scatter index

%Nash Sutcliffe efficiency coefficient
NSE=1-((sum((x-y).^2))/(sum((x-xmean).^2))); %Nash Sutcliffe efficiency coefficient

%Index of agreement 
d=1-(sum((x-y).^2))/(sum((abs(y-xmean)+abs(x-xmean)).^2)); %Index of agreement

%Bias
%Bias=((sum(y(:,1)))/N-(sum(x(:,1)))/N); %Bias
Bias=(ymean-xmean); %Bias

%Normalized mean bias
%NMBias=((sum(y(:,1)))/N-(sum(x(:,1)))/N)/((sum(x(:,1)))/N); %Normalized mean bias
NMBias=(ymean-xmean)/(xmean); %Normalized mean bias

%Relative error
RE=(y-x)./x; %Relative error

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1

    val=[R2 r RMSE MAE SI NSE d Bias NMBias xmean xstd ymean ystd NTotal N NNaNInf];
    name={'R2','r','RMSE','MAE','SI','NSE','d','Bias','NMBias','x mean','x st. dev.','y mean','y st. dev.','N Input','N Used','N NaN Inf'};
    for i=1:length(val)
        fprintf('%14s   %g\n',name{i},val(i));
    end

    %Plotting
    scatter(x,y)
    hold on
    plot([min(x(:,1));max(x(:,1))],[min(x(:,1));max(x(:,1))])

    xlabel('x Data')
    ylabel('y Data')
    legend('Data','1:1 Line')

end

%--------------------------------------------------------------------------
