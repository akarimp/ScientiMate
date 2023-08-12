function [c, yPredicted] = curvefit2d(x, y, MathExpression, coefIniGuess, fitmethod, dispout)
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

curvefit2d
==========

.. code:: MATLAB

    [c, yPredicted] = curvefit2d(x, y, MathExpression, coefIniGuess, fitmethod, dispout)

Description
-----------

Fit curve to 2 dimensinal input dataset

Inputs
------

x
    x data
y
    y data
MathExpression
    | Right hand side of the relation y=f(x) as 'f(x)'
    | Example: y=c(1).*x.^2+c(2).*x+c(3) then MathExpression='c(1).*x.^2+c(2).*x+c(3)'
    | Example: y=c(1).*exp(x)+c(2) then MathExpression='c(1).*exp(x)+c(2)'
    | Example: y=c(1).*sin(x)+c(2) then MathExpression='c(1).*sin(x)+c(2)'
    | Desired coefficients to be found should be as : c(1), c(2),...,c(n)
coefIniGuess
    | Initial guess for desired coefficients to be found
    | coefIniGuess=[guess1,guess2,...,guessn]
    | guess1 is initial guess for c(1),...,guessn is initial guess for c(n) 
fitmethod
    | Fitting method: 
    | 'lsq': curve-fitting using nonlinear least-squares  
    | 'fmin': curve-fitting by minimizing a sum of squared errors
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

c
    Desired coefficients to be found
yPredicted
    Predicted value from fitted curve

Examples
--------

.. code:: MATLAB

    x(:,1)=linspace(0,10,100);
    y(:,1)=0.5.*x.^2+2.*x+10+(-10+(10-(-10))).*rand(100,1);
    [c,yPredicted]=curvefit2d(x,y,'c(1).*x.^2+c(2).*x+c(3)',[1,1,1],'lsq','yes');

    x(:,1)=linspace(0,10,100);
    y(:,1)=5.*x.^2+7+(-200+(200-(-200))).*rand(100,1);
    [c,yPredicted]=curvefit2d(x,y,'c(1).*x.^c(2)+c(3)',[1,2,10],'lsq','yes');

    x(:,1)=linspace(0,10,100);
    y(:,1)=0.5.*exp(x)+100+(-200+(200-(-200))).*rand(100,1);
    [c,yPredicted]=curvefit2d(x,y,'c(1).*exp(x)+c(2)',[1,1],'fmin','yes');

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
    case 4
        fitmethod='fmin'; dispout='no';
    case 5
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
%Fitting curve 

if strcmp(fitmethod,'lsq')==1
    
    %Constructing mathematical expression for defining function
    MathExpressionforFun=strcat('@(c,x)','(',MathExpression,')'); 
    
    %Define Anonymous Functions
    %fun=@(c,x)(c(1).*x.^2+c(2).*x+c(3))
    fun = str2func(MathExpressionforFun);

    %Calculating coefficient for input function
    c = lsqcurvefit(fun,coefIniGuess,x,y);

elseif strcmp(fitmethod,'fmin')==1 

    %Constructing mathematical expression for defining function
    MathExpressionforFun=strcat('@(c)','(sum((',MathExpression,'-y).^2))'); %Sum of squared error

    %Define Anonymous Functions
    % fun = @(x)(sum((x(1)*Windvelall.^x(2)+x(3))-Turbdespikeall).^2)
    fun = str2func(MathExpressionforFun);

    %Calculating coefficient for input function
    c = fminsearch(fun,coefIniGuess);

end

%--------------------------------------------------------------------------
%Calculating Predicted value from fitted curve

MathExpressionforyPredicted=strcat('@(c,x)','(',MathExpression,')');
fun = str2func(MathExpressionforyPredicted);
yPredicted=fun(c,x);

%--------------------------------------------------------------------------
%Calculating statistical properties for y and yPredicted data

ymean=mean(y(:,1)); %Mean value of y data
ystd=std(y(:,1)); %Standard deviation of y data

yPredictedmean=mean(yPredicted(:,1)); %Mean value of yPredicted data
yPredictedstd=std(yPredicted(:,1)); %Standard deviation of yPredicted data

%--------------------------------------------------------------------------
%Calculating goodness of fit parameters (quantitavie evaluation of model) 

%Pearson's correlation coefficient
r=(sum((y-ymean).*(yPredicted-yPredictedmean)))/(sqrt((sum((y-ymean).^2))*(sum((yPredicted-yPredictedmean).^2)))); %Pearson's correlation coefficient 

%Coefficient of determination
%R2=1-(sum((y-yPredicted).^2))/(sum((y-ymean).^2)); %Coefficient of determination
R2=r^2; %Coefficient of determination

%Root-mean-square error
RMSE=sqrt((sum((yPredicted-y).^2))/N); %Root-mean-square error

%Mean absolute error 
MAE=(sum(abs(yPredicted-y)))/N; %Mean absolute error

%Scatter index
SI=RMSE/((sum(y(:,1)))/N); %Scatter index

%Nash Sutcliffe efficiency coefficient
NSE=1-((sum((y-yPredicted).^2))/(sum((y-ymean).^2))); %Nash Sutcliffe efficiency coefficient

%Index of agreement 
d=1-(sum((y-yPredicted).^2))/(sum((abs(yPredicted-ymean)+abs(y-ymean)).^2)); %Index of agreement

%Bias
%Bias=((sum(yPredicted(:,1)))/N-(sum(y(:,1)))/N); %Bias
Bias=(yPredictedmean-ymean); %Bias

%Normalized mean bias
%NMBias=((sum(yPredicted(:,1)))/N-(sum(y(:,1)))/N)/((sum(y(:,1)))/N); %Normalized mean bias
NMBias=(yPredictedmean-ymean)/(ymean); %Normalized mean bias

%Relative error
RE=(yPredicted-y)./y; %Relative error

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1
    
    for i=1:length(c)
    val=[c(i)];
    name={strcat('c',num2str(i))};
        fprintf('%14s   %g\n',name{1},val(1));
    end

    val=[RMSE R2 NTotal N NNaNInf];
    name={'RMSE','R2','N Input','N Used','N NaN Inf'};
    for i=1:length(val)
        fprintf('%14s   %g\n',name{i},val(i));
    end
    
    %Plotting
    scatter(x,y)
    hold on
    plot(x,yPredicted)
    xlabel('x')
    ylabel('y')
    legend('Data','Fitted Curve')

end

%--------------------------------------------------------------------------
