function [c, zPredicted] = curvefit3d(x, y, z, MathExpression, coefIniGuess, fitmethod, dispout)
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

curvefit3d
==========

.. code:: MATLAB

    [c, zPredicted] = curvefit3d(x, y, z, MathExpression, coefIniGuess, fitmethod, dispout)

Description
-----------

Fit curve to 3 dimensinal input dataset

Inputs
------

x
    x data
y
    y data
z
    z data
MathExpression
    | Right hand side of the relation z=f(x,y) as 'f(x,y)'
    | Example: z=c(1).*x.^2+c(2).*y+c(3) then MathExpression='c(1).*x.^2+c(2).*y+c(3)'
    | Example: z=c(1).*exp(x)+c(2).*sin(y) then MathExpression='c(1).*exp(x)+c(2).*sin(y)'
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
zPredicted
    Predicted value from fitted curve

Examples
--------

.. code:: MATLAB

    [xx,yy]=meshgrid(linspace(0,10,100),linspace(0,10,100));
    zz=xx.^2+yy.^2+10+(-10+(10-(-10))).*rand(100,100);
    x(:,1)=xx(:);
    y(:,1)=yy(:);
    z(:,1)=zz(:);
    [c,zPredicted]=curvefit3d(x,y,z,'c(1).*x.^2+c(2).*y.^2',[1,1],'fmin','yes');

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
    case 5
        fitmethod='fmin'; dispout='no';
    case 6
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

if isrow(z)==1
    z=z';
end

%--------------------------------------------------------------------------
%Removing NaN and Inf values

%Total number of elements in dataset
NTotal=length(x(:,1));

%Finding NaN and Inf in x dataset
Indx1=find(isnan(x(:,1))~=1 & isinf(x(:,1))~=1);
x1=x(Indx1,1);
y1=y(Indx1,1);
z1=z(Indx1,1);

%Finding NaN and Inf in y dataset
Indx2=find(isnan(y1(:,1))~=1 & isinf(y1(:,1))~=1);
x2=x1(Indx2,1);
y2=y1(Indx2,1);
z2=z1(Indx2,1);

%Finding NaN and Inf in z dataset
Indx3=find(isnan(z2(:,1))~=1 & isinf(z2(:,1))~=1);
x3=x2(Indx3,1);
y3=y2(Indx3,1);
z3=z2(Indx3,1);

%Datasets without NaN and Inf
x=x3;
y=y3;
z=z3;

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
    %fun=@(c,x)(c(1).*x.^2+c(2).*y+c(3))
    fun = str2func(MathExpressionforFun);

    %Calculating coefficient for input function
    c = lsqcurvefit(fun,coefIniGuess,x,z);

elseif strcmp(fitmethod,'fmin')==1 

    %Constructing mathematical expression for defining function
    MathExpressionforFun=strcat('@(c)','(sum((',MathExpression,'-z).^2))'); %Sum of squared error

    %Define Anonymous Functions
    % fun = @(x)(sum((x(1)*Windvelall.^x(2)+x(3))-Turbdespikeall).^2)
    fun = str2func(MathExpressionforFun);

    %Calculating coefficient for input function
    c = fminsearch(fun,coefIniGuess);

end

%--------------------------------------------------------------------------
%Calculating Predicted value from fitted curve

MathExpressionforzPredicted=strcat('@(c,x,y)','(',MathExpression,')');
fun = str2func(MathExpressionforzPredicted);
zPredicted=fun(c,x,y);

%--------------------------------------------------------------------------
%Calculating statistical properties for z and zPredicted data

zmean=mean(z(:,1)); %Mean value of z data
zstd=std(z(:,1)); %Standard deviation of z data

zPredictedmean=mean(zPredicted(:,1)); %Mean value of zPredicted data
zPredictedstd=std(zPredicted(:,1)); %Standard deviation of zPredicted data

%--------------------------------------------------------------------------
%Calculating goodness of fit parameters (quantitavie evaluation of model) 

%Pearson's correlation coefficient
r=(sum((z-zmean).*(zPredicted-zPredictedmean)))/(sqrt((sum((z-zmean).^2))*(sum((zPredicted-zPredictedmean).^2)))); %Pearson's correlation coefficient 

%Coefficient of determination
%R2=1-(sum((z-zPredicted).^2))/(sum((z-zmean).^2)); %Coefficient of determination
R2=r^2; %Coefficient of determination

%Root-mean-square error
RMSE=sqrt((sum((zPredicted-z).^2))/N); %Root-mean-square error

%Mean absolute error 
MAE=(sum(abs(zPredicted-z)))/N; %Mean absolute error

%Scatter index
SI=RMSE/((sum(z(:,1)))/N); %Scatter index

%Nash Sutcliffe efficiency coefficient
NSE=1-((sum((z-zPredicted).^2))/(sum((z-zmean).^2))); %Nash Sutcliffe efficiency coefficient

%Index of agreement 
d=1-(sum((z-zPredicted).^2))/(sum((abs(zPredicted-zmean)+abs(z-zmean)).^2)); %Index of agreement

%Bias
%Bias=((sum(zPredicted(:,1)))/N-(sum(z(:,1)))/N); %Bias
Bias=(zPredictedmean-zmean); %Bias

%Normalized mean bias
%NMBias=((sum(zPredicted(:,1)))/N-(sum(z(:,1)))/N)/((sum(z(:,1)))/N); %Normalized mean bias
NMBias=(zPredictedmean-zmean)/(zmean); %Normalized mean bias

%Relative error
RE=(zPredicted-z)./z; %Relative error

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
    [xx,yy]=meshgrid(linspace(min(x(:,1)),max(x(:,1)),25),linspace(min(x(:,1)),max(x(:,1)),25));
    zPredictedforplot=fun(c,xx,yy);
    scatter3(x,y,z)
    hold on
    surf(xx,yy,zPredictedforplot)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    legend('Data','Fitted Curve')

end

%--------------------------------------------------------------------------
