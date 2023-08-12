function dataoverview(x, y, z)
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

dataoverview
============

.. code:: MATLAB

    dataoverview(x, y, z)

Description
-----------

Display an overview of the input data

Inputs
------

x
    x data
y=[];
    y data (Optional)
z=[];
    z data (Optional)

Outputs
-------


Examples
--------

.. code:: MATLAB

    x(:,1)=(-0.1+(0.1-(-0.1))).*randn(1024*2,1);
    dataoverview(x);

    x(:,1)=(-0.1+(0.1-(-0.1))).*randn(1024*2,1);
    y(:,1)=x+(-0.01+(0.01-(-0.01))).*randn(1024*2,1);
    dataoverview(x,y);

    x(:,1)=(-0.1+(0.1-(-0.1))).*randn(1024*2,1);
    y(:,1)=(-0.2+(0.2-(-0.2))).*randn(1024*2,1);
    z=x.^2+y.^2+0.1+(-0.1+(0.1-(-0.1))).*rand(1024*2,1);
    dataoverview(x,y,z);

References
----------

* https://www.mathworks.com/help/matlab/ref/fprintf.html

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
        y=[]; z=[];
    case 2
        z=[];
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(x)==1
    x=x';
end

if isempty(y)~=1 & isrow(y)==1
    y=y';
end

if isempty(z)~=1 & isrow(z)==1
    z=z';
end

%--------------------------------------------------------------------------
%Removing NaN and Inf values

%Total number of elements in dataset
NTotal=length(x(:,1));

if numel(x)>0 & numel(y)==0 & numel(z)==0 %Only x is available

    %Finding NaN and Inf in x dataset
    Indx1=find(isnan(x(:,1))~=1 & isinf(x(:,1))~=1);
    x1=x(Indx1,1);

    %Datasets without NaN and Inf
    x=x1;

elseif numel(x)>0 & numel(y)>0 & numel(z)==0 %x and y are available

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

elseif numel(x)>0 & numel(y)>0 & numel(z)>0 %x, y and z are available

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

end

%Total number of elements in dataset after removal of NaN and Inf
N=length(x(:,1));

%Total number of NaN and Inf elements in dataset
NNaNInf=NTotal-N;

%--------------------------------------------------------------------------
%Calculating statistical properties for x, y and z data

if numel(x)>0 & numel(y)==0 & numel(z)==0 %Only x is available

    xmean=mean(x(:,1)); %Mean value of x data
    xmin=min(x(:,1)); %Minimum value of x data
    xmax=max(x(:,1)); %Maximum value of x data
    xstd=std(x(:,1)); %Standard deviation of x data
    xSE=xstd/sqrt(length(x(:,1))); %Standard error of x data
    xCV=xstd/xmean; %Coefficient of variation of x data

elseif numel(x)>0 & numel(y)>0 & numel(z)==0 %x and y are available

    xmean=mean(x(:,1)); %Mean value of x data
    xmin=min(x(:,1)); %Minimum value of x data
    xmax=max(x(:,1)); %Maximum value of x data
    xstd=std(x(:,1)); %Standard deviation of x data
    xSE=xstd/sqrt(length(x(:,1))); %Standard error of x data
    xCV=xstd/xmean; %Coefficient of variation of x data

    ymean=mean(y(:,1)); %Mean value of y data
    ymin=min(y(:,1)); %Minimum value of y data
    ymax=max(y(:,1)); %Maximum value of y data
    ystd=std(y(:,1)); %Standard deviation of y data
    ySE=ystd/sqrt(length(y(:,1))); %Standard error of y data
    yCV=ystd/ymean; %Coefficient of variation of y data

elseif numel(x)>0 & numel(y)>0 & numel(z)>0 %x, y and z are available

    xmean=mean(x(:,1)); %Mean value of x data
    xmin=min(x(:,1)); %Minimum value of x data
    xmax=max(x(:,1)); %Maximum value of x data
    xstd=std(x(:,1)); %Standard deviation of x data
    xSE=xstd/sqrt(length(x(:,1))); %Standard error of x data
    xCV=xstd/xmean; %Coefficient of variation of x data

    ymean=mean(y(:,1)); %Mean value of y data
    ymin=min(y(:,1)); %Minimum value of y data
    ymax=max(y(:,1)); %Maximum value of y data
    ystd=std(y(:,1)); %Standard deviation of y data
    ySE=ystd/sqrt(length(y(:,1))); %Standard error of y data
    yCV=ystd/ymean; %Coefficient of variation of y data

    zmean=mean(z(:,1)); %Mean value of z data
    zmin=min(z(:,1)); %Minimum value of z data
    zmax=max(z(:,1)); %Maximum value of z data
    zstd=std(z(:,1)); %Standard deviation of z data
    zSE=zstd/sqrt(length(z(:,1))); %Standard error of z data
    zCV=zstd/zmean; %Coefficient of variation of z data

end

%--------------------------------------------------------------------------
%Calculating goodness of fit parameters (quantitavie evaluation of model) 

if numel(x)>0 & numel(y)>0 & numel(z)==0 %x and y are available

    %x-y data
    %Pearson's correlation coefficient
    rxy=(sum((x-xmean).*(y-ymean)))/(sqrt((sum((x-xmean).^2))*(sum((y-ymean).^2)))); %Pearson's correlation coefficient 

    %Coefficient of determination
    R2xy=rxy^2; %Coefficient of determination

elseif numel(x)>0 & numel(y)>0 & numel(z)>0 %x, y and z are available

    %x-y data
    %Pearson's correlation coefficient
    rxy=(sum((x-xmean).*(y-ymean)))/(sqrt((sum((x-xmean).^2))*(sum((y-ymean).^2)))); %Pearson's correlation coefficient 

    %Coefficient of determination
    R2xy=rxy^2; %Coefficient of determination

    %x-z data
    %Pearson's correlation coefficient
    rxz=(sum((x-xmean).*(z-zmean)))/(sqrt((sum((x-xmean).^2))*(sum((z-zmean).^2)))); %Pearson's correlation coefficient 

    %Coefficient of determination
    R2xz=rxz^2; %Coefficient of determination

    %y-z data
    %Pearson's correlation coefficient
    ryz=(sum((y-ymean).*(z-zmean)))/(sqrt((sum((y-ymean).^2))*(sum((z-zmean).^2)))); %Pearson's correlation coefficient 

    %Coefficient of determination
    R2yz=ryz^2; %Coefficient of determination

end

%--------------------------------------------------------------------------
%Calculating probability density distribution

if numel(x)>0 & numel(y)==0 & numel(z)==0 %Only x is available

    %Probability of x data
    binedgex(:,1)=linspace(min(x(:,1)),max(x(:,1)),11); %Bin edges for x data
    binwidthx(:,1)=diff(binedgex(:,1)); %Bin width for x data, length(binwidth) is one element less than length(binedge) 
    bincenterx(:,1)=binedgex(1:end-1,1)+binwidthx./2; %Bin center for x data

    [NPointsx(:,1),bincenteroutx(:,1)]=hist(x,bincenterx); %Calculating number of data points in each bin for x data
    histareax=NPointsx.*binwidthx; %Area under each element of histogram for x data
    fdensityx=NPointsx./(sum(histareax)); %Calculating probability density distribution for x data
    %Note: sum(fdensityx.*binwidthx)=1

elseif numel(x)>0 & numel(y)>0 & numel(z)==0 %x and y are available

    %Probability of x data
    binedgex(:,1)=linspace(min(x(:,1)),max(x(:,1)),11); %Bin edges for x data
    binwidthx(:,1)=diff(binedgex(:,1)); %Bin width for x data, length(binwidth) is one element less than length(binedge) 
    bincenterx(:,1)=binedgex(1:end-1,1)+binwidthx./2; %Bin center for x data

    [NPointsx(:,1),bincenteroutx(:,1)]=hist(x,bincenterx); %Calculating number of data points in each bin for x data
    histareax=NPointsx.*binwidthx; %Area under each element of histogram for x data
    fdensityx=NPointsx./(sum(histareax)); %Calculating probability density distribution for x data
    %Note: sum(fdensityx.*binwidthx)=1

    %Probability of y data
    binedgey(:,1)=linspace(min(y(:,1)),max(y(:,1)),11); %Bin edges for y data
    binwidthy(:,1)=diff(binedgey(:,1)); %Bin width for y data, length(binwidth) is one element less than length(binedge) 
    bincentery(:,1)=binedgey(1:end-1,1)+binwidthy./2; %Bin center for y data

    [NPointsy(:,1),bincenterouty(:,1)]=hist(y,bincentery); %Calculating number of data points in each bin for y data
    histareay=NPointsy.*binwidthy; %Area under each element of histogram for y data
    fdensityy=NPointsy./(sum(histareay)); %Calculating probability density distribution for y data
    %Note: sum(fdensityy.*binwidthy)=1

elseif numel(x)>0 & numel(y)>0 & numel(z)>0 %x, y and z are available

    %Probability of x data
    binedgex(:,1)=linspace(min(x(:,1)),max(x(:,1)),11); %Bin edges for x data
    binwidthx(:,1)=diff(binedgex(:,1)); %Bin width for x data, length(binwidth) is one element less than length(binedge) 
    bincenterx(:,1)=binedgex(1:end-1,1)+binwidthx./2; %Bin center for x data

    [NPointsx(:,1),bincenteroutx(:,1)]=hist(x,bincenterx); %Calculating number of data points in each bin for x data
    histareax=NPointsx.*binwidthx; %Area under each element of histogram for x data
    fdensityx=NPointsx./(sum(histareax)); %Calculating probability density distribution for x data
    %Note: sum(fdensityx.*binwidthx)=1

    %Probability of y data
    binedgey(:,1)=linspace(min(y(:,1)),max(y(:,1)),11); %Bin edges for y data
    binwidthy(:,1)=diff(binedgey(:,1)); %Bin width for y data, length(binwidth) is one element less than length(binedge) 
    bincentery(:,1)=binedgey(1:end-1,1)+binwidthy./2; %Bin center for y data

    [NPointsy(:,1),bincenterouty(:,1)]=hist(y,bincentery); %Calculating number of data points in each bin for y data
    histareay=NPointsy.*binwidthy; %Area under each element of histogram for y data
    fdensityy=NPointsy./(sum(histareay)); %Calculating probability density distribution for y data
    %Note: sum(fdensityy.*binwidthy)=1

    %Probability of z data
    binedgez(:,1)=linspace(min(z(:,1)),max(z(:,1)),11); %Bin edges for z data
    binwidthz(:,1)=diff(binedgez(:,1)); %Bin width for z data, length(binwidth) is one element less than length(binedge) 
    bincenterz(:,1)=binedgez(1:end-1,1)+binwidthz./2; %Bin center for z data

    [NPointsz(:,1),bincenteroutz(:,1)]=hist(z,bincenterz); %Calculating number of data points in each bin for z data
    histareaz=NPointsz.*binwidthz; %Area under each element of histogram for z data
    fdensityz=NPointsz./(sum(histareaz)); %Calculating probability density distribution for z data
    %Note: sum(fdensityz.*binwidthz)=1

end

%--------------------------------------------------------------------------
%Displaying results

if numel(x)>0 & numel(y)==0 & numel(z)==0 %Only x is available

    disp('--------------------------------------------------')
    val=[NTotal N NNaNInf];
    name={'N Input','N Used','N NaN Inf'};
    for i=1:length(val)
        fprintf('     %-15s     %-g\n',name{i},val(i));
    end

    disp('--------------------------------------------------')
    val=[xmean xmin xmax xstd xSE xCV];
    name={'x mean','x min','x max','x st. dev.','x st. err.','x coef. of var.'};
    for i=1:length(val)
        fprintf('     %-15s     %-+g\n',name{i},val(i));
    end

    disp('--------------------------------------------------')

    %Plotting
    subplot(2,1,1)
    plot(x)

    xlabel('Data Points')
    ylabel('x Data')

    subplot(2,3,5)
    %bar(bincenterx,fdensityx)
    plot(bincenterx,fdensityx)
    hold on
    plot([xmean;xmean],[min(fdensityx);max(fdensityx)])
    plot([xmean+xstd;xmean+xstd],[min(fdensityx);max(fdensityx)])
    plot([xmean-xstd;xmean-xstd],[min(fdensityx);max(fdensityx)])
    title('Probability Density Distribution, x Data')
    ylabel('Probability Density, f')
    xlabel('x Data')
    %legend('x pdf','x mean','mean+std','mean-std')

elseif numel(x)>0 & numel(y)>0 & numel(z)==0 %x and y are available

    disp('--------------------------------------------------')
    val=[NTotal N NNaNInf];
    name={'N Input','N Used','N NaN Inf'};
    for i=1:length(val)
        fprintf('     %-15s     %-g\n',name{i},val(i));
    end

    disp('--------------------------------------------------')
    val=[xmean xmin xmax xstd xSE xCV];
    name={'x mean','x min','x max','x st. dev.','x st. err.','x coef. of var.'};
    for i=1:length(val)
        fprintf('     %-15s     %-+g\n',name{i},val(i));
    end

    disp('--------------------------------------------------')
    val=[ymean ymin ymax ystd ySE yCV];
    name={'y mean','y min','y max','y st. dev.','y st. err.','y coef. of var.'};
    for i=1:length(val)
        fprintf('     %-15s     %-+g\n',name{i},val(i));
    end

    disp('--------------------------------------------------')
    val=[R2xy];
    name={'R2 x-y'};
    for i=1:length(val)
        fprintf('     %-15s     %-+g\n',name{i},val(i));
    end

    disp('--------------------------------------------------')

    %Plotting
    subplot(1,3,1)
    scatter(x,y)

    xlabel('x Data')
    ylabel('y Data')

    subplot(1,3,2)
    %bar(bincenterx,fdensityx)
    plot(bincenterx,fdensityx)
    hold on
    plot([xmean;xmean],[min(fdensityx);max(fdensityx)])
    plot([xmean+xstd;xmean+xstd],[min(fdensityx);max(fdensityx)])
    plot([xmean-xstd;xmean-xstd],[min(fdensityx);max(fdensityx)])
    title('Probability Density Distribution, x Data')
    xlabel('x Data')
    ylabel('Probability Density, f')
    %legend('x pdf','x mean','mean+std','mean-std')

    subplot(1,3,3)
    %bar(bincentery,fdensityy)
    plot(bincentery,fdensityy)
    hold on
    plot([ymean;ymean],[min(fdensityy);max(fdensityy)])
    plot([ymean+ystd;ymean+ystd],[min(fdensityy);max(fdensityy)])
    plot([ymean-ystd;ymean-ystd],[min(fdensityy);max(fdensityy)])
    title('Probability Density Distribution, y Data')
    xlabel('y Data')
    ylabel('Probability Density, f')
    %legend('y pdf','y mean','mean+std','mean-std')

elseif numel(x)>0 & numel(y)>0 & numel(z)>0 %x, y and z are available

    disp('--------------------------------------------------')
    val=[NTotal N NNaNInf];
    name={'N Input','N Used','N NaN Inf'};
    for i=1:length(val)
        fprintf('     %-15s     %-g\n',name{i},val(i));
    end

    disp('--------------------------------------------------')
    val=[xmean xmin xmax xstd xSE xCV];
    name={'x mean','x min','x max','x st. dev.','x st. err.','x coef. of var.'};
    for i=1:length(val)
        fprintf('     %-15s     %-+g\n',name{i},val(i));
    end

    disp('--------------------------------------------------')
    val=[ymean ymin ymax ystd ySE yCV];
    name={'y mean','y min','y max','y st. dev.','y st. err.','y coef. of var.'};
    for i=1:length(val)
        fprintf('     %-15s     %-+g\n',name{i},val(i));
    end

    disp('--------------------------------------------------')
    val=[zmean zmin zmax zstd zSE zCV];
    name={'z mean','z min','z max','z st. dev.','z st. err.','z coef. of var.'};
    for i=1:length(val)
        fprintf('     %-15s     %-+g\n',name{i},val(i));
    end

    disp('--------------------------------------------------')
    val=[R2xy R2xz R2yz];
    name={'R2 x-y','R2 x-z','R2 y-z'};
    for i=1:length(val)
        fprintf('     %-15s     %-+g\n',name{i},val(i));
    end

    disp('--------------------------------------------------')

    %Plotting
    subplot(2,3,1)
    scatter(x,z)

    xlabel('x Data')
    ylabel('z Data')

    subplot(2,3,2)
    scatter(y,z)

    xlabel('y Data')
    ylabel('z Data')

    subplot(2,3,3)
    scatter3(x,y,z)

    xlabel('x Data')
    ylabel('y Data')
    zlabel('z Data')

    subplot(2,3,4)
    %bar(bincenterx,fdensityx)
    plot(bincenterx,fdensityx)
    hold on
    plot([xmean;xmean],[min(fdensityx);max(fdensityx)])
    plot([xmean+xstd;xmean+xstd],[min(fdensityx);max(fdensityx)])
    plot([xmean-xstd;xmean-xstd],[min(fdensityx);max(fdensityx)])
    title('Probability Density Distribution, x Data')
    xlabel('x Data')
    ylabel('Probability Density, f')
    %legend('x pdf','x mean','mean+std','mean-std')

    subplot(2,3,5)
    %bar(bincentery,fdensityy)
    plot(bincentery,fdensityy)
    hold on
    plot([ymean;ymean],[min(fdensityy);max(fdensityy)])
    plot([ymean+ystd;ymean+ystd],[min(fdensityy);max(fdensityy)])
    plot([ymean-ystd;ymean-ystd],[min(fdensityy);max(fdensityy)])
    title('Probability Density Distribution, y Data')
    xlabel('y Data')
    ylabel('Probability Density, f')
    %legend('y pdf','y mean','mean+std','mean-std')

    subplot(2,3,6)
    %bar(bincenterz,fdensityz)
    plot(bincenterz,fdensityz)
    hold on
    plot([zmean;zmean],[min(fdensityz);max(fdensityz)])
    plot([zmean+zstd;zmean+zstd],[min(fdensityz);max(fdensityz)])
    plot([zmean-zstd;zmean-zstd],[min(fdensityz);max(fdensityz)])
    title('Probability Density Distribution, z Data')
    xlabel('z Data')
    ylabel('Probability Density, f')
    %legend('z pdf','z mean','mean+std','mean-std')

end

%--------------------------------------------------------------------------
