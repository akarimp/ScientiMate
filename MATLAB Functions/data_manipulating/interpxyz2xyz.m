function [zPoint] = interpxyz2xyz(x, y, z, xPoint, yPoint, interpMethod, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-10-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

interpxyz2xyz
=============

.. code:: MATLAB

    [zPoint] = interpxyz2xyz(x, y, z, xPoint, yPoint, interpMethod, dispout)

Description
-----------

| Interpolate 2d scattered data on given point(s) by down-sampling the input data 
| For down-sampling, the first point of the k-nearest neighbors calculated from Euclidean distance is used
| interpxyz2xyz is more suitable for large dataset
| griddata is more efficient than interpxyz2xyz for regular size dataset

Inputs
------

x
    Coordinate of data in x direction, as an 1D array
y
    Coordinate of data in y direction, as an 1D array
z
    Value of data at (x,y) as z(x,y), as an 1D array
xPoint
    Coordinate (in x direction) of point that nearest point to that is desired to be found 
yPoint
    Coordinate (in y direction) of point that nearest point to that is desired to be found 
interpMethod='nearest';
    | Interpolation method 
    | 'linear': Use default or 'linear' method to interpolate
    | 'nearest': Use nearest neighbor method to interpolate
dispout='no';
    | Define to display outputs or not ('yes': display, 'no': not display)
    | '2d': 2 dimensional scatter plot 
    | 'surface': 3 dimensional surface plot 
    | 'no': not display 

Outputs
-------

zPoint
    Value of interpolated data at (xPoint,yPoint) as z(xPoint,yPoint) 

Examples
--------

.. code:: MATLAB

    x=10.*rand(100,1);
    y=10.*rand(100,1);
    z=100.*rand(100,1);
    xPoint=[2.5;5;7.5];
    yPoint=[3;6;9];
    [zPoint]=interpxyz2xyz(x,y,z,xPoint,yPoint,'linear','2d');

    x=10.*rand(100,1);
    y=10.*rand(100,1);
    z=100.*rand(100,1);
    [xgrid,ygrid]=meshgrid(linspace(min(x(:)),max(x(:)),100),linspace(min(y(:)),max(y(:)),100));
    [zgrid]=interpxyz2xyz(x,y,z,xgrid,ygrid,'nearest','no');

    x=10.*rand(100,1);
    y=10.*rand(100,1);
    z=y.*sin(x)-x.*cos(y);
    xPoint=10.*rand(10,1);
    yPoint=10.*rand(10,1);
    [zPoint]=interpxyz2xyz(x,y,z,xPoint,yPoint,'nearest','no');

    x=10.*rand(100,1);
    y=10.*rand(100,1);
    z=y.*sin(x)-x.*cos(y);
    [xgrid,ygrid]=meshgrid(linspace(min(x(:)),max(x(:)),100),linspace(min(y(:)),max(y(:)),100));
    [zgrid]=interpxyz2xyz(x,y,z,xgrid,ygrid,'nearest','surface');

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
        interpMethod='nearest'; dispout='no';
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

if isrow(xPoint)==1
    xPoint=xPoint';
end

if isrow(yPoint)==1, 
    yPoint=yPoint';
end

%--------------------------------------------------------------------------
%Distance calculation method
distCalcMethod='2d';
%    | Distance calculation method 
%    | '2d': use 2d array
%    | '1d': use 1d array
%    | 'pdist2': use pdist2 function
%    | 'vector': use vectorized distance 

%--------------------------------------------------------------------------
%Downsampling input data by finding the nearest neighbors

numofneighbors=1; %Number of nearest neighbors to (xPoint,yPoint) that are desired to be found

%M=length(xPoint(:,1));
%N=length(xPoint(1,:));
[M,N]=size(xPoint);

%First method, using 2d array
if strcmp(distCalcMethod,'2d')==1

    xnearest=ones(M,N); %Pre-assigning array
    ynearest=ones(M,N); %Pre-assigning array
    znearest=ones(M,N); %Pre-assigning array
    for i=1:M
        for j=1:N
            Dist(:,1)=sqrt((x-xPoint(i,j)).^2+(y-yPoint(i,j)).^2); %Calculating distance of (x,y) to (xPoint,yPoint)
            %[DistSort,Indx]=sort(Dist,'ascend'); %Sorting distances
            [minDist,Indx]=min(Dist); %Finding minimum distances
            xnearest(i,j)=x(Indx(1,1),1); %x of the nearest point to (xPoint,yPoint)
            ynearest(i,j)=y(Indx(1,1),1); %y of the nearest point to (xPoint,yPoint)
            znearest(i,j)=z(Indx(1,1),1); %Value of the nearest point to (xPoint,yPoint)
        end
    end

%Second method, using 1d array
elseif strcmp(distCalcMethod,'1d')==1

    %Reshaping grid data to 1D vector
    xPoint1d=reshape(xPoint,[numel(xPoint),1]);
    yPoint1d=reshape(yPoint,[numel(xPoint),1]);

    xnearest1d=ones(length(xPoint1d(:,1)),1); %Pre-assigning array
    ynearest1d=ones(length(xPoint1d(:,1)),1); %Pre-assigning array
    znearest1d=ones(length(xPoint1d(:,1)),1); %Pre-assigning array
    for i=1:length(xPoint1d(:,1))
        Dist(:,1)=sqrt((x-xPoint1d(i,1)).^2+(y-yPoint1d(i,1)).^2); %Calculating distance of (X,Y) to (XPoint,YPoint)
        %[DistSort,Indx]=sort(Dist,'ascend'); %Sorting distances
        [minDist,Indx]=min(Dist); %Finding minimum distances
        xnearest1d(i,1)=x(Indx(1,1),1); %x of the nearest point to (xPoint,yPoint)
        ynearest1d(i,1)=y(Indx(1,1),1); %y of the nearest point to (xPoint,yPoint)
        znearest1d(i,1)=z(Indx(1,1),1); %Value of the nearest point to (xPoint,yPoint)
    end

    xnearest=reshape(xnearest1d,[M,N]);
    ynearest=reshape(ynearest1d,[M,N]);
    znearest=reshape(znearest1d,[M,N]);

%Third method, using pdist2 function
elseif strcmp(distCalcMethod,'pdist2')==1

    xnearest=ones(M,N); %Pre-assigning array
    ynearest=ones(M,N); %Pre-assigning array
    znearest=ones(M,N); %Pre-assigning array
    for j=1:N
        X=[x,y];
        Y=[(xPoint(:,j)),(yPoint(:,j))];
        Dist=pdist2(X,Y); %Calculating distance of (x,y) to (xPoint,yPoint)
        %Dist=pdist2([x,y],[(xPoint(i,:))',(yPoint(i,:))']); %Calculating distance of (x,y) to (xPoint,yPoint)
        %[DistSort,Indx]=sort(Dist,'ascend'); %Sorting distances
        [minDist,Indx]=min(Dist); %Finding minimum distances
        xnearest(:,j)=x(Indx(1,:),1); %x of the nearest point to (xPoint,yPoint)
        ynearest(:,j)=y(Indx(1,:),1); %y of the nearest point to (xPoint,yPoint)
        znearest(:,j)=z(Indx(1,:),1); %Value of the nearest point to (xPoint,yPoint)
    end

%Fourth method, using vectorized distance 
elseif strcmp(distCalcMethod,'vector')==1

    xnearest=ones(M,N); %Pre-assigning array
    ynearest=ones(M,N); %Pre-assigning array
    znearest=ones(M,N); %Pre-assigning array
    [M1,N1]=size([x,y]);
    [M2,N2]=size([(xPoint(:,1)),(yPoint(:,1))]);
    for j=1:N
        X=[x,y];
        Y=[(xPoint(:,j)),(yPoint(:,j))];
        %[M1,N1]=size(X);
        %[M2,N2]=size(Y);
        
        X2=sum(X.^2,2);
        Y2=sum(Y.^2,2);
        Dist=sqrt(ones(M1,1)*Y2'+X2*ones(1,M2)-2*X*Y'); %Calculating distance of (x,y) to (xPoint,yPoint)
        
        %X2=sum(([x,y]).^2,2);
        %Y2=sum(([(xPoint(i,:))',(yPoint(i,:))']).^2,2);
        %Dist=sqrt(ones(M1,1)*Y2'+X2*ones(1,M2)-2*([x,y])*([(xPoint(i,:))',(yPoint(i,:))'])'); %Calculating distance of (x,y) to (xPoint,yPoint)

        %[DistSort,Indx]=sort(Dist,'ascend'); %Sorting distances
        [minDist,Indx]=min(Dist); %Finding minimum distances
        xnearest(:,j)=x(Indx(1,:),1); %x of the nearest point to (xPoint,yPoint)
        ynearest(:,j)=y(Indx(1,:),1); %y of the nearest point to (xPoint,yPoint)
        znearest(:,j)=z(Indx(1,:),1); %Value of the nearest point to (xPoint,yPoint)
    end

end

%--------------------------------------------------------------------------
%Interpolating data into (xPoint,yPoint)

%Interpolating data into grid using default or linear method
if strcmp(interpMethod,'linear')==1
    [zPoint]=griddata(xnearest,ynearest,znearest,xPoint,yPoint,interpMethod);

    %Replacing NaN data point resulted from interpolation by nearest point
    if sum(isnan(zPoint(:)))>0
        zPoint(isnan(zPoint)==1)=znearest(isnan(zPoint)==1); 
    end

%Interpolating data into grid using nearest neighbor method
elseif strcmp(interpMethod,'nearest')==1
    zPoint=znearest; 

end

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'2d')==1
    
    %Plotting data
    scatter(x,y)
    hold on
    scatter(xPoint,yPoint,'filled')
    scatter(xnearest,ynearest,'filled')
    
    xlabel('x')
    ylabel('y')
    legend('Input Data','Interpolated Data','KNN')
    
elseif strcmp(dispout,'surface')==1
    
    %Plotting data
    scatter3(x,y,z,'r','filled')
    hold on
    surf(xPoint,yPoint,zPoint)
    
    colorbar
    xlabel('x')
    ylabel('y')
    zlabel('z')
    legend('Input Data','Interpolated Data')
    
end

%--------------------------------------------------------------------------
