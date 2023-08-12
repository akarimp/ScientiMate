function [indxknn, distknn] = findknn(x, y, xpoint, ypoint, numofneighbors, distCalcMethod, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-07-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

findknn
=======

.. code:: MATLAB

    [indxknn, distknn] = findknn(x, y, xpoint, ypoint, numofneighbors, distCalcMethod, dispout)

Description
-----------

Find k-nearest neighbors using Euclidean distance

Inputs
------

x
    Coordinate of data in x direction
y
    Coordinate of data in y direction
xpoint
    Coordinate (in x direction) of point that nearest point to that is disired to be found 
ypoint
    Coordinate (in y direction) of point that nearest point to that is disired to be found 
numofneighbors=1;
    Number of nearest neighbors to (xpoint,ypoint) that are desired to be found
distCalcMethod='1d';
    | Distance calculation method 
    | '1d': use 1d array
    | 'pdist2': Use 2d distance function
    | 'vector': Use vectorized distance 
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

indxknn
    | Index of nearest neighbors points
    | returns M*N array where M=length(xpoint) and N=numofneighbors
    | nth row associated with nth point in (xpoint,ypoint)
distknn
    | Distance of nearest neighbors points
    | returns M*N array where M=length(xpoint) and N=numofneighbors
    | mth row associated with mth point in (xpoint,ypoint)
    | nth column associated with nth nearest neighbors

Examples
--------

.. code:: MATLAB

    x(:,1)=10.*rand(100,1);
    y(:,1)=10.*rand(100,1);
    xpoint=mean(x);
    ypoint=mean(y);
    [indxknn,distknn]=findknn(x,y,xpoint,ypoint,1,'1d','yes');

    x(:,1)=10.*rand(100,1);
    y(:,1)=10.*rand(100,1);
    xpoint=[2.5;5;7.5];
    ypoint=[3;6;9];
    [indxknn,distknn]=findknn(x,y,xpoint,ypoint,10,'1d','yes');

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
        numofneighbors=1; distCalcMethod='1d'; dispout='no';
    case 5
        distCalcMethod='1d'; dispout='no';
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

if isrow(xpoint)==1
    xpoint=xpoint';
end

if isrow(ypoint)==1
    ypoint=ypoint';
end

%--------------------------------------------------------------------------
%Calculating distance from (x,y) to (xpoint,ypoint)

M=length(x(:,1));
N=length(xpoint(:,1));

%Calculating distance
%First method, using 1d array
if strcmp(distCalcMethod,'1d')==1

    %indxknn=ones(N,numofneighbors); %Pre-assigning array
    %distknn=ones(N,numofneighbors); %Pre-assigning array
    distxy=zeros(M,N); %Pre-assigning array
    for i=1:N
        distxy(:,i)=sqrt((x-xpoint(i,1)).^2+(y-ypoint(i,1)).^2); %Calculating distance from (x,y) to (xpoint,ypoint)
        %[distxysort,Indx]=sort(distxy,'ascend'); %Sorting distances
        %xsort=x(Indx,1); %Sorting X based on distance to (xpoint,ypoint)
        %ysort=y(Indx,1); %Sorting Y based on distance to (xpoint,ypoint)
        %indxknn(i,:)=Indx(1:numofneighbors,1); %Storing index of first numofneighbors points
        %distknn(i,:)=distxysort(1:numofneighbors,1); %Storing distance of first numofneighbors points
    end

%Second method, using pdist2 function
elseif strcmp(distCalcMethod,'pdist2')==1

    X=[x,y];
    Y=[xpoint,ypoint];
    distxy=pdist2(X,Y); %Calculating distance from (x,y) to (xpoint,ypoint)

%Third method, using vectorized distance 
elseif strcmp(distCalcMethod,'vector')==1

    X=[x,y];
    Y=[xpoint,ypoint];
    
    X2=sum(X.^2,2);
    Y2=sum(Y.^2,2);
    distxy=sqrt(ones(M,1)*Y2'+X2*ones(1,N)-2*X*Y'); %Calculating distance from (x,y) to (xpoint,ypoint)
    
end

%--------------------------------------------------------------------------
%Finding k-nearest neighbors

numofneighbors(numofneighbors<1)=1; %Minimum value for numofneighbors is one

if numofneighbors==1
    [distxymin,Indx]=min(distxy); %Finding minimum distance
    indxknn=(Indx(1,:))'; %Storing index of first numofneighbors points
    distknn=(distxymin(1,:))'; %Storing distance of first numofneighbors points
else
    [distxysort,Indx]=sort(distxy,'ascend'); %Sorting distances
    indxknn=(Indx(1:numofneighbors,:))'; %Storing index of first numofneighbors points
    distknn=(distxysort(1:numofneighbors,:))'; %Storing distance of first numofneighbors points
end

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1
    
    %Plotting data
    scatter(x,y)
    hold on
    scatter(xpoint,ypoint,'filled')
    for i=1:N
        scatter(x(indxknn(i,1:numofneighbors),1),y(indxknn(i,1:numofneighbors),1),'filled')
    end
    
    xlabel('x')
    ylabel('y')
    legend('Data','Point','KNN')
    
end

%--------------------------------------------------------------------------
