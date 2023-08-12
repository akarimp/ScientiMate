function [xReplaced, NaN_Indx] = replacemissing2d(x, what2replace, gridsize_x, gridsize_y, interpMethod, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2020-02-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

replacemissing2d
================

.. code:: MATLAB

    [xReplaced, NaN_Indx] = replacemissing2d(x, what2replace, gridsize_x, gridsize_y, interpMethod, dispout)

Description
-----------

Replace missing data points in 2d array

Inputs
------

x
    Input data
what2replace='both';
    | What needs to be replaced
    | 'NaN': replacing NaN data points
    | 'Inf': replacing Inf data points
    | 'both': replacing NaN and Inf data points
    | Number: replacing data points equal to Number
gridsize_x=1;
    | Grid size (distance between grid points) in x direction
    | Leave gridsize_x=1 if you do not have it
gridsize_y=1;
    | Grid size (distance between grid points) in y direction
    | Leave gridsize_y=1 if you do not have it
interpMethod='nearest';
    | Interpolation method
    | 'linear': Use default or 'linear' method to interpolate
    | 'nearest': Use nearest neighbor method to interpolate
    | 'knn': Use nearest neighbor method to interpolate (Use 'knn' for large array)
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

xReplaced
    Replaced data
NaN_Indx
    Logical index of replaced points

Examples
--------

.. code:: MATLAB

    x=[[1,0,3];[2,5,NaN];[3,NaN,1];[5,7,2]];
    [xReplaced, NaN_Indx] = replacemissing2d(x, 'NaN', 1, 1, 'nearest', 'yes');

    xgrid=randn(100,50);
    xgrid(randi(100,20,1),randi(50,20,1))=NaN;
    [xReplaced, NaN_Indx] = replacemissing2d(xgrid, 'NaN', 1, 1, 'knn', 'yes');

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
        what2replace='both'; gridsize_x=1; gridsize_y=1; interpMethod='nearest'; dispout='no';
    case 2
        gridsize_x=1; gridsize_y=1; interpMethod='nearest'; dispout='no';
    case 3
        gridsize_y=1; interpMethod='nearest'; dispout='no';
    case 4
        interpMethod='nearest'; dispout='no';
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
%Replacing missing data points

%Assigning empty to NaN_Indx
NaN_Indx=[];

%Defining missing data points
if strcmp(what2replace,'NaN')==1
    %Define NaN data points
    NaN_Indx=(isnan(x)==1);
    Valid_Indx=(isnan(x)~=1);
elseif strcmp(what2replace,'Inf')==1
    %Define Inf data points
    NaN_Indx=(isinf(x)==1);
    Valid_Indx=(isinf(x)~=1);
elseif strcmp(what2replace,'both')==1
    %Define NaN and Inf data points
    NaN_Indx=(isnan(x)==1 | isinf(x)==1);
    Valid_Indx=(isnan(x)~=1 | isinf(x)~=1);
elseif isnumeric(what2replace)==1
    %Define NaN and Inf data points
    NaN_Indx=(x==what2replace);
    Valid_Indx=(x~=what2replace);
end

%Assigning x to xReplaced
xReplaced=x;

%Create 2d i and j positions
[M,N]=size(x);
[samples_ii,samples_jj]=meshgrid(1:N,1:M);
samples_ii=gridsize_x.*samples_ii; %Map position to have correct grid length
samples_jj=gridsize_y.*samples_jj; %Map position to have correct grid length

%Replacing missing points
if strcmp(interpMethod,'knn')==1

    if sum(NaN_Indx(:))>0

        [NaN_row,NaN_col]=find(NaN_Indx==1); %Find row and columns that have NaN values
    
        %for i=1:M
            %for j=1:N
        for i=1:length(NaN_row(:,1))
            Indx_ii=NaN_row(i,1); %Index of NaN point in x direction
            Indx_jj=NaN_col(i,1); %Index of NaN point in y direction
            if NaN_Indx(Indx_ii,Indx_jj)==1
                Dist=sqrt((samples_ii-samples_ii(Indx_ii,Indx_jj)).^2+(samples_jj-samples_jj(Indx_ii,Indx_jj)).^2); %Calculating distance of (x,y) to (xPoint,yPoint)
                Dist(NaN_Indx==1)=max(Dist(:)); %make sure a NaN point is not replace by another NaN point (NaN point distance to itself is zero)
                [minDist,MinDist_Indx]=min(Dist(:)); %Finding minimum distances
                x_1d=xReplaced(:);
                xReplaced(Indx_ii,Indx_jj)=x_1d(MinDist_Indx,1); %Value of the nearest point to (Indx_ii,Indx_jj)
            end
        end

    end

else

    if sum(NaN_Indx(:))>0
        xReplaced(NaN_Indx==1)=griddata(samples_ii(Valid_Indx==1),samples_jj(Valid_Indx==1),xReplaced(Valid_Indx==1),samples_ii(NaN_Indx==1),samples_jj(NaN_Indx==1),interpMethod);

        %Replacing NaN data point resulted from default method with ones from nearest method
        if sum(isnan(xReplaced(:)))>0
            [samples_ii,samples_jj,xReplacednearest]=griddata(samples_ii,samples_jj,xReplaced,samples_ii,samples_jj,'nearest');
            xReplaced(isnan(xReplaced)==1)=xReplacednearest(isnan(xReplaced)==1); 
        end

    end
    
end

%Assigning a replaced data to x for next repeat
x=xReplaced;


%--------------------------------------------------------------------------
%Defining replaced points

if sum(NaN_Indx(:))>0
    NumberReplacedPoint=sum(NaN_Indx(:));
    PercentReplacedPoint=sum(NaN_Indx(:))/numel(x)*100;
else
    NumberReplacedPoint=0;
    PercentReplacedPoint=0;
end

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1

    disp(['Total number of pints replaced = ',num2str(NumberReplacedPoint)])
    disp(['Percent of replaced points     = ',num2str(PercentReplacedPoint),' %'])
    
    subplot(2,1,1)
    imagesc(xInput)
    set(gca,'YDir','normal'); %Correct y axis direction
    title('Input Data')

    subplot(2,1,2)
    imagesc(xReplaced)
    set(gca,'YDir','normal'); %Correct y axis direction
    hold on
    scatter(samples_ii(NaN_Indx)./gridsize_x,samples_jj(NaN_Indx)./gridsize_y,[],'r','filled')
    title('Replaced Data')
    
end

%--------------------------------------------------------------------------
