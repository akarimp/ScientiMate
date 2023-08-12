function [x_ds, y_ds] = downsamplexy(x, y, RetainRatio)
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

downsamplexy
============

.. code:: MATLAB

    function [x_ds, y_ds] = downsamplexy(x, y, RetainRatio)

Description
-----------

Downsample x and y data and retain given ratio

Inputs
------

x
    x data
y
    y data
RetainRatio=0.5;
    | Define percentage of data to retain, value between 0 and 1
    | Example: RetainRatio=0.8; means 80% of data are retained

Outputs
-------

x_ds
    Downsample x data
y_ds
    Downsample y data

Examples
--------

.. code:: MATLAB

    x(:,1)=10.*rand(1000,1);
    y=x.^2;
    [x_ds, y_ds]=downsamplexy(x, y, 0.7);

    xgrid=(-90-(-91)).*rand(1000,500)+(-91);
    ygrid=xgrid.^2;
    [x_ds, y_ds]=downsamplexy(xgrid, ygrid, 0.3);

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
        RetainRatio=0.5;
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
%Down sampling the data

%Check RetainRatio
RetainRatio(RetainRatio<=0 | RetainRatio>1)=1;

if ndims(x)==2

    %Data size
    [Sz1,Sz2]=size(x);

    if Sz1==1 | Sz2==1
 
        %Defining index number to down sample data
        dsIndx_dim1(:,1)=fix(linspace(1,Sz1,fix(Sz1*RetainRatio)));

        %Retaining down sampled data
        x_ds=x(dsIndx_dim1,1); %x data
        y_ds=y(dsIndx_dim1,1); %y data

    else

        %Defining index number to down sample data
        dsIndx_dim1(:,1)=fix(linspace(1,Sz1,fix(Sz1*RetainRatio)));
        dsIndx_dim2(:,1)=fix(linspace(1,Sz2,fix(Sz2*RetainRatio)));

        %Retaining down sampled data
        x_ds=x(dsIndx_dim1,dsIndx_dim2); %x data
        y_ds=y(dsIndx_dim1,dsIndx_dim2); %y data

    end

elseif ndims(x)==3

    %Data size
    [Sz1,Sz2,Sz3]=size(x);

    %Defining index number to down sample data
    dsIndx_dim1(:,1)=fix(linspace(1,Sz1,fix(Sz1*RetainRatio)));
    dsIndx_dim2(:,1)=fix(linspace(1,Sz2,fix(Sz2*RetainRatio)));
    dsIndx_dim3(:,1)=fix(linspace(1,Sz3,fix(Sz3*RetainRatio)));

    %Retaining down sampled data
    x_ds=x(dsIndx_dim1,dsIndx_dim2,dsIndx_dim3); %x data
    y_ds=y(dsIndx_dim1,dsIndx_dim2,dsIndx_dim3); %y data

elseif ndims(x)==4

    %Data size
    [Sz1,Sz2,Sz3,Sz4]=size(x);

    %Defining index number to down sample data
    dsIndx_dim1(:,1)=fix(linspace(1,Sz1,fix(Sz1*RetainRatio)));
    dsIndx_dim2(:,1)=fix(linspace(1,Sz2,fix(Sz2*RetainRatio)));
    dsIndx_dim3(:,1)=fix(linspace(1,Sz3,fix(Sz3*RetainRatio)));
    dsIndx_dim4(:,1)=fix(linspace(1,Sz4,fix(Sz4*RetainRatio)));

    %Retaining down sampled data
    x_ds=x(dsIndx_dim1,dsIndx_dim2,dsIndx_dim3,dsIndx_dim4); %x data
    y_ds=y(dsIndx_dim1,dsIndx_dim2,dsIndx_dim3,dsIndx_dim4); %y data

end

%--------------------------------------------------------------------------
