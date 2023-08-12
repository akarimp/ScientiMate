function [x_ds, y_ds, z_ds] = downsamplexyz(x, y, z, RetainRatio)
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

downsamplexyz
=============

.. code:: MATLAB

    function [x_ds, y_ds, z_ds] = downsamplexyz(x, y, z, RetainRatio)

Description
-----------

Downsample x, y,  and z data and retain given ratio

Inputs
------

x
    x data
y
    y data
z
    z data
RetainRatio=0.5;
    | Define percentage of data to retain, value between 0 and 1
    | Example: RetainRatio=0.8; means 80% of data are retained

Outputs
-------

x_ds
    Downsample x data
y_ds
    Downsample y data
z_ds
    Downsample z data

Examples
--------

.. code:: MATLAB

    x(:,1)=10.*rand(1000,1);
    y(:,1)=10.*rand(1000,1);
    z=x.^2+y.^2;
    [x_ds, y_ds, z_ds]=downsamplexyz(x, y, z, 0.7);

    x(:,1)=(-90-(-91)).*rand(1000,1)+(-91);
    y(:,1)=(31-(30)).*rand(1000,1)+(30);
    [xgrid,ygrid]=meshgrid(linspace(min(x),max(x),1000),linspace(min(y),max(y),500));
    zgrid=xgrid.^2+ygrid.^2;
    [x_ds, y_ds, z_ds]=downsamplexyz(xgrid, ygrid, zgrid, 0.3);

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
    case 3
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

if isrow(z)==1
    z=z';
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
        z_ds=z(dsIndx_dim1,1); %z data

    else

        %Defining index number to down sample data
        dsIndx_dim1(:,1)=fix(linspace(1,Sz1,fix(Sz1*RetainRatio)));
        dsIndx_dim2(:,1)=fix(linspace(1,Sz2,fix(Sz2*RetainRatio)));

        %Retaining down sampled data
        x_ds=x(dsIndx_dim1,dsIndx_dim2); %x data
        y_ds=y(dsIndx_dim1,dsIndx_dim2); %y data
        z_ds=z(dsIndx_dim1,dsIndx_dim2); %z data

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
    z_ds=z(dsIndx_dim1,dsIndx_dim2,dsIndx_dim3); %z data

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
    z_ds=z(dsIndx_dim1,dsIndx_dim2,dsIndx_dim3,dsIndx_dim4); % data

end

%--------------------------------------------------------------------------
