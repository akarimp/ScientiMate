function [cmap_ncolor] = gencolormap(cmapcolors, ncolor, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-08-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

gencolormap
===========

.. code:: MATLAB

    [cmap_ncolor] = gencolormap(cmapcolors, ncolor, dispout)

Description
-----------

Generate a colormap from input colors

Inputs
------

cmapcolors=[[255,255,255];[0,0,0]];
    | Colors that are used to generate colormap
    | Colors should be defined as [M*3] array in RGB color format
    | At least two colors should be defined, i.e. M should be equal or larger than 2
    | All values should be between 0 and 255
    | Any available colormap name such as 'cool', 'winter', etc also can be used
ncolor=256;
    Number of colors in generated colormap
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

cmap_ncolor
    | Colormap with ncolor number of colors in RGB color format between 0 and 1
    | To convert 0-1 scale to 0-255 scale, multiply cmap_ncolor values by 255

Examples
--------

.. code:: MATLAB

    %Blue-Red
    [cmap_ncolor]=gencolormap([[41, 128, 185];[192, 57, 43]],256,'yes');

    %Blue-White-Red
    [cmap_ncolor]=gencolormap([[41, 128, 185];[255,255,255];[192, 57, 43]],256,'yes');

    %Blue-White-Green
    [cmap_ncolor]=gencolormap([[33, 150, 243];[255,255,255];[76, 175, 80]],256,'yes');

    %Red-Brown sequential
    cmapcolors=[[192,57,43];[155,89,182];[41,128,185];[39,174,96];[241,196,15];[211,84,0]];
    [cmap_ncolor]=gencolormap(cmapcolors,10,'yes');

    %cool colormap
    [cmap_ncolor]=gencolormap('cool',256,'yes');

References
----------

Colormap

* http://colorbrewer2.org
* http://matplotlib.org/cmocean/
* https://matplotlib.org/users/colormaps.html
* http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml
* https://www.giss.nasa.gov/tools/panoply/colorbars/
* http://jdherman.github.io/colormap/

Color

* http://htmlcolorcodes.com

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
    case 0
    cmapcolors=[[255,255,255];[0,0,0]]; ncolor=256; dispout='no';
    case 1
    ncolor=256; dispout='no';
    case 2
    dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

%--------------------------------------------------------------------------
%Colormaps

%User input colormap
if isnumeric(cmapcolors)==1

    zclmap=cmapcolors./255; %Converting RGB to values between 0 and 1
    zclmap(zclmap<0 )=0;
    zclmap(zclmap>1 )=1;

%Pre-defined colormap such as 'cool', 'winter', etc
else

    zclmap=colormap(cmapcolors);

end

%--------------------------------------------------------------------------
%Interpolating a colormap to include ncolor colors

cmap_ncolor=zeros(ncolor,3);
cmap_ncolor(:,1:3)=interp1(linspace(1,ncolor,length(zclmap(:,1))),zclmap(:,1:3),linspace(1,ncolor,ncolor));

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1
    
    [xgrid,ygrid]=meshgrid(linspace(1,ncolor,ncolor),linspace(1,ncolor,ncolor));
    zgrid=ygrid;

    imagesc(xgrid(1,:),ygrid(:,1),zgrid)
    set(gca,'YDir','normal'); 

    %Setting colormap
    colormap(cmap_ncolor)

    %Plotting colorbar
    colorbar

end

%--------------------------------------------------------------------------
