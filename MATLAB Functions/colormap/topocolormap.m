function [cmap_ncolor, cmap_water, cmap_land] = topocolormap(zcolormap, ncolor, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2022-05-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

topocolormap
============

.. code:: MATLAB

    [cmap_ncolor, cmap_water, cmap_land] = topocolormap(zcolormap, ncolor, dispout)

Description
-----------

Export a topographic colormap

Inputs
------

zcolormap='topocmap';
    | Colormap for z data
    | Topographic (Water/Land) colormaps:
    | 'topocmap': colormap developed by Arash Karimpour
    | 'topocmaprelief': colormap developed by Arash Karimpour
    | 'topocmapocean': colormap developed by Arash Karimpour
    | 'topocmapearth': colormap developed by Arash Karimpour
    | 'blueocean': colormap developed by Arash Karimpour
    | 'blueoceansea': colormap developed by Arash Karimpour
    | 'greenearth': colormap developed by Arash Karimpour
    | 'greenearthland': colormap developed by Arash Karimpour
    | 'blgrtopocmap': colormap developed by Arash Karimpour
    | 'blrdtopocmap': colormap developed by Arash Karimpour
    | 'grayearth': colormap developed by Arash Karimpour
    | 'etopo1': ETOPO1 colormap, https://www.ngdc.noaa.gov/mgg/global/global.html
    | 'gmtglobe': GMT_globe colormap, https://www.giss.nasa.gov/tools/panoply/colorbars/
    | 'gmtrelief': GMT_relief colormap, https://www.giss.nasa.gov/tools/panoply/colorbars/
    | 'aendekerk': Colormap from  Florian Aendekerk, http://www.mathworks.com/matlabcentral/fileexchange/63590-landseacolormap-m-
    | Any other available color map such as 'cool', 'winter', etc can be used
    | Colormap can be defined by user as [n*3] array in RGB color format between 0 and 255
ncolor=256;
    Number of colors to be used in colormap
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

cmap_ncolor
    | Colormap for z levels with ncolor number of colors in RGB color format between 0 and 1
    | To convert 0-1 scale to 0-255 scale, multiply cmap_ncolor values by 255
cmap_water
    Colormap for water in RGB color format between 0 and 1
cmap_land
    Colormap for land in RGB color format between 0 and 1

Examples
--------

.. code:: MATLAB

    [cmap_ncolor,cmap_water,cmap_land]=topocolormap('topocmap',256,'yes');
    [cmap_ncolor,cmap_water,cmap_land]=topocolormap('cool',256,'yes');

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
    zcolormap='topocmap'; ncolor=256; dispout='no';
    case 1
    ncolor=256; dispout='no';
    case 2
    dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

%--------------------------------------------------------------------------
%Colormaps

%'topocmap' colormap
%Colors from: http://htmlcolorcodes.com
if strcmp(zcolormap,'topocmap')==1

    %Water colormap
    topocmap_water=[[21, 67, 96];[27, 79, 114];[40, 116, 166];[52, 152, 219];[133, 193, 233];[214, 234, 248]]./255;

    %Land colormap
    topocmapgreen=[[56, 142, 60]]./255; %Green: 1/30 land colormap
    topocmapyellow=[[252, 243, 207];[249, 231, 159];[247, 220, 111];[241, 196, 15]]./255; %Yellow: 3/30 land colormap
    topocmapbrown=[[243, 156, 18];[230, 126, 34];[211, 84, 0];[186, 74, 0];[160, 64, 0];[135, 54, 0];[110, 44, 0]]./255;
    topocmapwhite=[[236, 240, 241]]./255; %White: 1/30 land colormap

    %Interpolating brwon color from 7 colors to 24 colors
    topocmapbrown1(:,1:3)=interp1(linspace(1,24,length(topocmapbrown(:,1))),topocmapbrown(:,1:3),linspace(1,24,24)); %Brown: 24/30 land colormap

    %Assigning water and land colormap
    cmap_water1=topocmap_water;
    cmap_land1=[topocmapgreen;topocmapyellow;topocmapbrown1;topocmapwhite];

    %Convert to 0.55% sea and 0.45% land 
    %Mariana Trench depth= 10994 m
    %Mount Everest height= 8848 m
    %(max depth)/(max depth+max height)~=0.55
    cmap_water(:,1:3)=interp1(linspace(1,55,length(cmap_water1(:,1))),cmap_water1(:,1:3),linspace(1,55,55));
    cmap_land(:,1:3)=interp1(linspace(1,45,length(cmap_land1(:,1))),cmap_land1(:,1:3),linspace(1,45,45));

    %Assigning negative (water) and positive (land) colormap
    cmap_z=[cmap_water;cmap_land];

%'topocmaprelief' colormap
%Colors from globalrelief.m
elseif strcmp(zcolormap,'topocmaprelief')==1

    %Water colormap
    topocmap_water=bone(128);
    topocmap_water=topocmap_water(65:128,:);
    
    %Land colormap
    topocmapsummer=summer(10);
    topocmapcopper=flipud(copper(54));
    
    %Assigning water and land colormap
    cmap_water1=topocmap_water;
    cmap_land1=[topocmapsummer;topocmapcopper];

    %Convert to 0.55% sea and 0.45% land 
    %Mariana Trench depth= 10994 m
    %Mount Everest height= 8848 m
    %(max depth)/(max depth+max height)~=0.55
    cmap_water(:,1:3)=interp1(linspace(1,55,length(cmap_water1(:,1))),cmap_water1(:,1:3),linspace(1,55,55));
    cmap_land(:,1:3)=interp1(linspace(1,45,length(cmap_land1(:,1))),cmap_land1(:,1:3),linspace(1,45,45));

    %Assigning negative (water) and positive (land) colormap
    cmap_z=[cmap_water;cmap_land];

%'topocmapocean' colormap
%Colors from: http://htmlcolorcodes.com
elseif strcmp(zcolormap,'topocmapocean')==1

    %Water colormap
    topocmap_water=[[21, 67, 96];[27, 79, 114];[40, 116, 166];[52, 152, 219];[133, 193, 233];[214, 234, 248]]./255;

    cmap_z=topocmap_water;

    %Assigning half of the colormap to negative (water) values and other half to positive (land) values 
    [M1,N1]=size(cmap_z);
    cmap_water=cmap_z(1:fix(M1/2),:);
    cmap_land=cmap_z(fix(M1/2)+1:end,:);

%'topocmapearth' colormap
%Colors from: http://htmlcolorcodes.com
elseif strcmp(zcolormap,'topocmapearth')==1

    %Land colormap
    topocmapgreen=[[56, 142, 60]]./255; %Green: 1/30 land colormap
    topocmapyellow=[[252, 243, 207];[249, 231, 159];[247, 220, 111];[241, 196, 15]]./255; %Yellow: 3/30 land colormap
    topocmapbrown=[[243, 156, 18];[230, 126, 34];[211, 84, 0];[186, 74, 0];[160, 64, 0];[135, 54, 0];[110, 44, 0]]./255;
    topocmapwhite=[[236, 240, 241]]./255; %White: 1/30 land colormap

    %Interpolating brwon color from 7 colors to 24 colors
    topocmapbrown1(:,1:3)=interp1(linspace(1,24,length(topocmapbrown(:,1))),topocmapbrown(:,1:3),linspace(1,24,24)); %Brown: 24/30 land colormap

    %Assigning water and land colormap
    cmap_land1=[topocmapgreen;topocmapyellow;topocmapbrown1;topocmapwhite];

    cmap_z=cmap_land1;

    %Assigning half of the colormap to negative (water) values and other half to positive (land) values 
    [M1,N1]=size(cmap_z);
    cmap_water=cmap_z(1:fix(M1/2),:);
    cmap_land=cmap_z(fix(M1/2)+1:end,:);

%'blueocean' colormap
%Colors from: http://htmlcolorcodes.com
elseif strcmp(zcolormap,'blueocean')==1

    %Water colormap
    blueoceanwater=[[52, 152, 219];[133, 193, 233]]./255;

    %Land colormap
    topocmapgreen=[[56, 142, 60]]./255; %Green: 1/30 land colormap
    topocmapyellow=[[252, 243, 207];[249, 231, 159];[247, 220, 111];[241, 196, 15]]./255; %Yellow: 3/30 land colormap
    topocmapbrown=[[243, 156, 18];[230, 126, 34];[211, 84, 0];[186, 74, 0];[160, 64, 0];[135, 54, 0];[110, 44, 0]]./255;
    topocmapwhite=[[236, 240, 241]]./255; %White: 1/30 land colormap

    %Interpolating brwon color from 7 colors to 24 colors
    topocmapbrown1(:,1:3)=interp1(linspace(1,24,length(topocmapbrown(:,1))),topocmapbrown(:,1:3),linspace(1,24,24)); %Brown: 24/30 land colormap

    %Assigning water and land colormap
    cmap_water1=blueoceanwater;
    cmap_land1=[topocmapgreen;topocmapyellow;topocmapbrown1;topocmapwhite];

    %Convert to 0.55% sea and 0.45% land 
    %Mariana Trench depth= 10994 m
    %Mount Everest height= 8848 m
    %(max depth)/(max depth+max height)~=0.55
    cmap_water(:,1:3)=interp1(linspace(1,55,length(cmap_water1(:,1))),cmap_water1(:,1:3),linspace(1,55,55));
    cmap_land(:,1:3)=interp1(linspace(1,45,length(cmap_land1(:,1))),cmap_land1(:,1:3),linspace(1,45,45));

    %Assigning negative (water) and positive (land) colormap
    cmap_z=[cmap_water;cmap_land];

%'blueoceansea' colormap
%Colors from: http://htmlcolorcodes.com
elseif strcmp(zcolormap,'blueoceansea')==1

    %Water colormap
    blueoceanwater=[[52, 152, 219];[133, 193, 233]]./255;

    %At least 4 colors requred in colormap to guarantee at least two colors for each of cmap_water and cmap_land colormaps
    blueoceanwater4colors(:,1:3)=interp1(linspace(1,4,length(blueoceanwater(:,1))),blueoceanwater(:,1:3),linspace(1,4,4));

    cmap_z=blueoceanwater4colors;

    %Assigning half of the colormap to negative (water) values and other half to positive (land) values 
    [M1,N1]=size(cmap_z);
    cmap_water=cmap_z(1:fix(M1/2),:);
    cmap_land=cmap_z(fix(M1/2)+1:end,:);

%'greenearth' colormap
%Colors from: http://htmlcolorcodes.com
elseif strcmp(zcolormap,'greenearth')==1

    %Mariana Trench depth= 10994 m
    %Mount Everest height= 8848 m
    %(max depth)/(max depth+max height)~=0.55

    %Water colormap
    topocmap_water=[[21, 67, 96];[27, 79, 114];[40, 116, 166];[52, 152, 219];[133, 193, 233];[214, 234, 248]]./255;

    %Land colormap
    greenearthgreen=[[67, 160, 71];[255, 241, 118]]./255; %Green

    %Assigning water and land colormap
    cmap_water1=topocmap_water;
    cmap_land1=greenearthgreen;

    %Convert to 0.55% sea and 0.45% land 
    %Mariana Trench depth= 10994 m
    %Mount Everest height= 8848 m
    %(max depth)/(max depth+max height)~=0.55
    cmap_water(:,1:3)=interp1(linspace(1,55,length(cmap_water1(:,1))),cmap_water1(:,1:3),linspace(1,55,55));
    cmap_land(:,1:3)=interp1(linspace(1,45,length(cmap_land1(:,1))),cmap_land1(:,1:3),linspace(1,45,45));

    %Assigning negative (water) and positive (land) colormap
    cmap_z=[cmap_water;cmap_land];

%'greenearthland' colormap
%Colors from: http://htmlcolorcodes.com
elseif strcmp(zcolormap,'greenearthland')==1

    %Mariana Trench depth= 10994 m
    %Mount Everest height= 8848 m
    %(max depth)/(max depth+max height)~=0.55

    %Land colormap
    greenearthgreen=[[67, 160, 71];[255, 241, 118]]./255; %Green

    %At least 4 colors requred in colormap to guarantee at least two colors for each of cmap_water and cmap_land colormaps
    greenearthgreen4colors(:,1:3)=interp1(linspace(1,4,length(greenearthgreen(:,1))),greenearthgreen(:,1:3),linspace(1,4,4));

    cmap_z=greenearthgreen4colors;

    %Assigning half of the colormap to negative (water) values and other half to positive (land) values 
    [M1,N1]=size(cmap_z);
    cmap_water=cmap_z(1:fix(M1/2),:);
    cmap_land=cmap_z(fix(M1/2)+1:end,:);

%'blgrtopocmap' colormap
%Colors from: http://htmlcolorcodes.com
elseif strcmp(zcolormap,'blgrtopocmap')==1

    %Mariana Trench depth= 10994 m
    %Mount Everest height= 8848 m
    %(max depth)/(max depth+max height)~=0.55

    %Water colormap
    blgrtopocmap_water=[[52, 152, 219];[133, 193, 233]]./255;

    %Land colormap
    blgrtopocmapgreen=[[67, 160, 71];[255, 241, 118]]./255; %Green

    %Assigning water and land colormap
    cmap_water1=blgrtopocmap_water;
    cmap_land1=blgrtopocmapgreen;

    %Convert to 0.55% sea and 0.45% land 
    %Mariana Trench depth= 10994 m
    %Mount Everest height= 8848 m
    %(max depth)/(max depth+max height)~=0.55
    cmap_water(:,1:3)=interp1(linspace(1,55,length(cmap_water1(:,1))),cmap_water1(:,1:3),linspace(1,55,55));
    cmap_land(:,1:3)=interp1(linspace(1,45,length(cmap_land1(:,1))),cmap_land1(:,1:3),linspace(1,45,45));

    %Assigning negative (water) and positive (land) colormap
    cmap_z=[cmap_water;cmap_land];

%'blrdtopocmap' colormap
%Colors from: http://htmlcolorcodes.com
elseif strcmp(zcolormap,'blrdtopocmap')==1

    %Mariana Trench depth= 10994 m
    %Mount Everest height= 8848 m
    %(max depth)/(max depth+max height)~=0.55

    %Water colormap
    blrdtopocmap_water=[[41, 128, 185];[212, 230, 241]]./255; %Blue

    %Land colormap
    blrdtopocmap_land=[[242, 215, 213];[192, 57, 43]]./255; %Red

    %Assigning water and land colormap
    cmap_water1=blrdtopocmap_water;
    cmap_land1=blrdtopocmap_land;

    %Convert to 0.55% sea and 0.45% land 
    %Mariana Trench depth= 10994 m
    %Mount Everest height= 8848 m
    %(max depth)/(max depth+max height)~=0.55
    cmap_water(:,1:3)=interp1(linspace(1,55,length(cmap_water1(:,1))),cmap_water1(:,1:3),linspace(1,55,55));
    cmap_land(:,1:3)=interp1(linspace(1,45,length(cmap_land1(:,1))),cmap_land1(:,1:3),linspace(1,45,45));

    %Assigning negative (water) and positive (land) colormap
    cmap_z=[cmap_water;cmap_land];

%'grayearth' colormap
%Colors from: http://htmlcolorcodes.com
elseif strcmp(zcolormap,'grayearth')==1

    %Mariana Trench depth= 10994 m
    %Mount Everest height= 8848 m
    %(max depth)/(max depth+max height)~=0.55

    %Water colormap
    topocmap_water=[[21, 67, 96];[27, 79, 114];[40, 116, 166];[52, 152, 219];[133, 193, 233];[214, 234, 248]]./255;

    %Land colormap
    grayearthgray=[[179, 182, 183];[66, 73, 73]]./255; %Gray

    %Assigning water and land colormap
    cmap_water1=topocmap_water;
    cmap_land1=grayearthgray;

    %Convert to 0.55% sea and 0.45% land 
    %Mariana Trench depth= 10994 m
    %Mount Everest height= 8848 m
    %(max depth)/(max depth+max height)~=0.55
    cmap_water(:,1:3)=interp1(linspace(1,55,length(cmap_water1(:,1))),cmap_water1(:,1:3),linspace(1,55,55));
    cmap_land(:,1:3)=interp1(linspace(1,45,length(cmap_land1(:,1))),cmap_land1(:,1:3),linspace(1,45,45));

    %Assigning negative (water) and positive (land) colormap
    cmap_z=[cmap_water;cmap_land];

%ETOPO1 colormap
%https://www.ngdc.noaa.gov/mgg/global/global.html
elseif strcmp(zcolormap,'etopo1')==1

    %http://svn.idldev.com/vis/trunk/src/color/cpt-city/ngdc/ETOPO1.cpt
    % Colormap used in the ETOPO1 global relief map:
    % http://ngdc.noaa.gov/mgg/global/global.html
    %
    % Above sea level is a modified version of GMT_globe.cpt, 
    % Designed by Lester M. Anderson (CASP, UK) lester.anderson@casp.cam.ac.uk,
    % Modified by Jesse Varner and Elliot Lim (NOAA/NGDC) with a smaller band of white for the highest elevations.
    % The ocean colors are adapted from GMT_haxby.cpt, popularized by Bill Haxby, LDEO
    % COLOR_MODEL = RGB

    %ETOPO1 elavation range
    %min elevation= -10898.0
    %max elevation= 8271.0

    etopo1=[[10,0,121];[26,0,137];[38,0,152];[27,3,166];[16,6,180];[5,9,193];[0,14,203];[0,22,210];[0,30,216];[0,39,223];...
        [12,68,231];[26,102,240];[19,117,244];[14,133,249];[21,158,252];[30,178,255];[43,186,255];[55,193,255];[65,200,255];[79,210,255];...
        [94,223,255];[138,227,255];[188,230,255];[51,102,0];[51,204,102];[187,228,146];[255,220,185];[243,202,137];[230,184,88];[217,166,39];[168,154,31];...
        [164,144,25];[162,134,19];[159,123,13];[156,113,7];[153,102,0];[162,89,89];[178,118,118];[183,147,147];[194,176,176];[204,204,204];...
        [229,229,229];[255,255,255]]./255;

    etopo1water=[[10,0,121];[26,0,137];[38,0,152];[27,3,166];[16,6,180];[5,9,193];[0,14,203];[0,22,210];[0,30,216];[0,39,223];...
        [12,68,231];[26,102,240];[19,117,244];[14,133,249];[21,158,252];[30,178,255];[43,186,255];[55,193,255];[65,200,255];[79,210,255];...
        [94,223,255];[138,227,255];[188,230,255]]./255;

    etopo1land=[[51,102,0];[51,204,102];[187,228,146];[255,220,185];[243,202,137];[230,184,88];[217,166,39];[168,154,31];...
        [164,144,25];[162,134,19];[159,123,13];[156,113,7];[153,102,0];[162,89,89];[178,118,118];[183,147,147];[194,176,176];[204,204,204];...
        [229,229,229];[255,255,255]]./255;

    %Assigning negative (water) and positive (land) colormap
    cmap_water=etopo1water;
    cmap_land=etopo1land;
    cmap_z=etopo1;

%GMT_globe colormap
%https://www.giss.nasa.gov/tools/panoply/colorbars/
elseif strcmp(zcolormap,'gmtglobe')==1

    %https://www.giss.nasa.gov/tools/panoply/colorbars/
    %https://www.giss.nasa.gov/tools/panoply/colorbars/gmt/GMT_globe.cpt
    % $Id: GMT_globe.cpt,v 1.1.1.1 2000/12/28 01:23:45 gmt Exp $
    %
    % Colormap using in global relief maps
    % Bathymetry colours manually redefined for blue-shade effect and
    % new topography colour scheme for use with Generic Mapping Tools.
    % Designed by Lester M. Anderson (CASP, UK) lester.anderson@casp.cam.ac.uk
    % COLOR_MODEL = RGB

    gmtglobe=[[153,0,255];[153,0,255];[153,0,255];[136,17,255];[119,34,255];[102,51,255];[85,68,255];[68,85,255];[51,102,255];[34,119,255];...
        [17,136,255];[0,153,255];[27,164,255];[54,175,255];[81,186,255];[108,197,255];[134,208,255];[161,219,255];[188,230,255];[215,241,255];...
        [241,252,255];[51,102,0];[51,204,102];[187,228,146];[255,220,185];[243,202,137];[230,184,88];[217,166,39];[168,154,31];[164,144,25];...
        [162,134,19];[159,123,13];[156,113,7];[153,102,0];[162,89,89];[178,118,118];[183,147,147];[194,176,176];[204,204,204];[229,229,229];...
        [242,242,242];[255,255,255]]./255;

    gmtglobewater=[[153,0,255];[153,0,255];[153,0,255];[136,17,255];[119,34,255];[102,51,255];[85,68,255];[68,85,255];[51,102,255];[34,119,255];...
        [17,136,255];[0,153,255];[27,164,255];[54,175,255];[81,186,255];[108,197,255];[134,208,255];[161,219,255];[188,230,255];[215,241,255];...
        [241,252,255]]./255;

    gmtglobeland=[[51,102,0];[51,204,102];[187,228,146];[255,220,185];[243,202,137];[230,184,88];[217,166,39];[168,154,31];[164,144,25];...
        [162,134,19];[159,123,13];[156,113,7];[153,102,0];[162,89,89];[178,118,118];[183,147,147];[194,176,176];[204,204,204];[229,229,229];...
        [242,242,242];[255,255,255]]./255;

    %Assigning negative (water) and positive (land) colormap
    cmap_water=gmtglobewater;
    cmap_land=gmtglobeland;
    cmap_z=gmtglobe;

%GMT_relief colormap
%https://www.giss.nasa.gov/tools/panoply/colorbars/
elseif strcmp(zcolormap,'gmtrelief')==1

    %https://www.giss.nasa.gov/tools/panoply/colorbars/
    %https://www.giss.nasa.gov/tools/panoply/colorbars/gmt/GMT_relief.cpt
    %	$Id: GMT_relief.cpt,v 1.1.1.1 2000/12/28 01:23:45 gmt Exp $
    %
    % Colortable for whole earth relief used in Wessel topomaps
    % Designed by P. Wessel and F. Martinez, SOEST
    % COLOR_MODEL = RGB

    gmtrelief=[[0,0,0];[0,5,25];[0,10,50];[0,80,125];[0,150,200];[86,197,184];[172,245,168];[211,250,211];[250,255,255];[70,120,50];[120,100,50];[146,126,60];...
        [198,178,80];[250,230,100];[250,234,126];[252,238,152];[252,243,177];[253,249,216];[255,255,255]]./255;

    gmtreliefwater=[[0,0,0];[0,5,25];[0,10,50];[0,80,125];[0,150,200];[86,197,184];[172,245,168];[211,250,211];[250,255,255]]./255;

    gmtreliefland=[[70,120,50];[120,100,50];[146,126,60];[198,178,80];[250,230,100];[250,234,126];[252,238,152];[252,243,177];[253,249,216];[255,255,255]]./255;

    %Assigning negative (water) and positive (land) colormap
    cmap_water=gmtreliefwater;
    cmap_land=gmtreliefland;
    cmap_z=gmtrelief;

%Colormap from  Florian Aendekerk
%http://www.mathworks.com/matlabcentral/fileexchange/63590-landseacolormap-m-
elseif strcmp(zcolormap,'aendekerk')==1

    aendekerkcmap=[[0.00000, 0.00000, 0.40000];[0.02524, 0.03810, 0.42857];[0.05048, 0.07619, 0.45714];[0.07571, 0.11429, 0.48571];[0.10095, 0.15238, 0.51429];...
        [0.12619, 0.19048, 0.54286];[0.15143, 0.22857, 0.57143];[0.17667, 0.26667, 0.60000];[0.20190, 0.30476, 0.62857];[0.22714, 0.34286, 0.65714];...
        [0.25238, 0.38095, 0.68571];[0.27762, 0.41905, 0.71429];[0.30286, 0.45714, 0.74286];[0.32810, 0.49524, 0.77143];[0.35333, 0.53333, 0.80000];...
        [0.37857, 0.57143, 0.82857];[0.40381, 0.60952, 0.85714];[0.42905, 0.64762, 0.88571];[0.45429, 0.68571, 0.91429];[0.47952, 0.72381, 0.94286];...
        [0.50476, 0.76190, 0.97143];[0.53000, 0.80000, 1.00000];[0.00000, 0.50000, 0.00000];[1.00000, 1.00000, 0.80000];[1.00000, 0.86667, 0.53333];...
        [1.00000, 0.73333, 0.26667];[1.00000, 0.60000, 0.00000];[1.00000, 0.60000, 0.00000];[0.95000, 0.56250, 0.00000];[0.90000, 0.52500, 0.00000];...
        [0.85000, 0.48750, 0.00000];[0.80000, 0.45000, 0.00000];[0.75000, 0.41250, 0.00000];[0.70000, 0.37500, 0.00000];[0.65000, 0.33750, 0.00000];...
        [0.60000, 0.30000, 0.00000];[0.55000, 0.26250, 0.00000];[0.50000, 0.22500, 0.00000];[0.45000, 0.18750, 0.00000];[0.40000, 0.15000, 0.00000]];

    aendekerkcmap_water=[[0.00000, 0.00000, 0.40000];[0.02524, 0.03810, 0.42857];[0.05048, 0.07619, 0.45714];[0.07571, 0.11429, 0.48571];[0.10095, 0.15238, 0.51429];...
        [0.12619, 0.19048, 0.54286];[0.15143, 0.22857, 0.57143];[0.17667, 0.26667, 0.60000];[0.20190, 0.30476, 0.62857];[0.22714, 0.34286, 0.65714];...
        [0.25238, 0.38095, 0.68571];[0.27762, 0.41905, 0.71429];[0.30286, 0.45714, 0.74286];[0.32810, 0.49524, 0.77143];[0.35333, 0.53333, 0.80000];...
        [0.37857, 0.57143, 0.82857];[0.40381, 0.60952, 0.85714];[0.42905, 0.64762, 0.88571];[0.45429, 0.68571, 0.91429];[0.47952, 0.72381, 0.94286];...
        [0.50476, 0.76190, 0.97143];[0.53000, 0.80000, 1.00000]];

    aendekerkcmap_land=[[0.00000, 0.50000, 0.00000];[1.00000, 1.00000, 0.80000];[1.00000, 0.86667, 0.53333];...
        [1.00000, 0.73333, 0.26667];[1.00000, 0.60000, 0.00000];[1.00000, 0.60000, 0.00000];[0.95000, 0.56250, 0.00000];[0.90000, 0.52500, 0.00000];...
        [0.85000, 0.48750, 0.00000];[0.80000, 0.45000, 0.00000];[0.75000, 0.41250, 0.00000];[0.70000, 0.37500, 0.00000];[0.65000, 0.33750, 0.00000];...
        [0.60000, 0.30000, 0.00000];[0.55000, 0.26250, 0.00000];[0.50000, 0.22500, 0.00000];[0.45000, 0.18750, 0.00000];[0.40000, 0.15000, 0.00000]];

    %Assigning negative (water) and positive (land) colormap
    cmap_water=aendekerkcmap_water;
    cmap_land=aendekerkcmap_land;
    cmap_z=aendekerkcmap;

%User input colormap
elseif isnumeric(zcolormap)==1

    cmap_z=zcolormap./255; %Converting RGB to values between 0 and 1
    cmap_z(cmap_z<0 )=0;
    cmap_z(cmap_z>1 )=1;

    %Assigning half of the colormap to negative (water) values and other half to positive (land) values 
    [M1,N1]=size(cmap_z);
    cmap_water=cmap_z(1:fix(M1/2),:);
    cmap_land=cmap_z(fix(M1/2)+1:end,:);

%Pre-defined colormap such as 'cool', 'winter', etc
else

    cmap_z=colormap(zcolormap);

    %Assigning half of the colormap to negative (water) values and other half to positive (land) values 
    [M1,N1]=size(cmap_z);
    cmap_water=cmap_z(1:fix(M1/2),:);
    cmap_land=cmap_z(fix(M1/2)+1:end,:);

end

%--------------------------------------------------------------------------
%Interpolating a colormap to include ncolor colors

cmap_ncolor=zeros(ncolor,3);
cmap_ncolor(:,1:3)=interp1(linspace(1,ncolor,length(cmap_z(:,1))),cmap_z(:,1:3),linspace(1,ncolor,ncolor));

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1
    
    [xgrid,ygrid]=meshgrid(linspace(1,ncolor,ncolor),linspace(1,ncolor,ncolor));
    zgrid=ygrid;

    imagesc(xgrid,ygrid,zgrid)
    set(gca,'YDir','normal'); %It is not required if flipud is used

    %Setting colormap
    colormap(cmap_ncolor)

    %Plotting colorbar
    colorbar

end

%--------------------------------------------------------------------------
