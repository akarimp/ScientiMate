function [CmapAKLines, CmapAKSet1, CmapAKPaired, CmapAKExtraLines, CmapAKTab10, CBrewerSet1, CBrewerPaired, Tab10] = seqcolormap(dispout)
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

seqcolormap
===========

.. code:: MATLAB

    [CmapAKLines, CmapAKSet1, CmapAKPaired, CmapAKExtraLines, CmapAKTab10, CBrewerSet1, CBrewerPaired, Tab10] = seqcolormap(dispout)

Description
-----------

Generate sequential colormap for drawing lines

Inputs
------

dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

CmapAKLines
    7 Colors similar to matlab Lines in RGB (between 0 and 1)
CmapAKSet1
    9 Colors similar to colorbrewer2.org Set1 in RGB (between 0 and 1)
CmapAKPaired
    12 Colors similar to colorbrewer2.org Paired in RGB (between 0 and 1)
CmapAKExtraLines
    15 extra line colors similar to cbrewer qualitative Paired in RGB (between 0 and 1)
CmapAKTab10
    10 Colors similar to classic Tableau 10 color palettes in RGB (between 0 and 1)
CBrewerSet1
    9 Colors Set1 from Cynthia Brewer's ColorBrewer colorbrewer2.org in RGB (between 0 and 1)
CBrewerPaired
    12 Colors Paired from Cynthia Brewer's ColorBrewer colorbrewer2.org in RGB (between 0 and 1)
Tab10
    10 Colors from classic Tableau 10 color palettes in RGB (between 0 and 1)

Examples
--------

.. code:: MATLAB

    [CmapAKLines,CmapAKSet1,CmapAKPaired,CmapAKExtraLines,CmapAKTab10,CBrewerSet1,CBrewerPaired,Tab10]=seqcolormap('yes');

References
----------

Colormap

* http://colorbrewer2.org
* http://matplotlib.org/cmocean/
* https://matplotlib.org/users/colormaps.html
* http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml
* https://www.giss.nasa.gov/tools/panoply/colorbars/
* http://jdherman.github.io/colormap/

Source

* http://www.beautycolorcode.com  sample: http://www.beautycolorcode.com/0072bd  (RGB and hex code)
* http://www.htmlcsscolor.com/
* http://www.colorhexa.com/
* http://colorbrewer2.org/
* http://onlinehelp.tableau.com/current/pro/desktop/en-us/help.htm#formatting_create_custom_colors.html

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
        dispout='no';
end

%--------------------------------------------------------------------------
%Colormaps

%7 Colors similar to matlab Lines
CmapAKLines=[[0,114,187];[226,88,34];[225,173,33];[128,55,144];[122,172,33];[102,183,225];[159,29,53]]./255;

%9 Colors similar to to colorbrewer2.org Set1
CmapAKSet1=[[206,22,32];[61,133,184];[75,163,81];[159,95,159];[255,140,0];[253,227,54];[161,82,38];[247,127,190];[169,157,157]]./255;

%12 Colors similar to to colorbrewer2.org Paired
CmapAKPaired=[[164,210,224];[49,110,160];[165,215,133];[34,139,34];[255,152,137];[206,22,32];[255,200,120];[255,140,0];[202,180,212];[107,63,160];[255,255,153];[177,89,47]]./255;

%15 extra line colors similar to cbrewer qualitative Paired
CmapAKExtraLines=[[164,210,224];[0,153,255];[0,127,255];[0,102,204];[0,0,255];[165,215,133];[75,163,81];[255,153,51];[255,102,0];[226,88,34];[177,89,47];[224,176,255];[202,180,212];[255,153,153];[230,32,32]]./255;

%10 Colors similar to classic Tableau 10 color palettes
CmapAKTab10=[[49,110,160];[255,117,24];[34,139,34];[227,38,54];[148,112,196];[139,80,75];[247,127,190];[128,128,128];[183,198,26];[3,180,200]]./255;

%Cynthia Brewer's ColorBrewer colorbrewer2.org (should be cited if used)
CBrewerSet1=[[228,26,28];[55,126,184];[77,175,74];[152,78,163];[255,127,0];[255,255,51];[166,86,40];[247,129,191];[153,153,153]]./255;
CBrewerPaired=[[166,206,227];[31,120,180];[178,223,138];[51,160,44];[251,154,153];[227,26,28];[253,191,111];[255,127,0];[202,178,214];[106,61,154];[255,255,153];[177,89,40]]./255;

%10 Colors from classic Tableau 10 color palettes
Tab10=[[31,119,180];[255,127,14];[44,160,44];[214,39,40];[148,103,189];[140,86,75];[227,119,194];[127,127,127];[188,189,34];[23,190,207]]./255;

% CmapAKLines=[0,114,187;255,128,0;225,173,33;128,55,144;124,159,47;102,183,225;159,29,53]./255;
% CmapAKExtraLines=[164,210,224;0,153,255;0,127,255;0,102,204;0,0,255;165,215,133;75,163,81;255,153,51;255,102,0;226,88,34;177,89,47;224,176,255;202,180,212;255,153,153;230,32,32]./255;

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1

    
    figure(1)
    %subplot(3,3,1)
    subplot(2,2,1)
    N=7; % Number of color in color set
    X = linspace(0,pi*3,1000);
    Y = bsxfun(@(x,n)sin(x+2*n*pi/N), X.', 1:N);
    % axes('NextPlot','replacechildren', 'ColorOrder',CmapAKLines);
    set(gca,'NextPlot','replacechildren', 'ColorOrder',CmapAKLines);
    plot(X,Y,'linewidth',5)
    ylim([-1.1 1.1]);
    legend('show')
    title('CmapAKLines similar to matlab lines')
    box on
    
    %figure(2)
    %subplot(3,3,2)
    subplot(2,2,2)
    N=9; % Number of color in color set
    X = linspace(0,pi*3,1000);
    Y = bsxfun(@(x,n)sin(x+2*n*pi/N), X.', 1:N);
    % axes('NextPlot','replacechildren', 'ColorOrder',CmapAKSet1);
    set(gca,'NextPlot','replacechildren', 'ColorOrder',CmapAKSet1);
    plot(X,Y,'linewidth',5)
    ylim([-1.1 1.1]);
    legend('show')
    title('CmapAKSet1 similar to colorbrewer Set1')
    box on

    %figure(3)
    %subplot(3,3,3)
    subplot(2,2,3)
    N=12; % Number of color in color set
    X = linspace(0,pi*3,1000);
    Y = bsxfun(@(x,n)sin(x+2*n*pi/N), X.', 1:N);
    % axes('NextPlot','replacechildren', 'ColorOrder',CmapAKPaired);
    set(gca,'NextPlot','replacechildren', 'ColorOrder',CmapAKPaired);
    plot(X,Y,'linewidth',5)
    ylim([-1.1 1.1]);
    legend('show')
    title('CmapAKPaired similar to colorbrewer Paired')
    box on

    %figure(4)
    %subplot(3,3,4)
    subplot(2,2,4)
    N=15; % Number of color in color set
    X = linspace(0,pi*3,1000);
    Y = bsxfun(@(x,n)sin(x+2*n*pi/N), X.', 1:N);
    % axes('NextPlot','replacechildren', 'ColorOrder',CmapAKExtraLines);
    set(gca,'NextPlot','replacechildren', 'ColorOrder',CmapAKExtraLines);
    plot(X,Y,'linewidth',5)
    ylim([-1.1 1.1]);
    legend('show')
    title('CmapAKExtraLines')
    box on

    %figure(5)
    figure(2)
    %subplot(3,3,5)
    subplot(2,2,1)
    N=10; % Number of color in color set
    X = linspace(0,pi*3,1000);
    Y = bsxfun(@(x,n)sin(x+2*n*pi/N), X.', 1:N);
    % axes('NextPlot','replacechildren', 'ColorOrder',CmapAKExtraLines);
    set(gca,'NextPlot','replacechildren', 'ColorOrder',CmapAKTab10);
    plot(X,Y,'linewidth',5)
    ylim([-1.1 1.1]);
    legend('show')
    title('CmapAKTab10 similar to classic Tableau 10 color palettes')
    box on

    %figure(6)
    %subplot(3,3,6)
    subplot(2,2,2)
    N=9; % Number of color in color set
    X = linspace(0,pi*3,1000);
    Y = bsxfun(@(x,n)sin(x+2*n*pi/N), X.', 1:N);
    % axes('NextPlot','replacechildren', 'ColorOrder',CBrewerSet1);
    set(gca,'NextPlot','replacechildren', 'ColorOrder',CBrewerSet1);
    plot(X,Y,'linewidth',5)
    ylim([-1.1 1.1]);
    legend('show')
    title('Set1 from colorbrewer2.org')
    box on

    %figure(7)
    %subplot(3,3,7)
    subplot(2,2,3)
    N=12; % Number of color in color set
    X = linspace(0,pi*3,1000);
    Y = bsxfun(@(x,n)sin(x+2*n*pi/N), X.', 1:N);
    % axes('NextPlot','replacechildren', 'ColorOrder',CBrewerPaired);
    set(gca,'NextPlot','replacechildren', 'ColorOrder',CBrewerPaired);
    plot(X,Y,'linewidth',5)
    ylim([-1.1 1.1]);
    legend('show')
    title('Paired from colorbrewer2.org')
    box on
    
    %figure(8)
    %subplot(3,3,8)
    subplot(2,2,4)
    N=10; % Number of color in color set
    X = linspace(0,pi*3,1000);
    Y = bsxfun(@(x,n)sin(x+2*n*pi/N), X.', 1:N);
    % axes('NextPlot','replacechildren', 'ColorOrder',CmapAKExtraLines);
    set(gca,'NextPlot','replacechildren', 'ColorOrder',Tab10);
    plot(X,Y,'linewidth',5)
    ylim([-1.1 1.1]);
    legend('show')
    title('Tab10 from classic Tableau 10 color palettes')
    box on

end

%--------------------------------------------------------------------------
%Color Names

%X = linspace(0,pi*3,1000);
%Y = bsxfun(@(x,n)sin(x+2*n*pi/N), X.', 1:N);

%Colors similar to Matlab Lines
% plot(X,Y,'Color',[0,114,187]./255,'linewidth',5) % French Blue #0072bb
% plot(X,Y,'Color',[226,88,34]./255,'linewidth',5) % Flame #e25822
% plot(X,Y,'Color',[225,173,33]./255,'linewidth',5) % Urobilin #e1ad21
% plot(X,Y,'Color',[128,55,144]./255,'linewidth',5) % Vivid Violet #803790
% plot(X,Y,'Color',[122,172,33]./255,'linewidth',5) % Lima #7aac21
% plot(X,Y,'Color',[102,183,225]./255,'linewidth',5) % Malibu #66b7e1
% plot(X,Y,'Color',[159,29,53]./255,'linewidth',5) % Vivid Burgundy #9f1d35

%Colors similar to colorbrewer2.org Set1
% plot(X,Y,'Color',[206,22,32]./255,'linewidth',5) % Fire Engine Red #ce1620
% plot(X,Y,'Color',[61,133,184]./255,'linewidth',5) % Curious Blue #3d85b8
% plot(X,Y,'Color',[75,163,81]./255,'linewidth',5) % Fruit Salad #4ba351
% plot(X,Y,'Color',[159,95,159]./255,'linewidth',5) % Violet Blue #9f5f9f
% plot(X,Y,'Color',[255,140,0]./255,'linewidth',5) % DarkOrange #ff8c00
% plot(X,Y,'Color',[253,227,54]./255,'linewidth',5) % Gorse #fde336
% plot(X,Y,'Color',[161,82,38]./255,'linewidth',5) % Rich Gold #a15226
% plot(X,Y,'Color',[247,127,190]./255,'linewidth',5) % Persian Pink #f77fbe
% plot(X,Y,'Color',[169,157,157]./255,'linewidth',5) % Nobel #a99d9d

%Colors similar to colorbrewer2.org Paired
% plot(X,Y,'Color',[164,210,224]./255,'linewidth',5) % French Pass #a4d2e0
% plot(X,Y,'Color',[49,110,160]./255,'linewidth',5) % Lochmara #316ea0
% plot(X,Y,'Color',[165,215,133]./255,'linewidth',5) % Feijoa #a5d785
% plot(X,Y,'Color',[34,139,34]./255,'linewidth',5) % ForestGreen #228b22
% plot(X,Y,'Color',[255,152,137]./255,'linewidth',5) % Mona Lisa #ff9889
% plot(X,Y,'Color',[206,22,32]./255,'linewidth',5) % Fire Engine Red #ce1620
% plot(X,Y,'Color',[255,200,120]./255,'linewidth',5) % Chardonnay #ffc878
% plot(X,Y,'Color',[255,140,0]./255,'linewidth',5) % DarkOrange #ff8c00
% plot(X,Y,'Color',[202,180,212]./255,'linewidth',5) % Prelude #cab4d4
% plot(X,Y,'Color',[107,63,160]./255,'linewidth',5) % Royal Purple #6b3fa0
% plot(X,Y,'Color',[255,255,153]./255,'linewidth',5) % Canary #ffff99
% plot(X,Y,'Color',[177,89,47]./255,'linewidth',5) % Fiery Orange #b1592f

%Colors similar to classic Tableau 10 color palettes
% plot(X,Y,'Color',[31,119,180]./255,'linewidth',5) % Lochmara #316EA0
% plot(X,Y,'Color',[255,117,24]./255,'linewidth',5) % Pumpkin #FF7518
% plot(X,Y,'Color',[34,139,34]./255,'linewidth',5) % ForestGreen #228B22
% plot(X,Y,'Color',[227,38,54]./255,'linewidth',5) % Alizarin #E32636
% plot(X,Y,'Color',[148,112,196]./255,'linewidth',5) % Lilac Bush #9470C4
% plot(X,Y,'Color',[139,80,75]./255,'linewidth',5) % Lotus #8B504B
% plot(X,Y,'Color',[247,127,190]./255,'linewidth',5) % Persian Pink #F77FBE
% plot(X,Y,'Color',[128,128,128]./255,'linewidth',5) % Gray #808080
% plot(X,Y,'Color',[183,198,26]./255,'linewidth',5) % Rio Grande #B7C61A
% plot(X,Y,'Color',[3,180,200]./255,'linewidth',5) % Iris Blue #03B4C8

%Other Colors
% plot(X,Y,'Color',[164,210,224]./255,'linewidth',5) % French Pass #a4d2e0
% plot(X,Y,'Color',[0,153,255]./255,'linewidth',5) % Blue #0099ff
% plot(X,Y,'Color',[0,127,255]./255,'linewidth',5) % Azure #007fff
% plot(X,Y,'Color',[0,102,204]./255,'linewidth',5) % Blue #0066cc
% plot(X,Y,'Color',[0,0,255]./255,'linewidth',5) % Blue #0000ff
% plot(X,Y,'Color',[165,215,133]./255,'linewidth',5) % Feijoa #a5d785
% plot(X,Y,'Color',[75,163,81]./255,'linewidth',5) % Fruit Salad #4ba351
% plot(X,Y,'Color',[255,153,51]./255,'linewidth',5) % Deep Saffron #ff9933
% plot(X,Y,'Color',[255,102,0]./255,'linewidth',5) % Safety Orange #ff6600
% plot(X,Y,'Color',[226,88,34]./255,'linewidth',5) % Flame #e25822
% plot(X,Y,'Color',[177,89,47]./255,'linewidth',5) % Fiery Orange #b1592f
% plot(X,Y,'Color',[224,176,255]./255,'linewidth',5) % Mauve #e0b0ff
% plot(X,Y,'Color',[202,180,212]./255,'linewidth',5) % Prelude #cab4d4
% plot(X,Y,'Color',[255,153,153]./255,'linewidth',5) % Light Salmon Pink #a5d785
% plot(X,Y,'Color',[230,32,32]./255,'linewidth',5) % Lust #e62020

%--------------------------------------------------------------------------

% subplot(1,3,1)
% plot(X,Y,'Color',cmap(2,:),'linewidth',5) % cbrewer Paired Blue
% subplot(1,3,2)
% plot(X,Y,'Color',[0,114,189]./255,'linewidth',5) % Matlab Lines Blue #0072bd
% subplot(1,3,3)
% plot(X,Y,'Color',[0,114,187]./255,'linewidth',5) % French Blue #0072bb

% subplot(1,3,1)
% plot(X,Y,'Color',cmap(8,:),'linewidth',5) % cbrewer Paired Orange
% subplot(1,3,2)
% plot(X,Y,'Color',[217,83,25]./255,'linewidth',5) % Matlab Lines Orange #d95319
% subplot(1,3,3)
% plot(X,Y,'Color',[255,128,0]./255,'linewidth',5) % Orange #ff8000

% subplot(1,3,1)
% plot(X,Y,'Color',cmap(7,:),'linewidth',5) % cbrewer Paired Yellow
% subplot(1,3,2)
% plot(X,Y,'Color',[237,177,32]./255,'linewidth',5) % Matlab Lines Yellow #edb120
% subplot(1,3,3)
% plot(X,Y,'Color',[225,173,33]./255,'linewidth',5) % Urobilin #e1ad21

% subplot(1,3,1)
% plot(X,Y,'Color',cmap(10,:),'linewidth',5) % cbrewer Paired Purple
% subplot(1,3,2)
% plot(X,Y,'Color',[126,47,142]./255,'linewidth',5) % Matlab Lines Purple #7e2f8e
% subplot(1,3,3)
% plot(X,Y,'Color',[128,55,144]./255,'linewidth',5) % Vivid Violet #803790

% subplot(1,3,1)
% plot(X,Y,'Color',cmap(4,:),'linewidth',5) % cbrewer Paired Green
% subplot(1,3,2)
% plot(X,Y,'Color',[119,172,48]./255,'linewidth',5) % Matlab Lines Green #77ac30
% subplot(1,3,3)
% plot(X,Y,'Color',[124,159,47]./255,'linewidth',5) % Sushi #7c9f2f

% subplot(1,3,1)
% plot(X,Y,'Color',cmap(1,:),'linewidth',5) % cbrewer Paired Light Blue
% subplot(1,3,2)
% plot(X,Y,'Color',[77,190,238]./255,'linewidth',5) % Matlab Lines Light Blue #4dbeee
% subplot(1,3,3)
% plot(X,Y,'Color',[102,183,225]./255,'linewidth',5) % Malibu #66b7e1

% subplot(1,3,1)
% plot(X,Y,'Color',cmap(6,:),'linewidth',5) % cbrewer Paired Red
% subplot(1,3,2)
% plot(X,Y,'Color',[162,20,47]./255,'linewidth',5) % Matlab Lines Red #a2142f
% subplot(1,3,3)
% plot(X,Y,'Color',[159,29,53]./255,'linewidth',5) % Vivid Burgundy #9f1d35

%--------------------------------------------------------------------------
