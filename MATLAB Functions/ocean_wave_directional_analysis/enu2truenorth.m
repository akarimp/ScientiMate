function [directionTN] = enu2truenorth(directionENU, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-05-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

enu2truenorth
=============

.. code:: MATLAB

    [directionTN] = enu2truenorth(directionENU, dispout)

Description
-----------

Convert mathematical direction (angle) from ENU (East North Up) coordinate system to compass direction with respect to true north

Inputs
------

directionENU
    | Direction (angle) in ENU (East North Up) coordinate system between 0 and 360 (Degree)
    | If coordinate system is ENU, then x is East and y is North  
dispout='no';
    | Define to display outputs or not ('yes': display, 'no': not display)
    | Note: inputs can be as a single value or a 1-D vertical array

Outputs
-------

directionTN
    | Direction (angle) in compass direction with respect to true north (Degree)
    | In true north coordinate system, wave comes from as:
    | 0 degree: from north, 90 degree: from east, 180 degree: from south, 270 degree: from west  

Examples
--------

.. code:: MATLAB

    [directionTN]=enu2truenorth(90,'yes');

    [directionTN]=enu2truenorth([15;30;45;60;90],'no');

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
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(directionENU)==1
    directionENU=directionENU';
end

%--------------------------------------------------------------------------
%Converting mathematical angle direction to compass direction with respect to true north

directionTN=270-directionENU; %Converting mathematical direction to compass direction
directionTN(directionTN<0)=directionTN(directionTN<0)+360;
directionTN(directionTN>360)=directionTN(directionTN>360)-360;

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1 

    %[ENUx,ENUy]=pol2cart(deg2rad(directionENU),1);
    %compass(ENUx,ENUy)

    polar([deg2rad(directionENU),deg2rad(directionENU)],[0,1])
    hold on
    polar([deg2rad(directionTN),deg2rad(directionTN)],[0,1])

    legend('ENU','True North')

end

%--------------------------------------------------------------------------
