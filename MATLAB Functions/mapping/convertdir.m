function [dirout] = convertdir(dirin, CalcMethod)
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

convertdir
==========

.. code:: MATLAB

    [dirout] = convertdir(dirin, CalcMethod)

Description
-----------

Convert direction from one system to another one

Inputs
------

dirin
    Direction to be converted in (Degree)
CalcMethod='metetotrig';
    | Direction conversion method 
    | 'metetotrig': Converting meteorological direction to trigonometric direction, trig=-mete+270;
    | 'trigtomete': Converting trigonometric direction to meteorological direction, mete=270-trig;
    | 'aztomete': Converting azimuth (bearing) to meteorological direction, mete=az+180;
    | 'metetoaz': Converting meteorological direction to azimuth (bearing), az=mete-180;
    | 'aztotrig': Converting azimuth (bearing) to trigonometric direction, trig=-az+90;
    | 'trigtoaz': Converting trigonometric direction to azimuth (bearing), az=90-trig;
    | meteorological direction:
    |     0 (degree): from North, 90 (degree): from East, 180 (degree): from South, 270 (degree): from West 
    | azimuth (bearing) direction which is measured clockwise from the north:
    |     0 (degree): toward North, 90 (degree): toward East, 180 (degree): toward South, 270 (degree): toward West 

Outputs
-------

dirout
    Converted direction in (Degree)

Examples
--------

.. code:: MATLAB

    metedir=[0;90;180;270;360]; 
    [dirout]=convertdir(metedir,'metetotrig');

    azimuth=[0;90;180;270;360]; 
    [dirout]=convertdir(azimuth,'aztomete');

    azimuth=[0;90;180;270;360]; 
    [dirout]=convertdir(azimuth,'aztotrig');

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
        CalcMethod='metetotrig';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(dirin)==1
    dirin=dirin';
end

%--------------------------------------------------------------------------
%Converting direction 

%Converting meteorological direction to trigonometric direction
if strcmp(CalcMethod,'metetotrig')==1
    
    dirout=-dirin+270;

%Converting trigonometric direction to meteorological direction
elseif strcmp(CalcMethod,'trigtomete')==1

    dirout=270-dirin;

%Converting azimuth (bearing) to meteorological direction
elseif strcmp(CalcMethod,'aztomete')==1

    dirout=dirin+180;

%Converting meteorological direction to azimuth (bearing)
elseif strcmp(CalcMethod,'metetoaz')==1

    dirout=dirin-180;

%Converting azimuth (bearing) to trigonometric direction
elseif strcmp(CalcMethod,'aztotrig')==1

    dirout=-dirin+90;

%Converting trigonometric direction to azimuth (bearing)
elseif strcmp(CalcMethod,'trigtoaz')==1

    dirout=90-dirin;

end

%Add 360 to all numbers to have them all positive
%Use mod(360) to take care of the ones larger than 360
dirout=mod((dirout+360),360); 

%--------------------------------------------------------------------------
