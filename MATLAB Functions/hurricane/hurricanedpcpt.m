function [dPcdt, dPcdthPahr] = hurricanedpcpt(Pc, dt, CalcMethod, dispout)
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

hurricanedpcpt
==============

.. code:: MATLAB

    [dPcdt, dPcdthPahr] = hurricanedpcpt(Pc, dt, CalcMethod, dispout)

Description
-----------

Calculate hurricane central pressure (Pc) intensity change over time (dPc/dt)

Inputs
------

Pc
    Hurricane central surface pressure in (Pa)
dt=6*3600;
    | Time interval between pressure data points in (s)
    | National Hurricane Center reports data every 6 hours 
CalcMethod='backward';
    | Calculation method 
    | 'forward': Calculate hurricane central pressure intensity change over time using forward difference method
    |     If CalcMethod='forward'; then last element is zero
    | 'backward': Calculate hurricane central pressure intensity change over time using backward difference method
    |     If CalcMethod='backward'; then first element is zero
    | 'central': Calculate hurricane central pressure intensity change over time using central difference method
    |     If CalcMethod='central'; then first and last elements are zero
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

dPcdt                           
    Hurricane central pressure (Pc) intensity change (dPc/dt) in (Pa/s)
dPcdthPahr
    Hurricane central pressure (Pc) intensity change (dPc/dt) in (hPa/hr)

Examples
--------

.. code:: MATLAB

    %Hurricane Katrina centeral pressure (Pa)
    Pc=[100800;100700;100700;100600;100300;100000;99700;99400;98800;98400;98300;98700;...
        97900;96800;95900;95000;94200;94800;94100;93000;90900;90200;90500;91300;...
        92000;92300;92800;94800;96100;97800;98500;99000;99400;99600];

    [dPcdt,dPcdthPahr]=hurricanedpcpt(Pc,6*3600,'backward','yes');

References
----------

Data

* www.nhc.noaa.gov/data/
* www.nhc.noaa.gov/data/hurdat/hurdat2-format-nencpac.pdf
* coast.noaa.gov/hurricanes
* www.aoml.noaa.gov/hrd/data_sub/re_anal.html

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
        dt=6*3600; CalcMethod='backward'; dispout='no';
    case 2
        CalcMethod='backward'; dispout='no';
    case 3
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(Pc)==1
    Pc=Pc';
end

%--------------------------------------------------------------------------
%Calculating hurricane central pressure (Pc) intensity change over time (dPc/dt)

dPcdt(:,1)=zeros(length(Pc(:,1)),1); %Pre-assigning array

%Calculating hurricane central pressure intensity change over time using forward difference method
if strcmp(CalcMethod,'forward')==1

    %Calculation range: (1:end-1,1), last element is zero
    dPcdt(1:end-1,1)=(Pc(2:end,1)-Pc(1:end-1,1))./dt; %Central pressure intensity change over time

%Calculating hurricane central pressure intensity change over time using backward difference method
elseif strcmp(CalcMethod,'backward')==1

    %Calculation range: (2:end,1), first element is zero
    dPcdt(2:end,1)=(Pc(2:end,1)-Pc(1:end-1,1))./dt; %Central pressure intensity change over time

%Calculating hurricane central pressure intensity change over time using centra difference method
elseif strcmp(CalcMethod,'central')==1

    %Calculation range: (2:end-1,1), first and last elements are zero
    dPcdt(2:end-1,1)=(Pc(3:end,1)-Pc(1:end-2,1))./(2*dt); %Central pressure intensity change over time

end

%Convert central pressure intensity change over time from (Pa/s) to (hPa/hr)
dPcdthPahr=dPcdt.*3600.*1e-2;

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1
    
    %Plotting data
    t(:,1)=[0:1:length(Pc(:,1))-1].*dt./3600;
    plot(t,dPcdt)
    
    xlabel('Time (hr)')
    ylabel('dPc/dp (Pa/s)')
    
end

%--------------------------------------------------------------------------
