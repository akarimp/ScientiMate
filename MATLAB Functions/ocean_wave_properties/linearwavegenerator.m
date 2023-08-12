function [Eta, t, Etaij] = linearwavegenerator(amin, amax, Tmin, Tmax, Phimin, Phimax, fs, duration, NoOfWave, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-01-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

linearwavegenerator
===================

.. code:: MATLAB

    [Eta, t, Etaij] = linearwavegenerator(amin, amax, Tmin, Tmax, Phimin, Phimax, fs, duration, NoOfWave, dispout)

Description
-----------

Generate linear waves

Inputs
------

amin
    Min wave amplitude in (m)
amax
    Max wave amplitude in (m)
Tmin
    Min wave mean period in (s)
Tmax
    Max wave mean period in (s)
Phimin=0;
    Min Phase (radian)
Phimax=2*pi;
    Max Phase (radian) 
fs=32;
    Sample generation frequency (Hz), number of data points in one second
duration=10;
    Duration time that data will be generated in (s)
NoOfWave=2;
    Number of waves to be combined with each other
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

Eta
    Water Surface Level Time Series in (m)
t
    Time in (s)
Etaij
    Separated Water Surface Level Time Series in (m)

Examples
--------

.. code:: MATLAB

    [Eta,t,Etaij]=linearwavegenerator(0.2,0.4,1,3,0,2*pi,32,10,2,'yes');

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
        Phimin=0; Phimax=2*pi; fs=32; duration=10; NoOfWave=2; dispout='no';
    case 5
        Phimax=2*pi; fs=32; duration=10; NoOfWave=2; dispout='no';
    case 6
        fs=32; duration=10; NoOfWave=2; dispout='no';
    case 7
        duration=10; NoOfWave=2; dispout='no';
    case 8
        NoOfWave=2; dispout='no';
    case 9
        dispout='no';
end

%--------------------------------------------------------------------------

sample=fs*duration; %number of sample in input file
dt=1/fs;
t(:,1)=linspace(dt,duration,sample); %time
a(:,1)=linspace(amin,amax,NoOfWave);
T(:,1)=linspace(Tmin,Tmax,NoOfWave);
Phi(:,1)=linspace(Phimin,Phimax,NoOfWave);

w=2.*pi./T; %Wave angular frequency

Etaij=zeros(length(t(:,1)),NoOfWave); %Pre-assigning array to make program run faster
for i=1:NoOfWave
    Etaij(:,i)=a(i,1).*cos(-w(i,1).*t+Phi(i,1));
end

Eta=sum(Etaij,2);

%for i=1:length(t(:,1))
%    Eta(i,1)=sum(Etaij(i,1:end));
%end

% -------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1
    
    subplot(2,1,1)
    for i=1:NoOfWave
        plot(t,Etaij(:,i))
        hold on
    end
    title('Generated Waves')
    xlabel('Time (s)')
    ylabel('Water Level (m)')
    xlim([t(1,1) t(end,1)])
    
    subplot(2,1,2)
    plot(t,Eta)
    hold on
    title('Generated Combined Wave')
    xlabel('Time (s)')
    ylabel('Water Level (m)')
    % legend('Water Level (m)')
    xlim([t(1,1) t(end,1)])
    
end

% -------------------------------------------------------------------------
