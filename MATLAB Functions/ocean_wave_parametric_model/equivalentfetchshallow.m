function [FetchEq, FetchhatEq] = equivalentfetchshallow(windvel, hmean, SustWindDur, CalcMethod)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-09-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

equivalentfetchshallow
======================

.. code:: MATLAB

    [FetchEq, FetchhatEq] = equivalentfetchshallow(windvel, hmean, SustWindDur, CalcMethod)

Description
-----------

Calculate an equivalent wind fetch for duration limited wave growth in shallow and intermediate water

Inputs
------

windvel
    | Wind velocity in (m/s)
    | Wind velocity should be measured (or represents velocity) at 10 m above surface
    | Wind velocity should be a 10 minutes averaged wind for all methods, except for for 'cem' and 'spmdeep' methods
    | For 'cem' and 'spmdeep' methods, wind velocity should be converted to duration of sustained wind by using gust factor
hmean
    Mean water depth along a wind fetch in (m)
SustWindDur
    Sustained wind duration in (second)
CalcMethod='karimpour';
    | Parametric wave model to be used for calculation 
    | 'spmshallow': Use method by Shore Protection Manual (SPM),
    |     U.S. Army Corps of Engineers (1984) in shallow and intermediate water water
    | 'karimpour': Use method by Karimpour et al. (2017)

Outputs
-------

FetchEq
    Equivalent wind fetch values for a duration limited wave growth in (m)
FetchhatEq
    | Dimensionless equivalent wind fetch values for a duration limited wave growth: Fetchhat=g*Fetch/U10^2
    | Note, g=9.81: gravitational acceleration
    |     U10: wind velocity
    |     Fetch: wind fetch

Examples
--------

.. code:: MATLAB

    windvel(:,1)=10.*rand(100,1);
    SustWindDur(:,1)=10000.*rand(100,1);
    hmean(:,1)=3.*rand(100,1);
    [FetchEq,FetchhatEq]=equivalentfetchshallow(windvel,hmean,SustWindDur,'karimpour');

References
----------

U.S. Army Corps of Engineers (2015). 
Coastal Engineering Manual. 
Engineer Manual 1110-2-1100, Washington, D.C.: U.S. Army Corps of Engineers.

Department of the Army, Waterways Experiment Station, Corps of Engineers, 
and Coastal Engineering Research Center (1984), 
Shore Protection Manual, Washington, 
D.C., vol. 1, 4th ed., 532 pp.

Karimpour, A., Chen, Q., & Twilley, R. R. (2017). 
Wind Wave Behavior in Fetch and Depth Limited Estuaries. 
Scientific reports, 7, 40654.

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
        CalcMethod='karimpour';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(windvel)==1
    windvel=windvel';
end

if isrow(hmean)==1
    hmean=hmean';
end

if isrow(SustWindDur)==1
    SustWindDur=SustWindDur';
end

%--------------------------------------------------------------------------
%Calculating an equivalent wind fetch for duration limited wave growth

%Calculating an equivalent wind fetch
%using Shore Protection Manual (SPM), U.S. Army Corps of Engineers (1984) in shallow and intermediate water
if strcmp(CalcMethod,'spmshallow')==1

    %Calculating adjusted wind velcity
    UA=0.71.*windvel.^1.23; %Wind stress factor or adjusted wind velcity
    
    %Calculating dimensionless numbers from dimensional number
    hhat=9.81.*hmean./(windvel.^2); %Dimensionless depth: hhat=g*hmean/U10^2
    hhatUA=9.81.*hmean./(UA.^2); %Dimensionless Depth, hhatUA=g.hmean/UA^2

    Tshat=(9.81.*SustWindDur./UA./5.37e2).^(3/7); %Dimensionless significant wave period, Tshat=g*Ts/UA

    FetchhatUA=(abs(atanh(Tshat./(7.54.*tanh(0.833.*(hhatUA).^(3/8))))).*(tanh(0.833.*(hhatUA).^(3/8)))./0.0379).^(3); %Dimensionless Fetch, FetchhatUA=g.Fetch/UA^2
    FetchEq=FetchhatUA.*(UA.^2)./9.81; %Equivalent Fetch
    FetchhatEq=9.81.*FetchEq./(windvel.^2); %Dimensionless Equivalent Fetch


%Calculating an equivalent wind fetch using Karimpour et al. (2017)
elseif strcmp(CalcMethod,'karimpour')==1

    hhat=9.81.*hmean./(windvel.^2); %Dimensionless depth: hhat=g*hmean/U10^2

    SustWindDurhat=9.81.*SustWindDur./windvel; %Dimensionless Sustained Wind Duration
    FetchhatEq=SustWindDurhat.*hhat./(2.59.*hhat.^(2/3)); %Dimensionless Equivalent Fetch: Fetchhat=g*Fetch/U10^2
    FetchEq=FetchhatEq.*(windvel.^2)./9.81; %Equivalent Fetch


end

%--------------------------------------------------------------------------