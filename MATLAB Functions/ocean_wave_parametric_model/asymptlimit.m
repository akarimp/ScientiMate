function [EhatAsympt, fphatAsympt, hhat] = asymptlimit(windvel, hmean, CalcMethod)
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

asymptlimit
===========

.. code:: MATLAB

    [EhatAsympt, fphatAsympt, hhat] = asymptlimit(windvel, hmean, CalcMethod)

Description
-----------

Calculate dimensionless asymptotic limits of wind wave growth in shallow and intermediate water

Inputs
------

windvel
    | Wind velocity in (m/s)
    | Wind velocity should be measured (or represents velocity) at 10 m above surface
    | Wind velocity should be a 10 minutes averaged wind for all methods, except for for 'spmshallow' methods
    | For 'spmshallow' methods, wind velocity should be converted to duration of sustained wind by using gust factor
Fetch
    Wind fetch in (m)
hmean
    Mean water depth along a wind fetch in (m)
CalcMethod='karimpour';
    | Parametric wave model to be used for calculation 
    | 'smb': Use SMB method by Munk (1947), and Bretschneider (1952, 1958)
    | 'young': Use method by Young and Verhagen (1996) and Young and Babanin (2006)
    | 'karimpour': Use method by Karimpour et al. (2017)

Outputs
-------

EhatAsympt
    Dimensionless asymptotic limit values for wave energy, Ehat=g^2*m0/U10^4
fphatAsympt
    Dimensionless asymptotic limit values for peak wave frequency, fphat=fp*U10/9.81 
hhat
    | Dimensionless water depth along a wind fetch: hhat=g*hmean/U10^2
    | Note, g=9.81: gravitational acceleration
    |     U10: wind velocity
    |     fp: peak wave frequency
    |     hmean: mean depth along a wind fetch

Examples
--------

.. code:: MATLAB

    windvel(:,1)=10.*rand(100,1);
    hmean(:,1)=3.*rand(100,1);
    [EhatAsympt,fphatAsympt,hhat]=asymptlimit(windvel,hmean,'karimpour');

References
----------

Bretschneider, C. L. (1952). 
Revised wave forecasting relationships. 
Coastal Engineering Proceedings, 1(2), 1.

Bretschneider, C. L. (1958). 
Revisions in wave forecasting: deep and shallow water. 
Coastal Engineering Proceedings, 1(6), 3.

Karimpour, A., Chen, Q., & Twilley, R. R. (2017). 
Wind Wave Behavior in Fetch and Depth Limited Estuaries. 
Scientific reports, 7, 40654.

Sverdrup, H. U., & Munk, W. H. (1947). 
Wind, sea, and swell: theory of relations for forecasting. 
U.S. Navy Department, Hydrographic Office, Publication No. 601, 44 pp. 

Young, I. R., & Verhagen, L. A. (1996). 
The growth of fetch limited waves in water of finite depth. Part 1. Total energy and peak frequency. 
Coastal Engineering, 29(1-2), 47-78.

Young, I. R., & Babanin, A. V. (2006). 
The form of the asymptotic depth‐limited wind wave frequency spectrum. 
Journal of Geophysical Research: Oceans, 111(C6).

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

%--------------------------------------------------------------------------
%Calculating dimensionless asymptotic limits

%Calculating dimensionless water depth: hhat=g*hmean/U10^2
hhat=9.81.*hmean./(windvel.^2);  %Dimensionless depth: hhat=g*hmean/U10^2

%Calculating Dimensionless Asymptotic Limits using SMB method (Munk, 1947; and Bretschneider, 1952, 1958)
if strcmp(CalcMethod,'smb')==1

    EhatAsympt=1.4e-3.*hhat.^1.5; %Asymptotic Limit for Ehat from SMB
    fphatAsympt=0.16.*hhat.^-0.375; %Asymptotic Limit for fphat from SMB

%Calculating wave properties using Young and Verhagen (1996) and Young and Babanin (2006)
elseif strcmp(CalcMethod,'young')==1

    EhatAsympt=1e-3.*hhat.^1.2; %Asymptotic Limit from Young and Babanin (2006)
    fphatAsympt=0.2.*hhat.^-0.375; %Asymptotic Limit from Young and Verhagen (1996)

%Calculating wave properties using Karimpour et al. (2017)
elseif strcmp(CalcMethod,'karimpour')==1

    EhatAsympt=3e-5.*tan(1.56255.*(tanh(3.356.*hhat)).^0.315); %Asymptotic Limit for Ehat
    fphatAsympt=0.133.*(tanh(0.832.*hhat)).^-0.716.*(tanh(2.623.*hhat)).^0.461; %Asymptotic Limit for fphat

end

%--------------------------------------------------------------------------
