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

hurricanewindvelmaxh08
======================

.. code:: MATLAB

    [Vmax, Vgmax] = hurricanewindvelmaxh08(y, Pc, dPcdt, Vt, Pn, SST, VgToVCoeff, dispout)

Description
-----------

Calculate hurricane maximum wind velocity at the gradient level using Holland (2008) method

Inputs
------

y
    y (latitude) of points which outputs are calculated at (Degree)
Pc
    Hurricane central surface pressure in (Pa)
dPcdt
    | Hurricane central pressure (Pc) intensity change (dPc/dt) in (hPa/hr)
    | Convert from (Pa) to (hPa): hPa=0.01*Pa;
Vt
    Hurricane translation velocity in (m/s)
Pn=101325;
    | Ambient surface pressure (external pressure) in (Pa)
    | Standard atmosphere pressure is 101325 (Pa) 
    | Typical value: Pn=101500 (Pa) for the western North Pacific, Pn= 101000 (Pa) for the North Atlantic
    | (Batke et al., 2014)
SST=28;
    Sea surface temperature in (C)
VgToVCoeff=0.8;
    | Coefficient to convert gradient wind velocity to wind velocity at 10 m above surface as: 
    | V=VgToVCoeff*Vg, if VgToVCoeff=1, then V=Vg
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

Vmax
    Maximum hurricane 1-min wind velocity at the surface level (at 10-m height) in (m/s)
Vgmax
    | Maximum hurricane 1-min wind velocity at the gradient level in (m/s)
    | Vgmax is calculated from Vmax as Vgmax=Vmax/VgToVCoeff
    | Gradient wind velocity can be converted to wind velocity at 10 m above surface by V=VgToVCoeff*Vg
    | VgToVCoeff can be approximated as VgToVCoeff=0.8
    | For detail on converting gradient (mean boundary-layer) wind velocity to velocity 10 m above surface see
    | e.g. Graham & Numm (1959), Young & Vinoth (2013), Phadke et al. (2003), Powell et al. (2003), Valamanesh et al. (2016), Wei et al. (2017)
    | Shore Protection Manual (SPM), U.S. Army Corps of Engineers (1984)
    | Coastal Engineering Manual (CEM), U.S. Army Corps of Engineers (2015)

Examples
--------

.. code:: MATLAB

    %EXAMPLE 1

    y=26.3; %(Degree)
    Pc=90200; %(Pa)
    Vt=5; %Hurricane translation velocity in (m/s)
    dPcdt=3; %Hurricane central pressure (Pc) intensity change in (hPa/hr)
    Pn=101325; %Ambient surface pressure (external pressure) in (Pa)
    SST=27.77; %Sea surface temperature in (C)
    VgToVCoeff=0.8;
    [Vmax,Vgmax]=hurricanewindvelmaxh08(y,Pc,dPcdt,Vt,Pn,SST,VgToVCoeff,'no');


    %EXAMPLE 2

    %Latitude of Hurricane Katrine best track
    lattrack=[23.1;23.4;23.8;24.5;25.4;26.0;26.1;26.2;26.2;26.0;25.9;25.4;...
        25.1;24.9;24.6;24.4;24.4;24.5;24.8;25.2;25.7;26.3;27.2;28.2;...
        29.3;29.5;30.2;31.1;32.6;34.1;35.6;37.0;38.6;40.1];

    %Hurricane Katrina centeral pressure (Pa)
    Pc=[100800;100700;100700;100600;100300;100000;99700;99400;98800;98400;98300;98700;...
        97900;96800;95900;95000;94200;94800;94100;93000;90900;90200;90500;91300;...
        92000;92300;92800;94800;96100;97800;98500;99000;99400;99600];

    %Hurricane Katrina translational velocity (m/s)
    Vt=[0.00000;3.23091;3.13105;3.86928;4.99513;4.82816;3.27813;2.81998;2.77140;2.53041;...
        1.05928;5.30662;3.60661;2.98269;3.61863;3.43691;3.28168;2.85849;3.20404;4.26279;...
        5.31340;5.18467;5.39195;5.46121;5.66270;1.02958;3.60354;4.63312;8.02540;8.01558;...
        8.12721;8.31580;10.75406;12.28350];
        
    %Hurricane central pressure (Pc) intensity change in (hPa/hr)
    dPcdt=[0.00000;-0.16667;0.00000;-0.16667;-0.50000;-0.50000;-0.50000;-0.50000;-1.00000;-0.66667;-0.16667;...
        0.66667;-1.33333;-1.83333;-1.50000;-1.50000;-1.33333;1.00000;-1.16667;-1.83333;-3.50000;-1.16667;...
        0.50000;1.33333;1.16667;0.50000;0.83333;3.33333;2.16667;2.83333;1.16667;0.83333;0.66667;...
        0.33333];

    Pn=101325; %Ambient surface pressure (external pressure) in (Pa)
    SST=27.77; %Sea surface temperature in (C)
    VgToVCoeff=0.8;

    [Vmax,Vgmax]=hurricanewindvelmaxh08(lattrack,Pc,dPcdt,Vt,Pn,SST,VgToVCoeff,'yes');

References
----------

Data

* www.nhc.noaa.gov/data/
* www.nhc.noaa.gov/data/hurdat/hurdat2-format-nencpac.pdf
* coast.noaa.gov/hurricanes
* www.aoml.noaa.gov/hrd/data_sub/re_anal.html

Batke, S. P., Jocque, M., & Kelly, D. L. (2014). 
Modelling hurricane exposure and wind speed on a mesoclimate scale: a case study from Cusuco NP, Honduras. 
PloS one, 9(3), e91306.

Department of the Army, Waterways Experiment Station, Corps of Engineers, 
and Coastal Engineering Research Center (1984), 
Shore Protection Manual, Washington, 
D.C., vol. 1, 4th ed., 532 pp.

Graham and Numm (1959) 
Meteorological Conditions Pertinent to Standard Project Hurricane, Atlantic and Gulf Coasts of United States.
National Hurricane Research Project. U.S. Weather Service, Report no. 33.

Harper, B. A., & Holland, G. J. (1999, January). 
An updated parametric model of the tropical cyclone. 
In Proc. 23rd Conf. Hurricanes and Tropical Meteorology.

Holland, G. (2008). 
A revised hurricane pressure–wind model. 
Monthly Weather Review, 136(9), 3432-3445.

Phadke, A. C., Martino, C. D., Cheung, K. F., & Houston, S. H. (2003). 
Modeling of tropical cyclone winds and waves for emergency management. 
Ocean Engineering, 30(4), 553-578.

Powell, M. D., Vickery, P. J., & Reinhold, T. A. (2003). 
Reduced drag coefficient for high wind speeds in tropical cyclones. 
Nature, 422(6929), 279.

U.S. Army Corps of Engineers (2015). 
Coastal Engineering Manual. 
Engineer Manual 1110-2-1100, Washington, D.C.: U.S. Army Corps of Engineers.

Valamanesh, V., Myers, A. T., Arwade, S. R., Hajjar, J. F., Hines, E., & Pang, W. (2016). 
Wind-wave prediction equations for probabilistic offshore hurricane hazard analysis. 
Natural Hazards, 83(1), 541-562.

Wei, K., Arwade, S. R., Myers, A. T., Valamanesh, V., & Pang, W. (2017). 
Effect of wind and wave directionality on the structural performance of non‐operational offshore wind turbines supported by jackets during hurricanes. 
Wind Energy, 20(2), 289-303.

Young, I. R., & Vinoth, J. (2013). 
An 'extended fetch' model for the spatial distribution of tropical cyclone wind–waves as observed by altimeter. 
Ocean Engineering, 70, 14-24.

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
