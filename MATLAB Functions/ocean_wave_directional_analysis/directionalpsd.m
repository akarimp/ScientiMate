function [Syy2d, f2d, theta] = directionalpsd(Syy, f, Wavedir, calcMethod, Windvel, dtheta, dispout)
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

directionalpsd
==============

.. code:: MATLAB

    [Syy2d, f2d, theta] = directionalpsd(Syy, f, Wavedir, calcMethod, Windvel, dtheta, dispout)

Description
-----------

| Calculate wave directional spectrum using parametric directional spreading function
| such as Mitsuyasu et al. (1975), Hasselmann et al. (1980) and Donelan et al. (1985)

Inputs
------

Syy
    One dimensional power spectral density (m^2/Hz)
f
    Frequency (Hz)
Wavedir=0;
    Mean wave direction between 0 and 360 (Degree)
CalcMethod='mitsuyasu';
    | Directional wave spectrum calculation method 
    | 'pierson': Pierson et al. (1955), 'cos2': D=1/pi*(cos((theta-theta_mean)/2))^2 
    | 'mitsuyasu': Mitsuyasu (1975), 'hasselmann': Hasselmann (1980), 'donelan': Donelan et al. (1985) 
    | 'flat': uniform distribution in all directions
Windvel=0;
    | Wind velocity at 10 meter above surface level in (m/s)
    | Wind velocity is required for Mitsuyasu (1975) and Hasselmann (1980) methods
    | For other methods use Windvel=0
dtheta=15;
    Direction interval at which directional spectrum calculated between 0 and 360 (Degree)
dispout='no';
    | Define to display outputs or not
    | '2d': 2 dimensional plot, 'surface': Surface plot, 'polar': Polar plot, 'no': not display 

Outputs
-------

Syy2d
    Directional wave power spectral density (m^2/Hz/Degree)
f2d
    Directional frequency (Hz)
theta
    Direction (Degree)

Examples
--------

.. code:: MATLAB

    N=2^11; %Total number of points
    fs=8; %Sampling frequency
    df=fs/N; %Frequency difference 
    f(:,1)=[0:df:fs/2]; %Frequency vector 
    f(1,1)=f(2,1)/2; %Assign temporarily non-zero value to fisrt element of f to prevent division by zero
    Syy=0.016.*9.81.^2./((2.*pi).^4.*(f.^5)).*exp(-1.25.*(0.33./f).^4); %Calculating Spectrum 
    f(1,1)=0;
    Syy(1,1)=0;
    [Syy2d,f2d,theta]=directionalpsd(Syy,f,45,'mitsuyasu',10,15,'polar');

References
----------

Banner, M. L. (1990). 
Equilibrium spectra of wind waves. 
Journal of Physical Oceanography, 20(7), 966-984.

Donelan, M. A., Hamilton, J., & Hui, W. (1985). 
Directional spectra of wind-generated waves. 
Philosophical Transactions of the Royal Society of London A: Mathematical, Physical and Engineering Sciences, 315(1534), 509-562.

Ewans, K. C. (1998). 
Observations of the directional spectrum of fetch-limited waves. 
Journal of Physical Oceanography, 28(3), 495-512.

Goda, Y. (1999). 
A comparative review on the functional forms of directional wave spectrum. 
Coastal Engineering Journal, 41(01), 1-20.

Hasselmann, D. E., Dunckel, M., & Ewing, J. A. (1980). 
Directional wave spectra observed during JONSWAP 1973. 
Journal of physical oceanography, 10(8), 1264-1280.

Hwang, P. A., & Wang, D. W. (2001). 
Directional distributions and mean square slopes in the equilibrium and saturation ranges of the wave spectrum. 
Journal of physical oceanography, 31(5), 1346-1360.

Mitsuyasu, H., Tasai, F., Suhara, T., Mizuno, S., Ohkusu, M., Honda, T., & Rikiishi, K. (1975). 
Observations of the directional spectrum of ocean WavesUsing a cloverleaf buoy. 
Journal of Physical Oceanography, 5(4), 750-760.

Pierson, W. J., Neumann, G., & James, R. W. (1955). 
Practical methods for observing and forecasting ocean waves by means of wave spectra and statistics. 
Publication 603, U.S. Navy Hydrographic Office, 284 pp. 

Sorensen, R. M. (2005). 
Basic coastal engineering (Vol. 10). 
Springer Science & Business Media.

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
        Wavedir=0; calcMethod='mitsuyasu'; Windvel=0; dtheta=15; dispout='no';
    case 3
        calcMethod='mitsuyasu'; Windvel=0; dtheta=15; dispout='no';
    case 4
        Windvel=0; dtheta=15; dispout='no';
    case 5
        dtheta=15; dispout='no';
    case 6
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(Syy)==1
    Syy=Syy';
end

if isrow(f)==1
    f=f';
end

%--------------------------------------------------------------------------

df=f(2,1)-f(1,1); %Frequency interval
Wavedirrad=deg2rad(Wavedir); %Converting to radian

%--------------------------------------------------------------------------
%Calculating wave properties

%Calculating spectral moments
m0=sum(Syy.*f.^0*df);
m1=sum(Syy.*f.^1*df);
m2=sum(Syy.*f.^2*df);

%Calculating wave properties
Hm0=4*sqrt(m0); %Zero-Moment wave height
Tm01=m0/m1; %Mean period
Tm02=(m0/m2)^0.5; %Zero crossing period

%Calculation peak period
[Syymax SyymaxIndx]=max(Syy(:,1));
Tp=1/f(SyymaxIndx,1); %Peak period
fp=1/Tp; %Peak frequency

%Calculating peak frequency from weighted integral (Young, 1995)
% fp=(sum(Syy.^5.*f.^1*df))./(sum(Syy.^5.*f.^0*df)); %Peak frequency

%--------------------------------------------------------------------------
%Generating direction vector

%thetadeg(:,1)=linspace(0,360,(360+dtheta)/dtheta);
thetadeg(:,1)=[0:dtheta:360]; %Direction vector from 0 to 360 Degree, equally spaced at dtheta
thetarad=deg2rad(thetadeg); %Converting to radian

%--------------------------------------------------------------------------
%Directioal Spectrum Calculation (Ewans, 1998; Goda, 1999; Hwang and Wang, 2001; Sorensen, 2005)

%Wave phase velocity (wave celerity) in a deep water at peak wave frequency
Cp=9.81/(2*pi*fp);

%Calculating directional spreading function
if strcmp(calcMethod,'pierson')==1 %Pierson et al. (1955)
    
    %Directional spreading function: D(f,theta)
    %St. Denis and Pierson (1953): D(f,theta)=2/pi*(cos(theta-theta_mean))^2

    Gs=2/pi;

    %Calculating directional spreading function
    f2d=zeros(length(f(:,1)),length(thetadeg(:,1))); %Pre-assigning vector
    D=zeros(length(f(:,1)),length(thetadeg(:,1))); %Pre-assigning vector
    for j=1:length(thetadeg(:,1))
        f2d(:,j)=f; %Two dimensional frequency vector
        D(:,j)=(Gs.*abs(cos(thetarad(j,1)-Wavedirrad)).^(2));
        if ((Wavedirrad>0) & (Wavedirrad<=pi/2)) 
            if (abs((thetarad(j,1)-1*pi)-(Wavedirrad))<pi/2) 
                D(:,j)=0; %D=0 for abs(theta)>pi/2
            end
        elseif ((Wavedirrad>pi/2) & (Wavedirrad<=3*pi/2)) 
            if (abs(thetarad(j,1)-Wavedirrad)>pi/2) 
                D(:,j)=0; %D=0 for abs(theta)>pi/2
            end
        elseif ((Wavedirrad>3*pi/2) & (Wavedirrad<=2*pi)) 
            if (abs((thetarad(j,1)+1*pi)-(Wavedirrad))<pi/2) 
                D(:,j)=0; %D=0 for abs(theta)>pi/2
            end
        end
    end

elseif strcmp(calcMethod,'cos2')==1 %Mitsuyasu et al. (1975) for s=1
    
    %Directional spreading function: D(f,theta)
    %Mitsuyasu et al. (1975): D(f,theta)=G(s)*(cos((theta-theta_mean)/2))^2s
    %If s=1, then D(f,theta)=G(s)*(cos((theta-theta_mean)/2))^2
    %or, D(f,theta)=1/pi*(cos((theta-theta_mean)/2))^2

    %Calculating spreading parameter, s
    s=1;

    Gs=2.^(2*s-1)./(pi).*(gamma(s+1)).^2./(gamma(2*s+1));

    %Calculating directional spreading function
    f2d=zeros(length(f(:,1)),length(thetadeg(:,1))); %Pre-assigning vector
    D=zeros(length(f(:,1)),length(thetadeg(:,1))); %Pre-assigning vector
    for j=1:length(thetadeg(:,1))
        f2d(:,j)=f; %Two dimensional frequency vector
        D(:,j)=(Gs.*abs(cos((thetarad(j,1)-Wavedirrad)/2)).^(2*s));
    end

elseif strcmp(calcMethod,'mitsuyasu')==1 %Mitsuyasu et al. (1975)

    Sp=11.5*(Windvel/Cp)^(-2.5);

    %Calculating spreading parameter, s
    s=zeros(length(f(:,1)),1); %Pre-assigning s vector
    s(f<fp)=Sp.*(f(f<fp)./fp).^5;
    s(f>=fp)=Sp.*(f(f>=fp)./fp).^(-2.5);

    Gs=2.^(2*s-1)./(pi).*(gamma(s+1)).^2./(gamma(2*s+1));

    %Calculating directional spreading function
    %f2d=zeros(length(f(:,1)),length(thetadeg(:,1))); %Pre-assigning vector
    %D=zeros(length(f(:,1)),length(thetadeg(:,1))); %Pre-assigning vector
    %for i=1:length(f(:,1))
    %    for j=1:length(thetadeg(:,1))
    %        f2d(i,j)=f(i,1); %Two dimensional frequency vector
    %        D(i,j)=real(Gs(i,1)*abs(cos((thetarad(j,1)-Wavedirrad)/2)).^(2*s(i,1)));
    %    end
    %end

    %Calculating directional spreading function
    f2d=zeros(length(f(:,1)),length(thetadeg(:,1))); %Pre-assigning vector
    D=zeros(length(f(:,1)),length(thetadeg(:,1))); %Pre-assigning vector
    for j=1:length(thetadeg(:,1))
        f2d(:,j)=f; %Two dimensional frequency vector
        D(:,j)=(Gs.*abs(cos((thetarad(j,1)-Wavedirrad)/2)).^(2*s));
    end

elseif strcmp(calcMethod,'hasselmann')==1 %Hasselmann et al. (1980)

    mu=-2.33-1.45*(Windvel/Cp-1.17);

    %Calculating spreading parameter, s
    s=zeros(length(f(:,1)),1); %Pre-assigning s vector
    s(f<1.05*fp)=6.97.*(f(f<1.05*fp)./fp).^4.06;
    s(f>=1.05*fp)=9.77.*(f(f>=1.05*fp)./fp).^mu;

    Gs=2.^(2*s-1)./(pi).*(gamma(s+1)).^2./(gamma(2*s+1));

    %Calculating directional spreading function
    f2d=zeros(length(f(:,1)),length(thetadeg(:,1))); %Pre-assigning vector
    D=zeros(length(f(:,1)),length(thetadeg(:,1))); %Pre-assigning vector
    for j=1:length(thetadeg(:,1))
        f2d(:,j)=f; %Two dimensional frequency vector
        D(:,j)=(Gs.*abs(cos((thetarad(j,1)-Wavedirrad)/2)).^(2*s));
    end

elseif strcmp(calcMethod,'donelan')==1 %Donelan et al. (1985)

    %Calculating spreading parameter, Beta,  Donelan et al. (1985) and Banner (1990)
    Beta=zeros(length(f(:,1)),1); %Pre-assigning vector
    Beta(f<=0.56*fp)=1.24;
    Beta(f>0.56*fp & f<=0.95*fp)=2.61.*(f(f>0.56*fp & f<=0.95*fp)./fp).^1.3;
    Beta(f>0.95*fp & f<=1.6*fp)=2.28.*(f(f>0.95*fp & f<=1.6*fp)./fp).^-1.3;
    Beta(f>1.6*fp)=10.^(-0.4+0.8393.*exp(-0.567.*log((f(f>1.6*fp)./fp).^2)));

    %Calculating directional spreading function
    f2d=zeros(length(f(:,1)),length(thetadeg(:,1))); %Pre-assigning vector
    D=zeros(length(f(:,1)),length(thetadeg(:,1))); %Pre-assigning vector
    for j=1:length(thetadeg(:,1))
        f2d(:,j)=f; %Two dimensional frequency vector
        D(:,j)=(0.5.*Beta.*(sech(Beta.*(thetarad(j,1)-Wavedirrad))).^(2));
        if ((thetarad(j,1)-Wavedirrad)>pi & (thetarad(j,1)-Wavedirrad)<2*pi) 
            D(:,j)=(0.5.*Beta.*(sech(Beta.*((thetarad(j,1)-2*pi)-Wavedirrad))).^(2));
        end
    end

elseif strcmp(calcMethod,'flat')==1 %Uniform distribution in all directions 
    
    Gs=1/(2*pi);

    %Calculating directional spreading function
    f2d=zeros(length(f(:,1)),length(thetadeg(:,1))); %Pre-assigning vector
    D=zeros(length(f(:,1)),length(thetadeg(:,1))); %Pre-assigning vector
    for j=1:length(thetadeg(:,1))
        f2d(:,j)=f; %Two dimensional frequency vector
        D(:,j)=Gs;
    end

end

%Normalizing directional spreading function
for i=1:length(f(:,1))
    D(i,:)=D(i,:)./(sum(D(i,:)).*dtheta);
end

%Calculating directional power specral density
for i=1:length(f(:,1))
    for j=1:length(thetadeg(:,1))
        Syy2d(i,j)=Syy(i,1)*D(i,j); %Directional spectrum
    end
end

%--------------------------------------------------------------------------
%Calculating 1-D spectrum by integrating from directional spectrum

for i=1:length(f(:,1))
    Syy1d(i,1)=sum(real(Syy2d(i,:))*dtheta);
    D1d(i,1)=sum(real(D(i,:))*dtheta);
end

%--------------------------------------------------------------------------
%Assigning output parameter

theta=thetadeg; %Direction vector

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'2d')==1 %2 dimensional plot

    val=[Hm0 fp Tp Tm01 Tm02 m0 m1 m2];
    name={'Hm0','fp','Tp','Tm01','Tm02','m0','m1','m2'};
    for i=1:length(val)
        fprintf('%14s   %g\n',name{i},val(i));
    end

    %Plotting
    plot(f2d,Syy2d)

    title('Power Spectral Density')
    xlabel('Frequency(Hz)')
    ylabel('Spectral Density(m^2/Hz)')

elseif strcmp(dispout,'surface')==1 %Surface plot

    val=[Hm0 fp Tp Tm01 Tm02 m0 m1 m2];
    name={'Hm0','fp','Tp','Tm01','Tm02','m0','m1','m2'};
    for i=1:length(val)
        fprintf('%14s   %g\n',name{i},val(i));
    end

    [thetaGrid,fGrid]=meshgrid(thetadeg(:,1),f(:,1));
    surf1=surf(fGrid,thetaGrid,Syy2d);
    %set(surf1,'xscale','log','zscale','log')
    %set(surf1,'FaceColor','interp','EdgeColor','none')
    set(surf1,'EdgeColor','none')
    colorbar
    title('Directional Power Spectral Density')
    xlabel('Frequency (Hz)')
    ylabel('Direction (Degree)')
    zlabel('Power Spectral Density (m^2/Hz/degree)')

elseif strcmp(dispout,'polar')==1 %Polar plot

    val=[Hm0 fp Tp Tm01 Tm02 m0 m1 m2];
    name={'Hm0','fp','Tp','Tm01','Tm02','m0','m1','m2'};
    for i=1:length(val)
        fprintf('%14s   %g\n',name{i},val(i));
    end

    [thetaGrid,fGrid]=meshgrid(thetadeg(:,1),f(:,1));
    [X,Y]=pol2cart(deg2rad(thetaGrid),fGrid);
    %h=polar([0 2*pi], [0 fs],':k');
    pol1=polar(X,Y);
    delete(pol1);
    hold on;
    contour(X,Y,Syy2d);
    colorbar
    title('Directional Power Spectral Density')
    %xlabel('Direction (Degree)')
    %ylabel('Frequency (Hz)')
    %zlabel('Spectral Density (m^2/Hz/degree)')

end

%--------------------------------------------------------------------------
