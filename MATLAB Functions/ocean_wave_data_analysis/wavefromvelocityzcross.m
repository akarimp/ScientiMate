function [Hs, Ts, Hz, Tz, Hrms, H, T, Eta, t] = wavefromvelocityzcross(Ux, Uy, fs, h, heightfrombed, Kuvmin, kCalcMethod, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-04-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

wavefromvelocityzcross
======================

.. code:: MATLAB

    [Hs, Ts, Hz, Tz, Hrms, H, T, Eta, t] = wavefromvelocityzcross(Ux, Uy, fs, h, heightfrombed, Kuvmin, kCalcMethod, dispout)

Description
-----------

Calculate wave properties from wave orbital velocity by using an upward zero crossing method

Inputs
------

Ux
    Wave horizontal orbital velocity data in x direction in (m/s)
Uy
    Wave horizontal orbital velocity data in y direction in (m/s)
fs
    Sampling frequency that data collected at in (Hz)
h
    Water depth in (m)
heightfrombed=0;
    Height from bed that data collected at in (m)
Kuvmin=0.15;
    | Minimum acceptable value for an orbital velocity converstion factor
    | If Kuvmin=0.15, it avoid wave amplification larger than 6 times (1/0.15)
kCalcMethod='beji';
    | Wave number calculation method 
    | 'hunt': Hunt (1979), 'beji': Beji (2013), 'vatankhah': Vatankhah and Aghashariatmadari (2013) 
    | 'goda': Goda (2010), 'exact': calculate exact value 
Rho=1000;
    Water density (kg/m^3)
dispout='no';
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

Hs
    Significant Wave Height (m)
Ts
    Significant Wave Period (second)
Hz
    Zero Crossing Mean Wave Height (m)
Tz
    Zero Crossing Mean Wave Period (second)
Hrms
    Root Mean Square Wave Height (m)
H
    Wave Height Data Series array (m)
T
    Wave Period Data Series array (second)
Eta
    Water surface elevation time series in (m)
t
    Time (s)

Examples
--------

.. code:: MATLAB

    fs=2; %Sampling frequency
    duration=1024; %Duration of the data
    N=fs*duration; %Total number of points
    df=fs/N; %Frequency difference 
    dt=1/fs; %Time difference, dt=1/fs
    t(:,1)=linspace(0,duration-dt,N); %Time
    Eta(:,1)=detrend(0.5.*cos(2*pi*0.2*t)+(-0.1+(0.1-(-0.1))).*rand(N,1));
    hfrombed=4;
    h=5;
    k=0.2;
    Ux=(pi/5).*(2.*Eta).*(cosh(k*hfrombed)/sinh(k*h)); 
    Uy=0.2.*Ux;
    [Hs,Ts,Hz,Tz,Hrms,H,T,Eta,t]=wavefromvelocityzcross(Ux,Uy,fs,5,4,0.15,'beji','yes');

References
----------

Beji, S. (2013). 
Improved explicit approximation of linear dispersion relationship for gravity waves. 
Coastal Engineering, 73, 11-12.

Goda, Y. (2010). 
Random seas and design of maritime structures. 
World scientific.

Hunt, J. N. (1979). 
Direct solution of wave dispersion equation. 
Journal of the Waterway Port Coastal and Ocean Division, 105(4), 457-459.

Vatankhah, A. R., & Aghashariatmadari, Z. (2013). 
Improved explicit approximation of linear dispersion relationship for gravity waves: A discussion. 
Coastal engineering, 78, 21-22.

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
        heightfrombed=0; Kuvmin=0.15; kCalcMethod='beji'; dispout='no';
    case 5
        Kuvmin=0.15; kCalcMethod='beji'; dispout='no';
    case 6
        kCalcMethod='beji'; dispout='no';
    case 7
        dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

if isrow(Ux)==1
    Ux=Ux';
end

if isrow(Uy)==1
    Uy=Uy';
end

%--------------------------------------------------------------------------

%Deterending input data
UxDetrended=detrend(Ux,'linear');
UyDetrended=detrend(Uy,'linear');
U=sqrt(UxDetrended.^2+UyDetrended.^2);
U=U.*sign(UxDetrended);
UDetrended=detrend(U,'linear');
%h(h<=0)=0.001;

%--------------------------------------------------------------------------

%Generating time vector
N=length(Ux(:,1)); %Length of a time series, total number of points in input file, N=fs*duration
duration=N/fs; %Duration time that data are collected (second), duration=N/fs
dt=1/fs; %Time difference between consecutive samples, dt=1/fs=duration/N
%t(:,1)=linspace(0,duration-dt,N); %Time from 0 to T-dt, equally spaced at dt
t(:,1)=[0:dt:duration-dt]; %Time from 0 to T-dt, equally spaced at dt

%--------------------------------------------------------------------------
%Detecting first and last wave (first and last complete crest-trough cycle)

%Detecting the start point of the first wave (first complete crest-trough cycle)
if UDetrended(1,1)==0 & UDetrended(2,1)>0 
    len3=1;
end

if (UDetrended(1,1)==0 & UDetrended(2,1)<0) | (UDetrended(1,1)~=0)
    for i=1:N-1
        if UDetrended(i,1)<0 & UDetrended(i+1,1)>0 
            len3=i;
            break
        end
    end
end

%Detecting the end point of last wave (last complete crest-trough cycle)
if UDetrended(end,1)==0 & UDetrended(end-1,1)<0 
    len4=N;
end

if (UDetrended(end,1)==0 & UDetrended(end-1,1)>0) | (UDetrended(end,1)~=0)
    for i=N-1:-1:1
        if UDetrended(i,1)<0 & UDetrended(i+1,1)>0 
            len4=i;
            break
        end
    end
end

%--------------------------------------------------------------------------
%Detecting zero crossing points

m=0;
n=0;

for i=len3:len4 
    
    %Detecting up-crossing zero-crossing points
    if i==1
        m=m+1;
        xupcross(m,1)=dt;
        yupcross(m,1)=0;
        xupcrossIndx(m,1)=i;
    end
    
    if i>1 & i<N
        if UDetrended(i,1)<0 & UDetrended(i+1,1)>0
            m=m+1;
            xupcross(m,1)=t(i,1)-(t(i+1,1)-t(i,1))/(UDetrended(i+1,1)-UDetrended(i,1))*UDetrended(i,1);
            yupcross(m,1)=0;
            xupcrossIndx(m,1)=i;
        end
    end
    
    if i==N
        m=m+1;
        xupcross(m,1)=t(end,1);
        yupcross(m,1)=0;
        xupcrossIndx(m,1)=i;
    end
    
    %Detecting down-crossing zero-crossing points
    if i>1 & i<N
        if UDetrended(i,1)>0 & UDetrended(i+1,1)<0
            n=n+1;
            xdowncross(n,1)=t(i,1)-(t(i+1,1)-t(i,1))/(UDetrended(i+1,1)-UDetrended(i,1))*UDetrended(i,1);
            ydowncross(n,1)=0;
            xdowncrossIndx(n,1)=i;
            if xdowncrossIndx(n,1)>=N
                xdowncrossIndx(n,1)=N-1;
            end
        end
    end
    
end

%--------------------------------------------------------------------------
%Detecting crest and trough

m=0;
n=0;

for i=xupcrossIndx(1,1):xupcrossIndx(end,1)
    
    %Detecting crest
    if i>1 & i<N
        
        if UDetrended(i,1)>UDetrended(i-1,1) & UDetrended(i,1)>UDetrended(i+1,1)
            
            %Check if after last crest, a trough is detected or not 
            %m==n means it is detected and the new y(i,1) is the crest of the next wave
            if UDetrended(i,1)>0
                
                if m==n
                    
                    m=m+1;
                    xmax(m,1)=t(i,1);
                    ymax(m,1)=UDetrended(i,1);
                    xmaxIndx(m,1)=i;
                    
                elseif m~=0
                        
                    %Replacingthe old crest location with new one if the new one is larger
                    if UDetrended(i,1)>ymax(m,1) & m~=0
                        
                        xmax(m,1)=t(i,1);
                        ymax(m,1)=UDetrended(i,1);
                        xmaxIndx(m,1)=i;
                            
                    end
                    
                end
                
            end
            
        end
        
    end
    
    %Detecting trough
    if i>1 & i<N
        
        if UDetrended(i,1)<UDetrended(i-1,1) & UDetrended(i,1)<UDetrended(i+1,1)
            
            if UDetrended(i,1)<0
                
                %Check if after last trough, a crest is detected or not 
                %n==m-1 means it is detected and the new y(i,1) is the trough of next wave
                if n==m-1
                    
                    n=n+1;
                    xmin(n,1)=t(i,1);
                    ymin(n,1)=UDetrended(i,1);
                    xminIndx(n,1)=i;
                    
                elseif n~=0
                        
                    %Replacingthe old crest location with new one if the new one is smaller
                    if UDetrended(i,1)<ymin(n,1)
                        
                        xmin(n,1)=t(i,1);
                        ymin(n,1)=UDetrended(i,1);
                        xminIndx(n,1)=i;
                            
                    end
                    
                end
                
            end
            
        end
        
    end
    
    
end

%Eta=UDetrended; %Water surface elevation time series

%--------------------------------------------------------------------------
%Calculating wave height and surface elevation of the wave crest and trough

len1=length(xmax(:,1));
len2=length(xmin(:,1));

if len1~=len2
    error('Numbers of crests and troughs are not the same')
end

for i=1:len1
    
    H(i,1)=ymax(i,1)-ymin(i,1); %Wave height
    xmean(i,1)=(xmax(i,1)+xmin(i,1))/2;
    Etac(i,1)=ymax(i,1); %Water surface elevation of the wave crest
    Etat(i,1)=ymin(i,1); %Water surface elevation of the wave trough
    
end
        
%--------------------------------------------------------------------------
%Calculating wave period 

for i=1:length(H(:,1))
    for j=1:length(xupcross(:,1))-1
        
        if xupcross(j,1)<xmean(i,1) & xupcross(j+1,1)>xmean(i,1)
            T(i,1)=xupcross(j+1,1)-xupcross(j,1);
        end
        
    end
end

%--------------------------------------------------------------------------
%==========================================================================
%--------------------------------------------------------------------------
%Calculating wave number (k)

f=1./T; %Wave frequency
w=2.*pi.*f; %Wave angular frequency, w=2*pi/T=2*pi*f

%Deep water
k0=(w.^2)./9.81; %Deep water wave number
k0h=k0.*h;

%Estimation of wave number (k) from Hunt (1979)
if strcmp(kCalcMethod,'hunt')==1
    % c2gh=(k0h+(1.0+0.6522.*k0h+0.4622.*k0h.^2+0.0864.*k0h.^4+0.0675.*k0h.^5).^(-1)).^(-1); %Calculating wave number from Hunt (1979)
    % k=w./(sqrt(c2gh.*9.81*h));
    d1=0.6666666667; d2=0.3555555556; d3=0.1608465608; d4=0.0632098765; d5=0.0217540484; d6=0.0065407983;
    kh=sqrt(k0h.^2+(k0h)./(1+d1.*k0h+d2.*k0h.^2+d3.*k0h.^3+d4.*k0h.^4+d5.*k0h.^5+d6.*k0h.^6)); %Calculating wave number from Hunt (1979)
    k=kh./h; %Calculating wave number from Hunt (1979)

%Estimation of wave number (k) from Beji (2013)
elseif strcmp(kCalcMethod,'beji')==1
    kh=k0h.*(1+k0h.^1.09.*exp(-(1.55+1.3.*k0h+0.216.*k0h.^2)))./sqrt(tanh(k0h)); %Calculating wave number from Beji (2013)
    k=kh./h; %Calculating wave number from Beji (2013)
    k(w==0)=0;

%Estimation of wave number (k) from Vatankhah and Aghashariatmadari (2013)
elseif strcmp(kCalcMethod,'vatankhah')==1
    kh=(k0h+k0h.^2.*exp(-(3.2+k0h.^1.65)))./sqrt(tanh(k0h))+k0h.*(1-exp(-(k0h.^0.132))).^(5.0532+2.158.*k0h.^1.505); %Calculating wave number from Vatankhah and Aghashariatmadari (2013)
    k=kh./h; %Calculating wave number from Vatankhah and Aghashariatmadari (2013)
    k(w==0)=0;

%Estimation of wave number (k) from Goad (2010)
elseif strcmp(kCalcMethod,'goda')==1
    kh=zeros(length(k0h(:,1)),1);
    kh(k0h>=1)=k0h(k0h>=1);
    kh(k0h<1)=(k0h(k0h<1)).^0.5;
    for i=1:3
        kh=kh-((kh-k0h.*coth(kh))./(1+k0h.*((coth(kh)).^2-1))); %Calculating wave number from Goda (2010)
    end
    k=kh./h; %Calculating wave number from Goda (2010)
    k(w==0)=0;

%Calculating exact wave number (k)
elseif strcmp(kCalcMethod,'exact')==1

    %Estimation of wave number (k) from Beji (2013)
    kh=k0h.*(1+k0h.^1.09.*exp(-(1.55+1.3.*k0h+0.216.*k0h.^2)))./sqrt(tanh(k0h)); %Calculating wave number from Beji (2013)
    kini=kh./h; %Initial value for k (Wave number from Beji (2013))
    kini(w==0)=0;

    %Calculating exact wave number (k)
    k=zeros(length(f(:,1)),1); %Pre-assigning array to make program run faster
    for i=1:length(f(:,1))
        if length(h(:,1))==1
            fun = @(x)(w(i,1)^2-(9.81*x*tanh(x*h))); %w^2=g*k*tanh(kh)
        else
            fun = @(x)(w(i,1)^2-(9.81*x*tanh(x*h(i,1)))); %w^2=g*k*tanh(kh)
        end
        k(i,1)=fzero(fun,kini(i,1)); %Wave number
    end

end

%--------------------------------------------------------------------------
%Calculating the velocity response factor, Kuv, (Wiberg, 2008)
Kuv=(2*pi./T).*cosh(k.*heightfrombed)./sinh(k.*h);
Kuv(1,1)=Kuv(2,1);

%Calculating KmaxL based on linear wave theory
kmaxL=pi/(h-heightfrombed); %Wave number associated with fmaxuvcorrL

%Calculating fmaxuvcorrL based on linear wave theory
fmaxuvcorrL=1/(2*pi)*sqrt(9.81*pi/(h-heightfrombed)*tanh(pi*h/(h-heightfrombed))); %fmaxuvcorrL based on linear wave theory

%Comparing Kuv with KuvminL calculated based on linear wave theory
KuvminL=(2*pi*fmaxuvcorrL)*cosh(kmaxL*heightfrombed)/sinh(kmaxL*h); %Minimum Limit for Kuv calculated based on linear wave theory
Kuv(Kuv<KuvminL)=KuvminL; %Check to avoid large amplification, Kuv should be larger than minimum Kuv calculated based on linear wave theory

%Check to avoid wave amplification larger than defined by Kuvmin
Kuv(Kuv<Kuvmin)=Kuvmin; 

%--------------------------------------------------------------------------
%Water surface elevation calculation

for i=1:length(H(:,1))

    for j=xupcrossIndx(i,1):xupcrossIndx(i+1,1)
        EtawithKuv(j,1)=UDetrended(j,1)/Kuv(i,1);
    end

end

%Data before very first wave, first Kuv is used due lack of enough information
for j=1:xupcrossIndx(1,1)
    EtawithKuv(j,1)=UDetrended(j,1)/Kuv(1,1); 
end

%Data after very last wave, last Kuv is used due lack of enough information
for j=xupcrossIndx(end,1):length(UDetrended(:,1))
    EtawithKuv(j,1)=UDetrended(j,1)/Kuv(end,1); 
end

Eta=EtawithKuv; %Water surface elevation time series

%--------------------------------------------------------------------------
%==========================================================================
%--------------------------------------------------------------------------
%Detecting crest and trough

m=0;
n=0;

for i=xupcrossIndx(1,1):xupcrossIndx(end,1)
    
    %Detecting crest
    if i>1 & i<N
        
        if EtawithKuv(i,1)>EtawithKuv(i-1,1) & EtawithKuv(i,1)>EtawithKuv(i+1,1)
            
            %Check if after last crest, a trough is detected or not 
            %m==n means it is detected and the new y(i,1) is the crest of the next wave
            if EtawithKuv(i,1)>0
                
                if m==n
                    
                    m=m+1;
                    xmax(m,1)=t(i,1);
                    ymax(m,1)=EtawithKuv(i,1);
                    xmaxIndx(m,1)=i;
                    
                elseif m~=0
                        
                    %Replacingthe old crest location with new one if the new one is larger
                    if EtawithKuv(i,1)>ymax(m,1) & m~=0
                        
                        xmax(m,1)=t(i,1);
                        ymax(m,1)=EtawithKuv(i,1);
                        xmaxIndx(m,1)=i;
                            
                    end
                    
                end
                
            end
            
        end
        
    end
    
    %Detecting trough
    if i>1 & i<N
        
        if EtawithKuv(i,1)<EtawithKuv(i-1,1) & EtawithKuv(i,1)<EtawithKuv(i+1,1)
            
            if EtawithKuv(i,1)<0
                
                %Check if after last trough, a crest is detected or not 
                %n==m-1 means it is detected and the new y(i,1) is the trough of next wave
                if n==m-1
                    
                    n=n+1;
                    xmin(n,1)=t(i,1);
                    ymin(n,1)=EtawithKuv(i,1);
                    xminIndx(n,1)=i;
                    
                elseif n~=0
                        
                    %Replacingthe old crest location with new one if the new one is smaller
                    if EtawithKuv(i,1)<ymin(n,1)
                        
                        xmin(n,1)=t(i,1);
                        ymin(n,1)=EtawithKuv(i,1);
                        xminIndx(n,1)=i;
                            
                    end
                    
                end
                
            end
            
        end
        
    end
    
    
end

%Eta=EtawithKuv; %Water surface elevation time series

%--------------------------------------------------------------------------
%Calculating wave height and surface elevation of the wave crest and trough

len1=length(xmax(:,1));
len2=length(xmin(:,1));

if len1~=len2
    error('Numbers of crests and troughs are not the same')
end

for i=1:len1
    
    H(i,1)=ymax(i,1)-ymin(i,1); %Wave height
    xmean(i,1)=(xmax(i,1)+xmin(i,1))/2;
    Etac(i,1)=ymax(i,1); %Water surface elevation of the wave crest
    Etat(i,1)=ymin(i,1); %Water surface elevation of the wave trough
    
end
        
%--------------------------------------------------------------------------
%Calculating wave period 

for i=1:length(H(:,1))
    for j=1:length(xupcross(:,1))-1
        
        if xupcross(j,1)<xmean(i,1) & xupcross(j+1,1)>xmean(i,1)
            T(i,1)=xupcross(j+1,1)-xupcross(j,1);
        end
        
    end
end

%--------------------------------------------------------------------------
%Calculating wave properties

numberofwaves=length(H(:,1));
Etarms=sqrt(sum(Eta.^2)/N);

[Hsort,HsortIndex]=sort(H,'descend');
HTop3rdIndx=round(1/3*length(Hsort(:,1))); %Index where top 3rd wave height starts
Hs=mean(Hsort(1:HTop3rdIndx,1)); %Zero-crossing significant wave height

Tsort=T(HsortIndex);
TTop3rdIndx=round(1/3*length(Tsort(:,1))); %Index where top 3rd wave period starts
Ts=mean(Tsort(1:TTop3rdIndx,1)); %Zero-crossing significant wave period

Hrms=sqrt(sum(H.^2)/length(H(:,1))); %Root Mean Square Wave Height
Hz=mean(H(:,1)); %Zero-crossing mean wave height
Tz=mean(T(:,1)); %Zero-crossing mean wave period

%Hs=sqrt(2)*sqrt(sum(H.^2)/length(H(:,1))); %Significant wave height, Hs=(2^0.5)*Hrms
%Hz=sum(H)/numberofwaves; %Zero-crossing mean wave height
%Tz=duration/numberofwaves; %Zero-crossing mean wave period

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1

    val=[Hs Ts Hz Tz Hrms];
    name={'Hs','Ts','Hz','Tz','Hrms'};
    for i=1:length(val)
        fprintf('%14s   %g\n',name{i},val(i));
    end

    %Plotting data time series
    plot(t,Eta)
    hold on
    %plot(t,UDetrended,':')

    scatter(xmax,ymax,'*')
    scatter(xmin,ymin,'o')
    
    scatter(xupcross,yupcross,'^')
    scatter(xdowncross,ydowncross,'v')

    for i=1:length(H(:,1))

        plot([xmean(i,1),xmean(i,1)],[(ymax(i,1)+ymin(i,1))/2-H(i,1)/2,(ymax(i,1)+ymin(i,1))/2+H(i,1)/2],':m')

    end

    xlabel('Time (s)')
    ylabel('Water Surface Elevatio (m)')
    legend('Time Series','Crest','Trough','Up Crossing','Down Crossing','Wave Height')
    %legend('Time Series','Time Series without Kp','Crest','Trough','Up Crossing','Down Crossing','Wave Height')

end

%--------------------------------------------------------------------------
