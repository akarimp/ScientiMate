def hurricanetranslationvel(xCenter, yCenter, dt=6*3600, distCalcMethod='gc', CalcMethod='backward', dispout='no'):
    """
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

    scientimate.hurricanetranslationvel
    ===================================

    .. code:: python

        Vt, VtAzmdir, VtTrigdir, distxy = scientimate.hurricanetranslationvel(xCenter, yCenter, dt=6*3600, distCalcMethod='gc', CalcMethod='backward', dispout='no')

    Description
    -----------

    Calculate hurricane center translational (forward motion) velocity

    Inputs
    ------

    xCenter
        x (longitude) of hurricane center (track)
    yCenter
        y (latitude) of hurricane center (track)
    dt=6*3600
        | Time interval between pressure data points in (s)
        | National Hurricane Center reports data every 6 hours 
    distCalcMethod='gc'
        | Distance calculation method 
        | 'cart': Distances are calculated on cartesian coordinate
        | 'gc': Distances are calculated on Great Circle based on Vincenty formula, Vincenty (1975)
        | Earth radius coonsidered as mean earth radius=6371000 m
    CalcMethod='backward'
        | Calculation method 
        | 'forward': Calculate hurricane central pressure intensity change over time using forward difference method
        |     If CalcMethod='forward'; then last element is zero
        | 'backward': Calculate hurricane central pressure intensity change over time using backward difference method
        |     If CalcMethod='backward'; then first element is zero
        | 'central': Calculate hurricane central pressure intensity change over time using central difference method
        |     If CalcMethod='central'; then first and last elements are zero
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    Vt
        Hurricane central translational velocity in (m/s)
    VtAzmdir
        | Hurricane center velocity azimuth (bearing) direction in (Degree)
        | azimuth (bearing) direction which is measured clockwise from the north:
        | 0 (degree): toward North, 90 (degree): toward East, 180 (degree): toward South, 270 (degree): toward West 
    VtTrigdir
        Hurricane center velocity trigonometric direction in (Degree)
    distxy
        Distance between hurricane locations of hurricane center in (m)

    Examples
    --------

    .. code:: python

        import scientimate as sm

        #Longitude of Hurricane Katrine best track
        longtrack=[-75.1,-75.7,-76.2,-76.5,-76.9,-77.7,-78.4,-79.0,-79.6,-80.1,-80.3,-81.3,\
            -82.0,-82.6,-83.3,-84.0,-84.7,-85.3,-85.9,-86.7,-87.7,-88.6,-89.2,-89.6,\
            -89.6,-89.6,-89.6,-89.6,-89.1,-88.6,-88.0,-87.0,-85.3,-82.9]

        #Latitude of Hurricane Katrine best track
        lattrack=[23.1,23.4,23.8,24.5,25.4,26.0,26.1,26.2,26.2,26.0,25.9,25.4,\
            25.1,24.9,24.6,24.4,24.4,24.5,24.8,25.2,25.7,26.3,27.2,28.2,\
            29.3,29.5,30.2,31.1,32.6,34.1,35.6,37.0,38.6,40.1]

        Vt,VtAzmdir,VtTrigdir,distxy=sm.hurricanetranslationvel(longtrack,lattrack,6*3600,'gc','backward','yes')

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
    """

    #--------------------------------------------------------------------------
    #CODE
    #--------------------------------------------------------------------------
    #Import required packages

    import numpy as np
    import scipy as sp
    if dispout=='yes':
        import matplotlib.pyplot as plt 

    #--------------------------------------------------------------------------
    #Convert inputs to numpy array

    #Changing type to numpy array
    def type2numpy(variable):
        if type(variable) is not str:
            if np.size(variable)==1:
                if ((type(variable) is list) or (type(variable) is np.ndarray)):
                    variable=np.array(variable)
                else:
                    variable=np.array([variable])
            elif np.size(variable)>1:
                if (type(variable).__module__)!='numpy':
                    variable=np.array(variable) 
        return variable
    
    xCenter=type2numpy(xCenter)
    yCenter=type2numpy(yCenter)

    #--------------------------------------------------------------------------
    #Pre-assigning array

    distxy=np.zeros(len(xCenter)) #Pre-assigning array
    thetarad=np.zeros(len(xCenter)) #Pre-assigning array
    deltasigma=np.zeros(len(xCenter)) #Pre-assigning array
    azimuthrad=np.zeros(len(xCenter)) #Pre-assigning array

    #--------------------------------------------------------------------------
    #Calculating distance between hurricane locations

    #Calculating distance using cartesian formula
    if distCalcMethod=='cart':
    
        #Calculating distance between hurricane locations using forward difference method
        if CalcMethod=='forward':
    
            for i in range(0,len(xCenter)-1,1):
    
                #Calculation range: (1:end-1,1), last element is zero
                distxy[i]=np.sqrt((xCenter[i+1]-xCenter[i])**2+(yCenter[i+1]-yCenter[i])**2) #Calculating distance from (x,y) to (x(1),y(1))
    
                #Calculating angle of the line between start and end points
                thetarad[i]=np.arctan2(yCenter[i+1]-yCenter[i],xCenter[i+1]-xCenter[i]) #Angle in radian
    
    
        #Calculating distance between hurricane locations using backward difference method
        elif CalcMethod=='backward':
    
            for i in range(1,len(xCenter),1):
    
                #Calculation range: (2:end,1), first element is zero
                distxy[i]=np.sqrt((xCenter[i]-xCenter[i-1])**2+(yCenter[i]-yCenter[i-1])**2) #Calculating distance from (x,y) to (x(1),y(1))
    
                #Calculating angle of the line between start and end points
                thetarad[i]=np.arctan2(yCenter[i]-yCenter[i-1],xCenter[i]-xCenter[i-1]) #Angle in radian
        
    
        #Calculating distance between hurricane locations using centra difference method
        elif CalcMethod=='central':
    
            for i in range(1,len(xCenter)-1,1):
    
                #Calculation range: (2:end-1,1), first and last elements are zero
                distxy[i]=np.sqrt((xCenter[i+1]-xCenter[i-1])**2+(yCenter[i+1]-yCenter[i-1])**2) #Calculating distance from (x,y) to (x(1),y(1))
    
                #Calculating angle of the line between start and end points
                thetarad[i]=np.arctan2(yCenter[i+1]-yCenter[i-1],xCenter[i+1]-xCenter[i-1]) #Angle in radian
    
    
        theta=np.rad2deg(thetarad) #Angle in degree
    
        #Add 360 to all numbers to have them all positive
        #Use mod(360) to take care of the ones larger than 360
        theta=((theta+360)%360) 
    
    
    #Calculating distance using Vincenty formula
    elif distCalcMethod=='gc':
    
        #Calculating distance between hurricane locations using forward difference method
        if CalcMethod=='forward':
    
            for i in range(0,len(xCenter)-1,1):
    
                #Calculation range: (1:end-1,1), last element is zero
    
                #Converting to radian
                lat1rad=np.deg2rad(yCenter[i])
                lon1rad=np.deg2rad(xCenter[i])
                lat2rad=np.deg2rad(yCenter[i+1])
                lon2rad=np.deg2rad(xCenter[i+1])
    
                deltalatrad21=lat2rad-lat1rad
                deltalonrad21=lon2rad-lon1rad
    
                deltasigma[i]=np.arctan2(np.sqrt((np.cos(lat2rad)*np.sin(deltalonrad21))**2+(np.cos(lat1rad)*np.sin(lat2rad)-np.sin(lat1rad)*np.cos(lat2rad)*np.cos(deltalonrad21))**2),np.sin(lat1rad)*np.sin(lat2rad)+np.cos(lat1rad)*np.cos(lat2rad)*np.cos(deltalonrad21)) #Central angle
    
                #Calculating azimuth (bearing) between start and end of the line
                azimuthrad[i]=np.arctan2(np.sin(deltalonrad21)*np.cos(lat2rad),np.cos(lat1rad)*np.sin(lat2rad)-np.sin(lat1rad)*np.cos(lat2rad)*np.cos(deltalonrad21)) #Azimuth (bearing) in radian
    
    
        #Calculating distance between hurricane locations using backward difference method
        elif CalcMethod=='backward':
    
            for i in range(1,len(xCenter),1):
    
                #Calculation range: (2:end,1), first element is zero
    
                #Converting to radian
                lat1rad=np.deg2rad(yCenter[i-1])
                lon1rad=np.deg2rad(xCenter[i-1])
                lat2rad=np.deg2rad(yCenter[i])
                lon2rad=np.deg2rad(xCenter[i])
    
                deltalatrad21=lat2rad-lat1rad
                deltalonrad21=lon2rad-lon1rad
    
                deltasigma[i]=np.arctan2(np.sqrt((np.cos(lat2rad)*np.sin(deltalonrad21))**2+(np.cos(lat1rad)*np.sin(lat2rad)-np.sin(lat1rad)*np.cos(lat2rad)*np.cos(deltalonrad21))**2),np.sin(lat1rad)*np.sin(lat2rad)+np.cos(lat1rad)*np.cos(lat2rad)*np.cos(deltalonrad21)) #Central angle
    
                #Calculating azimuth (bearing) between start and end of the line
                azimuthrad[i]=np.arctan2(np.sin(deltalonrad21)*np.cos(lat2rad),np.cos(lat1rad)*np.sin(lat2rad)-np.sin(lat1rad)*np.cos(lat2rad)*np.cos(deltalonrad21)) #Azimuth (bearing) in radian

    
        #Calculating distance between hurricane locations using centra difference method
        elif CalcMethod=='central':
    
            for i in range(1,len(xCenter)-1,1):
    
                #Calculation range: (2:end-1,1), first and last elements are zero
    
                #Converting to radian
                lat1rad=np.deg2rad(yCenter[i-1])
                lon1rad=np.deg2rad(xCenter[i-1])
                lat2rad=np.deg2rad(yCenter[i+1])
                lon2rad=np.deg2rad(xCenter[i+1])
    
                deltalatrad21=lat2rad-lat1rad
                deltalonrad21=lon2rad-lon1rad
    
                deltasigma[i]=np.arctan2(np.sqrt((np.cos(lat2rad)*np.sin(deltalonrad21))**2+(np.cos(lat1rad)*np.sin(lat2rad)-np.sin(lat1rad)*np.cos(lat2rad)*np.cos(deltalonrad21))**2),np.sin(lat1rad)*np.sin(lat2rad)+np.cos(lat1rad)*np.cos(lat2rad)*np.cos(deltalonrad21)) #Central angle
    
                #Calculating azimuth (bearing) between start and end of the line
                azimuthrad[i]=np.arctan2(np.sin(deltalonrad21)*np.cos(lat2rad),np.cos(lat1rad)*np.sin(lat2rad)-np.sin(lat1rad)*np.cos(lat2rad)*np.cos(deltalonrad21)) #Azimuth (bearing) in radian
    
    
        REarth=6371000 #Earth radius in (m), mean earth radius=6371000 m
        arclen=REarth*deltasigma #Total distance of the line 
        distxy=arclen.copy()
    
        azimuth=np.rad2deg(azimuthrad) #Azimuth (bearing) in degree
    
        #Add 360 to all numbers to have them all positive
        #Use mod(360) to take care of the ones larger than 360
        azimuth=((azimuth+360)%360) 
    

    #--------------------------------------------------------------------------
    #Calculating hurricane central translational velocity

    #Calculating hurricane central translational velocity using forward difference method
    if CalcMethod=='forward':
    
        #Calculation range: (1:end-1,1), last element is zero
        Vt=distxy/dt #Central translational velocity
    
    #Calculating hurricane central translational velocity using backward difference method
    elif CalcMethod=='backward':
    
        #Calculation range: (2:end,1), first element is zero
        Vt=distxy/dt #Central translational velocity
    
    #Calculating hurricane central translational velocity using centra difference method
    elif CalcMethod=='central':
    
        #Calculation range: (2:end-1,1), first and last elements are zero
        Vt=distxy/(2*dt) #Central translational velocity
    

    #--------------------------------------------------------------------------
    #Calculating hurricane central translational direction

    #Calculating direction using cartesian formula
    if distCalcMethod=='cart':
    
        #Converting trigonometric direction to azimuth (bearing)
        VtAzmdir=90-theta #Hurricane center velocity azimuth (bearing) direction in (Degree)
    
        VtTrigdir=theta #Hurricane center velocity trigonometric direction in (Degree)
    
    #Calculating direction using Vincenty formula
    elif distCalcMethod=='gc':
    
        VtAzmdir=azimuth #Hurricane center velocity azimuth (bearing) direction in (Degree)
    
        #Converting azimuth (bearing) to trigonometric direction
        VtTrigdir=-azimuth+90 #Hurricane center velocity trigonometric direction in (Degree)
    
    
    #Add 360 to all numbers to have them all positive
    #Use mod(360) to take care of the ones larger than 360
    VtAzmdir=((VtAzmdir+360)%360) 
    VtTrigdir=((VtTrigdir+360)%360) 

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
        
        #Plotting data
        t=np.arange(0,len(xCenter),1)*dt/3600
    
        plt.subplot(2,1,1)
        plt.plot(t,Vt)
        
        plt.xlabel('Time (hr)')
        plt.ylabel('Translational Velocity (m/s)')
        
        plt.subplot(2,1,2)
        plt.plot(t,VtAzmdir)
        
        plt.xlabel('Time (hr)')
        plt.ylabel('Velocity Azimuth (Degree)')
    

    #--------------------------------------------------------------------------
    #Outputs
    return Vt, VtAzmdir, VtTrigdir, distxy

    #--------------------------------------------------------------------------
