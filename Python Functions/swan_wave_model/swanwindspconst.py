def swanwindspconst(windvel, winddir, winddirtype='mete', windvelmin=0, savedata='no', outfilename='swanwind.wnd', outfilelocation=None):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2017-11-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.swanwindspconst
    ===========================

    .. code:: python

        swanwind = scientimate.swanwindspconst(windvel, winddir, winddirtype='mete', windvelmin=0, savedata='no', outfilename='swanwind.wnd', outfilelocation=None)

    Description
    -----------

    Generate SWAN wind file for spatially constant wind

    Inputs
    ------

    windvel
        | Wind velocity
        | If size of windvel>1, then it is considered as a time series
        | 1st element is 1st time step, 2nd element is 2nd time step, ...
    winddir
        Wind direction in (Degree)
    winddirtype='mete'
        | Define wind direction type
        | 'mete': meteorological wind direction 
        |     Meteorological direction represents a direction wind comes from and is measured counter-clockwise from the North
        |     0 (degree): from North, 90 (degree): from East, 180 (degree): from South, 270 (degree): from West
        | 'trig': trigonometric wind direction
    windvelmin=0
        Minimum allowed wind velocity
    savedata='no'
        | Define if save data in a file or not
        | 'no': does not save 
        | 'yes': save data as ascii file
    outfilename='swanwind.wnd'
        | Name of output file between ' ' mark, example: 'swanwind.wnd'
        | outfilename should have '.wnd' extension
    outfilelocation=pwd
        Location of output file between ' ' mark, example: 'C:\' in MATLAB, or 'C:/' in Python

    Outputs
    -------

    swanwind
        | Spatially constant wind data formated for SWAN
        | Note: Wind at each time step is assigned into 4 points,
        |     assuming the wind domain is defined by 4 points, one at each corner

    Examples
    --------

    .. code:: python

        import scientimate as sm

        windvel=[10.5,10.6,10.55] #Data for 3 time steps
        winddir=[30,32,28] #Data for 3 time steps
        winddirtype='mete'
        windvelmin=2.5
        savedata='no'
        outfilename='swanwind.wnd'
        outfilelocation=None
        swanwind=sm.swanwindspconst(windvel,winddir,winddirtype,windvelmin,savedata,outfilename,outfilelocation)

    References
    ----------

    Booij, N. R. R. C., Ris, R. C., & Holthuijsen, L. H. (1999). 
    A third‐generation wave model for coastal regions: 1. Model description and validation. 
    Journal of geophysical research: Oceans, 104(C4), 7649-7666.

    SWAN Team. (2007). S
    WAN user manual. 
    Delft University of Technology. The Netherlands.

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
    import os

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
    
    windvel=type2numpy(windvel)
    winddir=type2numpy(winddir)
    windvelmin=type2numpy(windvelmin)

    #--------------------------------------------------------------------------
    #Assign default values

    if outfilelocation is None: outfilelocation=os.getcwd()

    #--------------------------------------------------------------------------
    #Defining required functions

    def cart2pol(x,y):
        rho=np.sqrt(x**2+y**2)
        theta=np.arctan2(y,x)
        return theta,rho
    
    def pol2cart(theta,rho):
        x=rho*np.cos(theta)
        y=rho*np.sin(theta)
        return x,y

    #--------------------------------------------------------------------------
    #Converting wind velocity into x and y components

    #Checking minimum value of wind velocity
    if (windvelmin<0): windvelmin=0
    
    #Checking for the minimum wind velocity
    windvel[windvel<windvelmin]=windvelmin
    
    #Meteorological direction
    if winddirtype=='mete':
        #Convert wind data from wind source to wind destination
        #Converting meteorological direction to trigonometric direction
        winddirtrig=-winddir+270
    
        #Add 360 to all numbers to have them all positive
        #Use mod(360) to take care of the ones larger than 360
        winddirtrig=((winddirtrig+360)%360) 
    
    #Trigonometric direction
    elif winddirtype=='trig':
        winddirtrig=winddir.copy()
    
    
    #Converting wind velocity into x and y components
    # windvelx=windvel.*cosd(90-winddirtrig);
    # windvely=windvel.*sind(90-winddirtrig);
    winddirtrigrad=np.deg2rad(winddirtrig) #Concverting wind direction to Radian
    windvelx,windvely=pol2cart(winddirtrigrad,windvel) #Calculating X and Y components of wind velocity

    #--------------------------------------------------------------------------
    #Generating SWAN wind file (.wnd)

    #If windvelx=[U1;U2;U3]; and windvely=[V1;V2;V3];
    #Then SWAN wind file (.wnd) for this windvel with three time steps has a format of:
    # U1 U1
    # U1 U1
    # V1 V1
    # V1 V1
    # U2 U2
    # U2 U2
    # V2 V2
    # V2 V2
    # U3 U3
    # U3 U3
    # V3 V3
    # V3 V3
    
    #Creating the wind
    M=len(windvel)
    swanwind=np.zeros((4*M,2))
    for i in range(0,len(windvel),1):
    
        #Wind velocity at time step i, first method
        swanwind[4*i:4*i+2,0:2]=windvelx[i] #x wind velocity at time step i
        swanwind[4*i+2:4*i+4,0:2]=windvely[i] #y wind velocity at time step i
    
        #Wind velocity at time step i, second method
        #x wind velocity at time step i
        #swanwind[4*i,0]=windvelx[i];
        #swanwind[4*i,1]=windvelx[i];
        #swanwind[4*i+1,0]=windvelx[i];
        #swanwind[4*i+1,1]=windvelx[i];
    
        #y wind velocity at time step i
        #swanwind[4*i+2,0]=windvely[i];
        #swanwind[4*i+2,1]=windvely[i];
        #swanwind[4*i+3,0]=windvely[i];
        #swanwind[4*i+3,1]=windvely[i];
    
    #--------------------------------------------------------------------------
    #Saving data

    if savedata=='yes':
    
        #Changing directory to saving directory
        currentFolder=os.getcwd()
        os.chdir(outfilelocation)
    
        #Saving data
        np.savetxt(outfilename,swanwind,delimiter=' ')
    
        #Changing directory to working directory
        os.chdir(currentFolder)

    #--------------------------------------------------------------------------
    #Outputs
    return swanwind

    #--------------------------------------------------------------------------
