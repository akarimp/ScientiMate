def hurricanedpcpt(Pc, dt=6*3600, CalcMethod='backward', dispout='no'):
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

    scientimate.hurricanedpcpt
    ==========================

    .. code:: python

        dPcdt, dPcdthPahr = scientimate.hurricanedpcpt(Pc, dt=6*3600, CalcMethod='backward', dispout='no')

    Description
    -----------

    Calculate hurricane central pressure (Pc) intensity change over time (dPc/dt)

    Inputs
    ------

    Pc
        Hurricane central surface pressure in (Pa)
    dt=6*3600
        | Time interval between pressure data points in (s)
        | National Hurricane Center reports data every 6 hours 
    CalcMethod='backward'
        | Calculation method 
        | 'forward': Calculate hurricane central pressure intensity change over time using forward difference method
        |            If CalcMethod='forward'; then last element is zero
        | 'backward': Calculate hurricane central pressure intensity change over time using backward difference method
        |            If CalcMethod='backward'; then first element is zero
        | 'central': Calculate hurricane central pressure intensity change over time using central difference method
        |            If CalcMethod='central'; then first and last elements are zero
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    dPcdt
        Hurricane central pressure (Pc) intensity change (dPc/dt) in (Pa/s)
    dPcdthPahr
        Hurricane central pressure (Pc) intensity change (dPc/dt) in (hPa/hr)

    Examples
    --------

    .. code:: python

        import scientimate as sm

        #Hurricane Katrina centeral pressure (Pa)
        Pc=[100800,100700,100700,100600,100300,100000,99700,99400,98800,98400,98300,98700,\
        97900,96800,95900,95000,94200,94800,94100,93000,90900,90200,90500,91300,\
        92000,92300,92800,94800,96100,97800,98500,99000,99400,99600]

        dPcdt,dPcdthPahr=sm.hurricanedpcpt(Pc,6*3600,'backward','yes')

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
    
    Pc=type2numpy(Pc)

    #--------------------------------------------------------------------------
    #Calculating hurricane central pressure (Pc) intensity change over time (dPc/dt)

    dPcdt=np.zeros(len(Pc)) #Pre-assigning array
    
    #Calculating hurricane central pressure intensity change over time using forward difference method
    if CalcMethod=='forward':
    
        #Calculation range: (1:end-1,1), last element is zero
        dPcdt[0:-2]=(Pc[1:-1]-Pc[0:-2])/dt #Central pressure intensity change over time
    
    #Calculating hurricane central pressure intensity change over time using backward difference method
    elif CalcMethod=='backward':
    
        #Calculation range: (2:end,1), first element is zero
        dPcdt[1:-1]=(Pc[1:-1]-Pc[0:-2])/dt #Central pressure intensity change over time
    
    #Calculating hurricane central pressure intensity change over time using centra difference method
    elif CalcMethod=='central':
    
        #Calculation range: (2:end-1,1), first and last elements are zero
        dPcdt[1:-2]=(Pc[2:-1]-Pc[0:-3])/(2*dt) #Central pressure intensity change over time
    
    
    #Convert central pressure intensity change over time from (Pa/s) to (hPa/hr)
    dPcdthPahr=dPcdt*3600*1e-2

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
        
        #Plotting data
        t=np.arange(0,len(Pc),1)*dt/3600
        plt.plot(t,dPcdt)
        
        plt.xlabel('Time (hr)')
        plt.ylabel('dPc/dp (Pa/s)')
        

    #--------------------------------------------------------------------------
    #Outputs
    return dPcdt, dPcdthPahr

    #--------------------------------------------------------------------------
