def swanwaterlevelspconst(waterlevel, savedata='no', outfilename='swanwaterlevel.wl', outfilelocation=None):
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

    scientimate.swanwaterlevelspconst
    =================================

    .. code:: python

        swanwaterlevel = scientimate.swanwaterlevelspconst(waterlevel, savedata='no', outfilename='swanwaterlevel.wl', outfilelocation=None)

    Description
    -----------

    | Generate SWAN water level file for spatially constant water level
    | This function can be used for any other scalar variable as well

    Inputs
    ------

    waterlevel
        | Water level (or any scalar variable)
        | If size of waterlevel>1, then it is considered as a time series
        | 1st element is 1st time step, 2nd element is 2nd time step, ...
    savedata='no'
        | Define if save data in a file or not
        | 'no': does not save 
        | 'yes': save data as ascii file
    outfilename='swanwaterlevel.wl'
        | Name of output file between ' ' mark, example: 'swanwaterlevel.wl'
        | outfilename should have '.wl' extension
        | In case of using other scalar varable than water level, use proper name and extension
    outfilelocation=pwd
        Location of output file between ' ' mark, example: 'C:\' in MATLAB, or 'C:/' in Python

    Outputs
    -------

    swanwaterlevel
        | Spatially constant water level data (or any scalar variable) formated for SWAN
        | Note: Water level at each time step is assigned into 4 points, 
        |     assuming the water level domain is defined by 4 points, one at each corner

    Examples
    --------

    .. code:: python

        import scientimate as sm

        waterlevel=[0.5,0.6,0.55] #Data for 3 time steps
        savedata='no'
        outfilename='swanwaterlevel.wl'
        outfilelocation=None
        swanwaterlevel=sm.swanwaterlevelspconst(waterlevel,savedata,outfilename,outfilelocation)

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
    
    waterlevel=type2numpy(waterlevel)

    #--------------------------------------------------------------------------
    #Assign default values

    if outfilelocation is None: outfilelocation=os.getcwd()

    #--------------------------------------------------------------------------
    #Generating SWAN water level file (.wl)

    #If waterlevel=[wl1;wl2;wl3];
    #Then SWAN water level file (.wl) for this waterlevel with three time steps has a format of:
    # wl1 wl1
    # wl1 wl1
    # wl2 wl2
    # wl2 wl2
    # wl3 wl3
    # wl3 wl3
    
    #Creating the water level
    M=len(waterlevel)
    swanwaterlevel=np.zeros((2*M,2))
    for i in range(0,len(waterlevel),1):
        
        #Water level at time step i, first method
        swanwaterlevel[2*i:2*i+2,0:2]=waterlevel[i]
    
        #Water level at time step i, second method
        #swanwaterlevel[2*i,0]=waterlevel[i];
        #swanwaterlevel[2*i,1]=waterlevel[i];
        #swanwaterlevel[2*i+1,0]=waterlevel[i];
        #swanwaterlevel[2*i+1,1]=waterlevel[i];
        
    #--------------------------------------------------------------------------
    #Saving data

    if savedata=='yes':
    
        #Changing directory to saving directory
        currentFolder=os.getcwd()
        os.chdir(outfilelocation)
    
        #Saving data
        np.savetxt(outfilename,swanwaterlevel,delimiter=' ')
    
        #Changing directory to working directory
        os.chdir(currentFolder)
    
    #--------------------------------------------------------------------------
    #Outputs
    return swanwaterlevel

    #--------------------------------------------------------------------------
