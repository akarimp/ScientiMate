def linearwavegenerator(amin, amax, Tmin, Tmax, Phimin=0, Phimax=2*3.1416, fs=32, duration=10, NoOfWave=2, dispout='no'):
    """
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

    scientimate.linearwavegenerator
    ===============================

    .. code:: python

        Eta, t, Etaij = scientimate.linearwavegenerator(amin, amax, Tmin, Tmax, Phimin=0, Phimax=2*3.1416, fs=32, duration=10, NoOfWave=2, dispout='no')

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
    Phimin=0
        Min Phase (radian)
    Phimax=2*pi
        Max Phase (radian) 
    fs=32
        Sample generation frequency (Hz), number of data points in one second
    duration=10
        Duration time that data will be generated in (s)
    NoOfWave=2
        Number of waves to be combined with each other
    dispout='no'
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

    .. code:: python

        import scientimate as sm
        import numpy as np

        Eta,t,Etaij=sm.linearwavegenerator(0.2,0.4,1,3,0,2*np.pi,32,10,2,'yes')

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
    """

    #--------------------------------------------------------------------------
    #CODE
    #--------------------------------------------------------------------------
    #Import required packages

    import numpy as np
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
    
    #x=type2numpy(x)

    #--------------------------------------------------------------------------

    sample=fs*duration #number of sample in input file
    dt=1/fs
    t=np.linspace(dt,duration,sample) #time
    a=np.linspace(amin,amax,NoOfWave)
    T=np.linspace(Tmin,Tmax,NoOfWave)
    Phi=np.linspace(Phimin,Phimax,NoOfWave)
    
    w=2*np.pi/T #Wave angular frequency
    
    Etaij=np.zeros((len(t),NoOfWave)) #Pre-assigning array to make program run faster
    for i in range(0,NoOfWave,1):
        Etaij[:,i]=a[i]*np.cos(-w[i]*t+Phi[i])

    Eta=np.sum(Etaij,axis=1)

#    Etaij=np.empty([len(t),0])
#    for i in range(0,NoOfWave,1):
#        Etaij=np.column_stack((Etaij,(a[i])*np.cos(-w[i]*t+Phi[i])))
#     
#    Eta=np.empty(len(t))
#    for i in range(0,len(t)):
#        Eta[i]=np.sum(Etaij[i,0:])

    # -------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
        
        plt.subplot(2,1,1)
        for i in range(0,NoOfWave,1):
            plt.plot(t,Etaij[:,i])
            #plt.hold='True'
        plt.title('Generated Waves')
        plt.xlabel('Time (s)')
        plt.ylabel('Water Level (m)')
        plt.xlim(t[0],t[-1])
         
        plt.subplot(2,1,2)
        plt.plot(t,Eta,label='Water Level (m)')
        #plt.hold='True'
        plt.title('Generated Combined Wave')
        plt.xlabel('Time (s)')
        plt.ylabel('Water Level (m)')
        plt.xlim(t[0],t[-1])        
 
    # -------------------------------------------------------------------------
    #Outputs
    return Eta, t, Etaij

    # -------------------------------------------------------------------------
