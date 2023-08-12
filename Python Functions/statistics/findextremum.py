def findextremum(x, y, winlen=3, dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2018-05-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.findextremum
    ========================

    .. code:: python

        xmin, ymin, xmax, ymax = scientimate.findextremum(x, y, winlen=3, dispout='no')

    Description
    -----------

    Find local extremum (minimum and maximum) in data

    Inputs
    ------

    x
        x data
    y
        y data
    winlen=3
        | Window length, defines a number of points in sliding window used for defining maximum and minimum
        | Example: winlen=5 means two points on each side of each data point is used in calculation
        | Using a larger value for winlen makes it less sensitive 
        | winlen should be an odd number equal or larger than 3
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    xmin
        x of minmum points
    ymin
        y of minmum points
    xmax
        x of maximum points
    ymax
        y of maximum points

    Examples
    --------

    .. code:: python

        import scientimate as sm
        import numpy as np

        x=np.linspace(0,30,1000)
        y=2*np.exp(-0.1*2*np.pi/5*x)*np.sin(np.sqrt(1-0.1**2)*2*np.pi/5*x)
        xmin,ymin,xmax,ymax=sm.findextremum(x,y,3,'yes')

        x=np.linspace(0,50,1000)
        y=np.sin(x)+0.1*np.random.rand(1000)
        xmin,ymin,xmax,ymax=sm.findextremum(x,y,15,'yes')

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
    import scipy as sp
    if dispout=='yes':
        import matplotlib.pyplot as plt 

    #--------------------------------------------------------------------------
    #Convert inputs to numpy array

    #Changing type to numpy array
    def type2numpy(variable):
        if np.size(variable)==1:
            if ((type(variable) is list) or (type(variable) is np.ndarray)):
                variable=np.array(variable)
            else:
                variable=np.array([variable])
        elif np.size(variable)>1:
            if (type(variable).__module__)!='numpy':
                variable=np.array(variable) 
        return variable
    
    x=type2numpy(x)
    y=type2numpy(y)

    #--------------------------------------------------------------------------
    #Checking initial inputs

    #Length of the input data
    L=len(x)
    
    #Checking winlen be smaller than the length of input data
    if winlen>=L:
        winlen=L-1 #Subtracting 1 guarantees winlen stays smaller than L if L is even number

    
    #Checking winlen be larger than 3
    if winlen<3:
        winlen=3

    
    #Checking if winlen is an odd number
    if (winlen%2)==0:
        winlen=winlen+1

    
    #Defining a number of points on each side of each point for defining maximum and minimum
    npoints=int(winlen/2)

    #--------------------------------------------------------------------------
    #Locating local minimum and maximum points

    #Pre-assigning variables
    xmin1=np.zeros(L)
    ymin1=np.zeros(L)
    xmax1=np.zeros(L)
    ymax1=np.zeros(L)
    
    #Locating local minimum and maximum points for each sliding window length
    for i in range(0,npoints,1):
        xmin1[i]=x[i]
        ymin1[i]=np.min(y[0:winlen])
        xmax1[i]=x[i]
        ymax1[i]=np.max(y[0:winlen])

    
    for i in range(npoints,L-npoints,1):
        xmin1[i]=x[i]
        ymin1[i]=np.min(y[i-npoints:i+npoints+1])
        xmax1[i]=x[i]
        ymax1[i]=np.max(y[i-npoints:i+npoints+1])

    
    for i in range(L-npoints,L,1):
        xmin1[i]=x[i]
        ymin1[i]=np.min(y[L-winlen:])
        xmax1[i]=x[i]
        ymax1[i]=np.max(y[L-winlen:])

    
    #Locating local minimum and maximum points
    xmin=[]
    ymin=[]
    xmax=[]
    ymax=[]
    for i in range(0,L,1):
        
        if ymin1[i]==y[i]:
            xmin.append(xmin1[i])
            ymin.append(ymin1[i])

        if ymax1[i]==y[i]:
            xmax.append(xmax1[i])
            ymax.append(ymax1[i])

    #xmin=np.array(xmin)
    #ymin=np.array(ymin)
    #xmax=np.array(xmax)
    #ymax=np.array(ymax)

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
    
        #Plotting data
        plt.plot(x,y,label='Data')
        plt.scatter(xmax,ymax,c='r',marker='*',label='Maximum')
        plt.scatter(xmin,ymin,c='m',label='Minimum')
        plt.xlabel('x data')
        plt.ylabel('y data')
        plt.legend()
    
    
    #--------------------------------------------------------------------------
    #Outputs
    return xmin, ymin, xmax, ymax

    #--------------------------------------------------------------------------
