def reckongc(lat1, lon1, arclen, azimuthdir, R=6371000, dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2017-07-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.reckongc
    ====================

    .. code:: python

        lat2, lon2 = scientimate.reckongc(lat1, lon1, arclen, azimuthdir, R=6371000, dispout='no')

    Description
    -----------

    Calculate end point (Latitude,Longitude) from start point (Latitude,Longitude) and distance and azimuth (bearing) using Great Circle

    Inputs
    ------

    lat1
        Latitude (y) of start point (first point) in (Degree)
    lon1
        Longitude (x) of start point (first point) in (Degree)
    arclen
        Total distance from start point to end point in (m)
    azimuthdir
        | Azimuth (bearing or compass direction) from start point to end point in (Degree)
        | 0 (degree): toward North, 90 (degree): toward East, 180 (degree): toward South, 270 (degree): toward West 
    R=6371000
        Earth radius in (m), mean earth radius=6371000 m
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    lat2
        Latitude (y) of end point (last point) in (Degree)
    lon2
        Longitude (x) of end point (last point) in (Degree)

    Examples
    --------

    .. code:: python

        import scientimate as sm

        lat1=29.5 #First point 
        lon1=-89.4 #First point 
        arclen=22239 #Arc length
        azimuthdir=0 #Azimuth
        lat2,lon2=sm.reckongc(lat1,lon1,arclen,azimuthdir)

        lat1=[29.5,29] #First point 
        lon1=[-89.4,-89] #First point 
        arclen=[22239,147410] #Arc length
        azimuthdir=[0,319.21] #Azimuth
        lat2,lon2=sm.reckongc(lat1,lon1,arclen,azimuthdir,6371000,'yes')

    References
    ----------

    | http://www.movable-type.co.uk/scripts/latlong.html
    | https://en.wikipedia.org/wiki/Great-circle_distance
    | https://en.wikipedia.org/wiki/Great-circle_navigation
    | http://www.codeguru.com/cpp/cpp/algorithms/article.php/c5115/Geographic-Distance-and-Azimuth-Calculations.htm

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
    
    lat1=type2numpy(lat1)
    lon1=type2numpy(lon1)
    arclen=type2numpy(arclen)
    azimuthdir=type2numpy(azimuthdir)

    #--------------------------------------------------------------------------
    #Generating (lon2,lat2) from (lon1,lat1) 

    #Converting to radian
    lat1rad=np.deg2rad(lat1)
    lon1rad=np.deg2rad(lon1)
    azimuthdirrad=np.deg2rad(azimuthdir)
    
    #Generating (lon2,lat2) from (lon1,lat1) 
    deltasigma=arclen/R #Central angle 
    lat2rad=np.arcsin(np.sin(lat1rad)*np.cos(deltasigma)+np.cos(lat1rad)*np.sin(deltasigma)*np.cos(azimuthdirrad))
    lon2rad=lon1rad+np.arctan2(np.sin(azimuthdirrad)*np.sin(deltasigma)*np.cos(lat1rad),np.cos(deltasigma)-np.sin(lat1rad)*np.sin(lat2rad))
    
    #Converting to degree
    lat2deg=np.rad2deg(lat2rad)
    lon2deg=np.rad2deg(lon2rad)

    #--------------------------------------------------------------------------
    #Assigning outputs

    lat2=lat2deg.copy()
    lon2=lon2deg.copy()

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':

        plt.scatter(lon1,lat1,label='Start Point')
        plt.scatter(lon2,lat2,label='End Point')
    
        plt.xlabel('Longitude (Degree)')
        plt.ylabel('Latitude (Degree)')
        plt.legend()

    #--------------------------------------------------------------------------
    #Outputs
    return lat2, lon2

    #--------------------------------------------------------------------------
