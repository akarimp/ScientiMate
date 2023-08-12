def gencolormap(cmapcolors=[[255,255,255],[0,0,0]], ncolor=256, dispout='no'):
    """
    .. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
    .. +                                                                        +
    .. + ScientiMate                                                            +
    .. + Earth-Science Data Analysis Library                                    +
    .. +                                                                        +
    .. + Developed by: Arash Karimpour                                          +
    .. + Contact     : www.arashkarimpour.com                                   +
    .. + Developed/Updated (yyyy-mm-dd): 2017-08-01                             +
    .. +                                                                        +
    .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    scientimate.gencolormap
    =======================

    .. code:: python

        cmap_ncolor, cmap_ncolor_plt = scientimate.gencolormap(cmapcolors=[[255,255,255],[0,0,0]], ncolor=256, dispout='no')

    Description
    -----------

    Generate a colormap from input colors

    Inputs
    ------

    cmapcolors=[[255,255,255],[0,0,0]]
        | Colors that are used to generate colormap
        | Colors should be defined as [M*3] array in RGB color format
        | At least two colors should be defined, i.e. M should be equal or larger than 2
        | All values should be between 0 and 255
        | Any available colormap name such as 'cool', 'winter', etc also can be used
    ncolor=256
        Number of colors in generated colormap
    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    cmap_ncolor
        | Colormap with ncolor number of colors in RGB color format between 0 and 1
        | To convert 0-1 scale to 0-255 scale, multiply cmap_ncolor values by 255
    cmap_ncolor_plt
        Colormap with ncolor number of colors in Matplotlib format

    Examples
    --------

    .. code:: python

        import scientimate as sm

        #Blue-Red
        cmap_ncolor,cmap_ncolor_plt=sm.gencolormap([[41, 128, 185],[192, 57, 43]],256,'yes')

        #Blue-White-Red
        cmap_ncolor,cmap_ncolor_plt=sm.gencolormap([[41, 128, 185],[255,255,255],[192, 57, 43]],256,'yes')

        #Blue-White-Green
        cmap_ncolor,cmap_ncolor_plt=sm.gencolormap([[33, 150, 243],[255,255,255],[76, 175, 80]],256,'yes')

        #Red-Brown sequential
        cmapcolors=[[192,57,43],[155,89,182],[41,128,185],[39,174,96],[241,196,15],[211,84,0]]
        cmap_ncolor,cmap_ncolor_plt=sm.gencolormap(cmapcolors,10,'yes')

        #cool colormap
        cmap_ncolor,cmap_ncolor_plt=sm.gencolormap('cool',256,'yes')

    References
    ----------

    Colormap

    * http://colorbrewer2.org
    * http://matplotlib.org/cmocean/
    * https://matplotlib.org/users/colormaps.html
    * http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml
    * https://www.giss.nasa.gov/tools/panoply/colorbars/
    * http://jdherman.github.io/colormap/

    Color

    * http://htmlcolorcodes.com

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
    import matplotlib as mpl 
    from matplotlib import colors 

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

    cmapcolors=type2numpy(cmapcolors)

    #--------------------------------------------------------------------------
    #Colormaps

    #User input colormap
    if ((type(cmapcolors) is list) or (type(cmapcolors) is np.ndarray)):
    
        zclmap=cmapcolors/255 #Converting RGB to values between 0 and 1
        zclmap[zclmap<0]=0
        zclmap[zclmap>1]=1
    
    #Pre-defined colormap such as 'cool', 'winter', etc
    else:
    
        pltcmap=mpl.cm.get_cmap(cmapcolors)
        zclmap=[]
        for i in range(pltcmap.N):
            zclmap.append(pltcmap(i)[0:3])
        zclmap=np.array(zclmap).astype('float') #Converting a list to numpy array
    

    #--------------------------------------------------------------------------
    #Interpolating a colormap to include ncolor colors

    cmap_ncolor=np.zeros((ncolor,3))
    for i in range(0,3,1):
        cmap_ncolor[:,i]=np.interp(np.linspace(1,ncolor,ncolor),np.linspace(1,ncolor,len(zclmap[:,0])),zclmap[:,i])

    #Generating matplotlib colormap 
    #Use matplotlib.colors.ListedColormap or matplotlib.colors.LinearSegmentedColormap
    #cmap_ncolor_plt=colors.LinearSegmentedColormap.from_list(name='my_colormap',colors=cmap_ncolor,N=ncolor)
    cmap_ncolor_plt=colors.ListedColormap(cmap_ncolor,name='my_colormap',N=None)

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':
    
        xgrid,ygrid=np.meshgrid(np.linspace(1,ncolor,ncolor),np.linspace(1,ncolor,ncolor))
        zgrid=ygrid.copy()
    
        plt.imshow(zgrid,extent=[np.min(xgrid),np.max(xgrid),np.min(ygrid),np.max(ygrid)],cmap=cmap_ncolor_plt,aspect='auto',origin='lower')
    
        #Plotting colorbar
        plt.colorbar()

    #--------------------------------------------------------------------------
    #Outputs
    return cmap_ncolor, cmap_ncolor_plt

    #--------------------------------------------------------------------------
