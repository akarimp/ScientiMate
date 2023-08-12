def seqcolormap(dispout='no'):
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

    scientimate.seqcolormap
    =======================

    .. code:: python

        CmapAKLines, CmapAKSet1, CmapAKPaired, CmapAKExtraLines, CmapAKTab10, CBrewerSet1, CBrewerPaired, Tab10 = scientimate.seqcolormap(dispout='no')

    Description
    -----------

    Generate sequential colormap for drawing lines

    Inputs
    ------

    dispout='no'
        Define to display outputs or not ('yes': display, 'no': not display)

    Outputs
    -------

    CmapAKLines
        7 Colors similar to matlab Lines in RGB (between 0 and 1)
    CmapAKSet1
        9 Colors similar to colorbrewer2.org Set1 in RGB (between 0 and 1)
    CmapAKPaired
        12 Colors similar to colorbrewer2.org Paired in RGB (between 0 and 1)
    CmapAKExtraLines
        15 extra line colors similar to cbrewer qualitative Paired in RGB (between 0 and 1)
    CmapAKTab10
        10 Colors similar to classic Tableau 10 color palettes in RGB (between 0 and 1)
    CBrewerSet1
        9 Colors Set1 from Cynthia Brewer's ColorBrewer colorbrewer2.org in RGB (between 0 and 1)
    CBrewerPaired
        12 Colors Paired from Cynthia Brewer's ColorBrewer colorbrewer2.org in RGB (between 0 and 1)
    Tab10
        10 Colors from classic Tableau 10 color palettes in RGB (between 0 and 1)

    Examples
    --------

    .. code:: python

        import scientimate as sm
        CmapAKLines,CmapAKSet1,CmapAKPaired,CmapAKExtraLines,CmapAKTab10,CBrewerSet1,CBrewerPaired,Tab10=sm.seqcolormap('yes')

    References
    ----------

    Colormap

    * http://colorbrewer2.org
    * http://matplotlib.org/cmocean/
    * https://matplotlib.org/users/colormaps.html
    * http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml
    * https://www.giss.nasa.gov/tools/panoply/colorbars/
    * http://jdherman.github.io/colormap/

    Source

    * http://www.beautycolorcode.com  sample: http://www.beautycolorcode.com/0072bd  (RGB and hex code)
    * http://www.htmlcsscolor.com/
    * http://www.colorhexa.com/
    * http://colorbrewer2.org/
    * http://onlinehelp.tableau.com/current/pro/desktop/en-us/help.htm#formatting_create_custom_colors.html

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
    #Loading a function or module

    #os.chdir('function or module directory')
    #import ModuleName
    #ModuleName.FunctionName(Par)
    #Var=ModuleName.FunctionName(Par)
    #Var1,Var2,...,VarN=ModuleName.FunctionName(Par)
    #Var1,_,...,_=ModuleName.FunctionName(Par)

    #or

    #from ModuleName import FunctionName
    #FunctionName(Par)
    #Var=FunctionName(Par)
    #Var1,Var2,...,VarN=FunctionName(Par)
    #Var1,_,...,_=FunctionName(Par)

    #Reloading module or function if it is updated
    #from imp import reload
    #reload(ModuleName)

    #--------------------------------------------------------------------------
    #Using descrete Matplotlib colormap

    #import matplotlib as mpl 

    #Colormap to array
    #pltcmap=mpl.cm.get_cmap(zcolormap)
    #for i in range(pltcmap.N):
    #    clmap=pltcmap(i)[:3]

    #Returns tulip contains 256 RGB
    #col1=mpl.cm.get_cmap('Set1')
    #col1=mpl.cm.jet

    #NCol=8 #Total number of output colormap
    #col1=mpl.cm.viridis(np.arange(256))
    #col1=col1[:,0:3]
    #col2=np.ones((NCol,3))
    #col2[:,0]=np.interp(np.linspace(0,255,NCol),np.arange(256),col1[:,0])
    #col2[:,1]=np.interp(np.linspace(0,255,NCol),np.arange(256),col1[:,1])
    #col2[:,2]=np.interp(np.linspace(0,255,NCol),np.arange(256),col1[:,2])

    #Return color array contains Ncol
    #Ncol=11 #Total number of output colormap
    #CM=mpl.cm.get_cmap('Paired') #Colormap tulip from 0 to 255, example CM(0) or CM(255)
    #col=CM(np.linspace(0,1,Ncol))
    #col=CM(np.int_(np.linspace(0,255,Ncol)))

    #bounds = np.linspace(np.min(FetchhatBS1/hhatBS1),np.max(FetchhatBS1/hhatBS1),Ncol)
    #norm = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=256)
    #fig1=plt.scatter(hhatBS1,fphatBS1,50,FetchhatBS1/hhatBS1,cmap=CM,norm=norm)
    #plt.colorbar()

    #--------------------------------------------------------------------------
    #Import required packages

    import numpy as np
    if dispout=='yes':
        import matplotlib.pyplot as plt 
        from cycler import cycler

    #--------------------------------------------------------------------------
    #Colormaps

    #7 Colors similar to matlab Lines
    CmapAKLines=np.array([[0,114,187],[226,88,34],[225,173,33],[128,55,144],[122,172,33],[102,183,225],[159,29,53]])/255
    
    #9 Colors similar to to colorbrewer2.org Set1
    CmapAKSet1=np.array([[206,22,32],[61,133,184],[75,163,81],[159,95,159],[255,140,0],[253,227,54],[161,82,38],[247,127,190],[169,157,157]])/255
    
    #12 Colors similar to to colorbrewer2.org Paired
    CmapAKPaired=np.array([[164,210,224],[49,110,160],[165,215,133],[34,139,34],[255,152,137],[206,22,32],[255,200,120],[255,140,0],[202,180,212],[107,63,160],[255,255,153],[177,89,47]])/255
    
    #15 extra line colors similar to cbrewer qualitative Paired
    CmapAKExtraLines=np.array([[164,210,224],[0,153,255],[0,127,255],[0,102,204],[0,0,255],[165,215,133],[75,163,81],[255,153,51],[255,102,0],[226,88,34],[177,89,47],[224,176,255],[202,180,212],[255,153,153],[230,32,32]])/255
    
    #10 Colors similar to classic Tableau 10 color palettes
    CmapAKTab10=np.array([[49,110,160],[255,117,24],[34,139,34],[227,38,54],[148,112,196],[139,80,75],[247,127,190],[128,128,128],[183,198,26],[3,180,200]])/255

    #Cynthia Brewer's ColorBrewer colorbrewer2.org (should be cited if used)
    CBrewerSet1=np.array([[228,26,28],[55,126,184],[77,175,74],[152,78,163],[255,127,0],[255,255,51],[166,86,40],[247,129,191],[153,153,153]])/255
    CBrewerPaired=np.array([[166,206,227],[31,120,180],[178,223,138],[51,160,44],[251,154,153],[227,26,28],[253,191,111],[255,127,0],[202,178,214],[106,61,154],[255,255,153],[177,89,40]])/255

    #10 Colors from classic Tableau 10 color palettes
    Tab10=np.array([[31,119,180],[255,127,14],[44,160,44],[214,39,40],[148,103,189],[140,86,75],[227,119,194],[127,127,127],[188,189,34],[23,190,207]])/255

    # CmapAKLines=[0,114,187;255,128,0;225,173,33;128,55,144;124,159,47;102,183,225;159,29,53]./255;
    # CmapAKExtraLines=[164,210,224;0,153,255;0,127,255;0,102,204;0,0,255;165,215,133;75,163,81;255,153,51;255,102,0;226,88,34;177,89,47;224,176,255;202,180,212;255,153,153;230,32,32]./255;

    #--------------------------------------------------------------------------
    #Displaying results

    if dispout=='yes':

        plt.figure(1)
        #plt.subplot(3,3,1)
        plt.subplot(2,2,1)
        N=7 # Number of color in color set
        X = np.linspace(0,np.pi*3,1000)
        Y = np.zeros((1000,N))
        for i in range(0,N):
            Y[:,i]=np.sin(X+2*i*np.pi/N)
        #mpl.rcParams['axes.color_cycle'] = CmapAKLines
        plt.rc('axes', prop_cycle=(cycler('color', CmapAKLines)))
        plt.plot(X,Y,linewidth=5,label='Color')
        plt.ylim(-1.1,1.1)
        plt.legend(loc='upper right',ncol=2)
        plt.title('CmapAKLines similar to matlab lines')

        #plt.figure(2)
        #plt.subplot(3,3,2)
        plt.subplot(2,2,2)
        N=9 # Number of color in color set
        X = np.linspace(0,np.pi*3,1000)
        Y = np.zeros((1000,N))
        for i in range(0,N):
            Y[:,i]=np.sin(X+2*i*np.pi/N)
        plt.rc('axes', prop_cycle=(cycler('color', CmapAKSet1)))
        plt.plot(X,Y,linewidth=5,label='Color')
        plt.ylim(-1.1,1.1)
        plt.legend(loc='upper right',ncol=2)
        plt.title('CmapAKSet1 similar to colorbrewer Set1')

        #plt.figure(3)
        #plt.subplot(3,3,3)
        plt.subplot(2,2,3)
        N=12 # Number of color in color set
        X = np.linspace(0,np.pi*3,1000)
        Y = np.zeros((1000,N))
        for i in range(0,N):
            Y[:,i]=np.sin(X+2*i*np.pi/N)
        plt.rc('axes', prop_cycle=(cycler('color', CmapAKPaired)))
        plt.plot(X,Y,linewidth=5,label='Color')
        plt.ylim(-1.1,1.1)
        plt.legend(loc='upper right',ncol=2)
        plt.title('CmapAKPaired similar to colorbrewer Paired')

        #plt.figure(4)
        #plt.subplot(3,3,4)
        plt.subplot(2,2,4)
        N=15 # Number of color in color set
        X = np.linspace(0,np.pi*3,1000)
        Y = np.zeros((1000,N))
        for i in range(0,N):
            Y[:,i]=np.sin(X+2*i*np.pi/N)
        plt.rc('axes', prop_cycle=(cycler('color', CmapAKExtraLines)))
        plt.plot(X,Y,linewidth=5,label='Color')
        plt.ylim(-1.1,1.1)
        plt.legend(loc='upper right',ncol=2)
        plt.title('CmapAKExtraLines')

        plt.figure(2)
        #plt.figure(5)
        #plt.subplot(3,3,5)
        plt.subplot(2,2,1)
        N=10 # Number of color in color set
        X = np.linspace(0,np.pi*3,1000)
        Y = np.zeros((1000,N))
        for i in range(0,N):
            Y[:,i]=np.sin(X+2*i*np.pi/N)
        plt.rc('axes', prop_cycle=(cycler('color', CmapAKTab10)))
        plt.plot(X,Y,linewidth=5,label='Color')
        plt.ylim(-1.1,1.1)
        plt.legend(loc='upper right',ncol=2)
        plt.title('CmapAKTab10 similar to classic Tableau 10 color palettes')

        #plt.figure(6)
        #plt.subplot(3,3,6)
        plt.subplot(2,2,2)
        N=9 # Number of color in color set
        X = np.linspace(0,np.pi*3,1000)
        Y = np.zeros((1000,N))
        for i in range(0,N):
            Y[:,i]=np.sin(X+2*i*np.pi/N)
        plt.rc('axes', prop_cycle=(cycler('color', CBrewerSet1)))
        plt.plot(X,Y,linewidth=5,label='Color')
        plt.ylim(-1.1,1.1)
        plt.legend(loc='upper right',ncol=2)
        plt.title('Set1 from colorbrewer2.org')

        #plt.figure(7)
        #plt.subplot(3,3,7)
        plt.subplot(2,2,3)
        N=12 # Number of color in color set
        X = np.linspace(0,np.pi*3,1000)
        Y = np.zeros((1000,N))
        for i in range(0,N):
            Y[:,i]=np.sin(X+2*i*np.pi/N)
        plt.rc('axes', prop_cycle=(cycler('color', CBrewerPaired)))
        plt.plot(X,Y,linewidth=5,label='Color')
        plt.ylim(-1.1,1.1)
        plt.legend(loc='upper right',ncol=2)
        plt.title('Paired from colorbrewer2.org')

        #plt.figure(8)
        #plt.subplot(3,3,8)
        plt.subplot(2,2,4)
        N=10 # Number of color in color set
        X = np.linspace(0,np.pi*3,1000)
        Y = np.zeros((1000,N))
        for i in range(0,N):
            Y[:,i]=np.sin(X+2*i*np.pi/N)
        plt.rc('axes', prop_cycle=(cycler('color', Tab10)))
        plt.plot(X,Y,linewidth=5,label='Color')
        plt.ylim(-1.1,1.1)
        plt.legend(loc='upper right',ncol=2)
        plt.title('Tab10 from classic Tableau 10 color palettes')

    #--------------------------------------------------------------------------
    #OUtputs
    return CmapAKLines, CmapAKSet1, CmapAKPaired, CmapAKExtraLines, CmapAKTab10, CBrewerSet1, CBrewerPaired, Tab10

    #--------------------------------------------------------------------------
