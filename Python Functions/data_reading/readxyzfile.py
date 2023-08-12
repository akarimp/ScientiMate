def readxyzfile(xyzfilename, xyzfilelocation=None, zscale=1, domain='all', xmin=-180, xmax=180, ymin=-90, ymax=90, savedata='no', outfilename='xyzdata.csv', outfilelocation=None):
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

    scientimate.readxyzfile
    =======================

    .. code:: python

        x, y, z = scientimate.readxyzfile(xyzfilename, xyzfilelocation=None, zscale=1, domain='all', xmin=-180, xmax=180, ymin=-90, ymax=90, savedata='no', outfilename='xyzdata.csv', outfilelocation=None)

    Description
    -----------

    | Read and extract x (longitude), y (latitude) and z (elevation) data from ASCII gridded (tabular) xyz file
    | Use readdatafile function for more options

    Inputs
    ------

    xyzfilename
        | Name of xyz file between ' ' mark, example: 'xyzfile.xyz'
        | xyz file should be in form of 3 coloumn format
    xyzfilelocation=pwd
        Location of xyz file between ' ' mark, example: 'C:\'
    zscale=1
        Scale z (elevation) data by factor of zscale
    domain='all'
        | Define a domain to be extracted from data
        | 'all': all xyz data in input file are extracted 
        | 'domain': only data within a defined domain are extracted
    xmin=-180
        Minimum x (longitude) of domain to be extracted
    xmax=180
        Maximum x (longitude) of domain to be extracted
    ymin=-90
        Minimum y (latitude) of domain to be extracted
    ymax=90
        Maximum y (latitude) of domain to be extracted
    savedata='no'
        | Define if save xyz data in a file or not in outfilelocation folder
        | 'no': does not save
        | 'yes': save xyz data as csv file
    outfilename='xyzdata.csv'
        | Name of output file between ' ' mark, example: 'xyzdata.csv'
        | outfilename should have '.csv' extension
    outfilelocation=pwd
        Location of output file between ' ' mark, example: 'C:\'

    Outputs
    -------

    x
        x (longitude) data extracted from xyz file
    y
        y (latitude) data extracted from xyz file
    z
        z (elevation) data extracted from xyz file

    Examples
    --------

    .. code:: python

        import scientimate as sm

        xyzfilename='xyzfile.xyz' #e.g. xyzfilename='PersianGulf_ETOPO1.xyz'
        xyzfilelocation='C:/' #e.g. xyzfilelocation='C:/datafolder'
        x,y,z=sm.readxyzfile(xyzfilename,xyzfilelocation)

        xyzfilename='xyzfile.xyz' #e.g. xyzfilename='PersianGulf_ETOPO1.xyz'
        xyzfilelocation='C:/' #e.g. xyzfilelocation='C:/datafolder'
        x,y,z=sm.readxyzfile(xyzfilename,xyzfilelocation,1,'all',-180,180,-90,90,'no')

    References
    ----------

    Geospatial data

    * https://www.mathworks.com/help/map/finding-geospatial-data.html
    * https://maps.ngdc.noaa.gov/viewers/wcs-client/
    * https://www.ngdc.noaa.gov/mgg/global/global.html
    * https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/
    * https://www.ngdc.noaa.gov/mgg/image/2minrelief.html
    * https://www.ngdc.noaa.gov/mgg/coastal/crm.html
    * https://viewer.nationalmap.gov/launch/
    * https://earthexplorer.usgs.gov
    * http://www.shadedrelief.com/cleantopo2/index.html

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
    
    #x=type2numpy(x)

    #--------------------------------------------------------------------------
    #Assigning default values

    if xyzfilelocation is None: xyzfilelocation=os.getcwd()
    if outfilename is None: outfilename=os.getcwd()

    #--------------------------------------------------------------------------
    #Reading xyz file

    #Reading all data
    currentFolder=os.getcwd()
    os.chdir(xyzfilelocation)

    #Importing data
    with open(xyzfilename,'r') as xyzfile:
        
        #Reading data in xyzfile, first method
        #xyzdatalist=list(xyzfile)
        ##xyzdatalist=xyzfile.readlines()
        #xyzdata=[]
        #for i in range(0,len(xyzdatalist),1):
        #    #row=xyzdatalist[i].strip()
        #    row=xyzdatalist[i].split()
        #    xyzdata.append(row)
           
        #Reading data in xyzfile, second method
        #x,y,z=[],[],[]
        xyzdata=[]
        for line in xyzfile:
            row=line.split()
            #x.append(row[0])
            #y.append(row[1])
            #z.append(row[2])
            xyzdata.append(row)
        
        #Converting a list to numpy array
        xyzdata=np.array(xyzdata).astype('float')
    
    #Reading all xyz data in input file
    if domain=='all':
    
        x=xyzdata[:,0] #x (longitude) data in (Degree)
        y=xyzdata[:,1] #y (latitude) data in (Degree)
        z=zscale*xyzdata[:,2] #z (elevation) data in (m)
    
        xyzdata1=np.zeros((len(x),3))
        xyzdata1[:,0]=x  #x (longitude) data in (Degree)
        xyzdata1[:,1]=y  #y (latitude) data in (Degree)
        xyzdata1[:,2]=z #z (elevation) data in (m)
    
    #Reading data within a defined domain
    elif domain=='domain':
        
        #All data
        x1=xyzdata[:,0] #x (longitude) data in (Degree)
        y1=xyzdata[:,1] #y (latitude) data in (Degree)
        z1=zscale*xyzdata[:,2] #z (elevation) data in (m)
        
        #Finding a location data that are within [xmin:xmax] range
        xIndx=np.int_((np.nonzero((x1>=xmin) & (x1<=xmax)))[0])
        
        #Retaining data that are within [xmin:xmax] range
        x2=x1[xIndx]   #xdata
        y2=y1[xIndx]   #ydata
        z2=z1[xIndx] #zdata
        
        #Finding a location data that are within [ymin:ymax] range
        yIndx=np.int_((np.nonzero((y2>=ymin) & (y2<=ymax)))[0])
        
        #Retaining data that are within [ymin:ymax] range
        x=x2[yIndx]   #x data
        y=y2[yIndx]   #y data
        z=z2[yIndx] #z data
        
        xyzdata1=np.zeros((len(x),3))
        xyzdata1[:,0]=x  #x (longitude) data in (Degree)
        xyzdata1[:,1]=y  #y (latitude) data in (Degree)
        xyzdata1[:,2]=z #z (elevation) data in (m)

    #--------------------------------------------------------------------------
    #Saving xyz data

    if savedata=='yes':
        os.chdir(outfilelocation)
        #np.savetxt('xyzdata.csv',xyzdata1,delimiter=',')
        #np.savetxt('xyzdata.txt',xyzdata1)
        np.savetxt(outfilename,xyzdata1,delimiter=',')

    #--------------------------------------------------------------------------
    #Changing directory to working directory

    os.chdir(currentFolder)

    #--------------------------------------------------------------------------
    #Outputs
    return x, y, z

    #--------------------------------------------------------------------------
