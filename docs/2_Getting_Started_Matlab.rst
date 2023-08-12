Getting Started (MATLAB)
========================

In order to use MATLAB version of ScientiMate library, first, one of the MATLAB or GNU Octave programming language along with required toolboxes should be installed (Refer to required packages for MATLAB/GNU Octave). 


Installing
----------

To use MATLAB version of ScientiMate library:

* Install MATLAB or GNU Octave

    * MATLAB: https://www.mathworks.com
    * GNU Octave: https://octave.org

* Download ScientiMate:

    * Version 1.1 (GitHub): https://github.com/akarimp/ScientiMate/releases/download/2.0/scientimate.zip
    
* Unzip ScientiMate in any location you choose such as "C:\\"
* Add ScientiMate folder to MATLAB or GNU Octave path


Add ScientiMate folder to MATLAB or GNU Octave path
---------------------------------------------------

You may access ScientiMate by copying ScientiMate files and its sub-folders to your desire working directory and then use them there.
However, a better option is to add ScientiMate folder to MATLAB or GNU Octave path. By doing that, you always have access to ScientiMate from any working directory.
Remember, you need to add ScientiMate to path only once.

To add ScientiMate folder to MATLAB or GNU Octave path, or remove it from the path, you may use one the following options:

* Add or remove the folder by using provided ``add_scientimate_to_path.m`` or ``remove_scientimate_from_path.m`` files
* Add or remove the folder manually in the Command Window

**Add ScientiMate folder to MATLAB or GNU Octave path by using add_scientimate_to_path.m**

To add ScientiMate to MATLAB or GNU Octave path:

* Open MATLAB or GNU Octave
* Change a current folder (working directory) to a folder that contains ScientiMate files, for example "C:\\scientimate", in MATLAB or GNU Octave.
* Run a file named ``add_scientimate_to_path.m`` in MATLAB or GNU Octave to add ScientiMate folder to MATLAB or GNU Octave path.

To remove ScientiMate from MATLAB or GNU Octave path:

* Open MATLAB or GNU Octave
* Change a current folder (working directory) to a folder that contains ScientiMate files, for example "C:\\scientimate", in MATLAB or GNU Octave.
* Run a file named ``remove_scientimate_from_path.m`` in MATLAB or GNU Octave to remove ScientiMate folder from MATLAB or GNU Octave path.

**Add ScientiMate folder to MATLAB or GNU Octave path manually in the Command Window**

For example, if ScientiMate files are in "C:\\scientimate" folder, then:

To add ScientiMate to MATLAB or GNU Octave path, run following commands in the Command Window:

.. code:: matlab

    ScientiMatePath = genpath('C:\scientimate'); %Generating path for ScientiMate folder and its sub-folders
    addpath(ScientiMatePath); %Add ScientiMate folder to path

To remove ScientiMate from MATLAB or GNU Octave path, run following commands in the Command Window:

.. code:: matlab

    ScientiMatePath = genpath('C:\scientimate'); %Generating path for ScientiMate folder and its sub-folders
    rmpath(ScientiMatePath); %Remove ScientiMate folder from path
    %restoredefaultpath; %Restore path to factory-installed state


Operating System
----------------

This code can be run on Microsoft Windows, Mac and Linux. However, make sure any given path is compatible with a running operating system. In particular, "\\" is used in Windows path, while "/" is used in Mac or Linux path. For example, if a path is "C:\\" on Windows machine, it would be "C:/" on Mac or Linux.


Required Programing Language
----------------------------

The MATLAB version of this library can be run by using MATLAB (https://www.mathworks.com) or GNU Octave (https://octave.org). 


Required Package for MATLAB
---------------------------

MATLAB users may need to install additional MATLAB Toolboxes such as Signal Processing Toolbox for some functions.


Required Package for GNU Octave
-------------------------------

GNU Octave users may need to install/load additional packages such as GNU Octave Signal package for some functions.
To find the list of the GNU Octave's pre-installed packages, run the following command in the Command Window:

.. code:: octave
    
    >> pkg list

For example, GNU Octave comes with Signal package but it needs to loaded every time GNU Octave starts. The Signal package can be loaded inside GNU Octave by running the following command in the Command Window (This should be done every time GNU Octave is opened):

.. code:: octave
    
    >> pkg load signal

If GNU Octave Signal Package is not already installed, it should be first installed from https://packages.octave.org, and then get loaded by running the following commands in the Command Window:

.. code:: octave

    >> pkg install "https://downloads.sourceforge.net/project/octave/Octave%20Forge%20Packages/Individual%20Package%20Releases/signal-1.4.5.tar.gz"
    >> pkg load signal


Quick Start
-----------

.. code:: matlab

    x(:,1)=linspace(1,10,10);
    y(:,1)=1+rand(10,1);
    y(:,2)=2+rand(10,1);
    plot2d(x,y,'line_confid','blue_red','large')
