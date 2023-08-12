%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2022-04-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

add_scientimate_to_path
=======================

.. code:: MATLAB

    add_scientimate_to_path

Description
-----------

Add ScientiMate folder to MATLAB or GNU Octave path

Inputs
------

Outputs
-------

Examples
--------

.. code:: MATLAB

    add_scientimate_to_path

References
----------

.. License & Disclaimer
.. --------------------
..
.. Copyright (c) 2022 Arash Karimpour
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
%}

%--------------------------------------------------------------------------
%CODE
%--------------------------------------------------------------------------
%Add ScientiMate folder to MATLAB or GNU Octave path

%[filepath, name, ext] = fileparts(pwd); %Get ScientiMate parent folder path
ScientiMateFolder = pwd; %Get ScientiMate folder path
ScientiMatePath = genpath(ScientiMateFolder); %Generating path for ScientiMate folder and its sub-folders
addpath(ScientiMatePath); %Add ScientiMate folder to path

disp('ScientiMate folder is added to path')
