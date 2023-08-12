"""
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2023-08-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
"""

#----------------------------------------------------------
#Import subdirectory modules
#It uses for module base import such as:
#import scientimate as sm
#from scientimate import plotting
#sm.plotting.plot2d(x,y,'line_confid','blue_red','large')
#----------------------------------------------------------

#Signal Processing
#-----------------
from .bartlettpsd import bartlettpsd
from .fftfrequency import fftfrequency
from .filtertimeseries import filtertimeseries
from .periodogrampsd import periodogrampsd
from .psd2timeseries import psd2timeseries
from .smoothsignal import smoothsignal
from .spectrogrampsd import spectrogrampsd
from .welchpsd import welchpsd
