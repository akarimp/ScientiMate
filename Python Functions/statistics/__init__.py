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

#Statistics
#----------
from .curvefit2d import curvefit2d
from .curvefit3d import curvefit3d
from .dataoverview import dataoverview
from .findextremum import findextremum
from .findknn import findknn
from .fitgoodness import fitgoodness
from .levelcrossing import levelcrossing
from .movingwindow import movingwindow
from .probability1d import probability1d
from .probability2d import probability2d
from .similaritymeasure import similaritymeasure
