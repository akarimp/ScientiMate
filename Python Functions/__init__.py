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

#Package Information
#-------------------
__version__ = "2.0"
__author__ = 'Arash Karimpour'
#__credits__ = 'xyz Laboratory'    

#----------------------------------------------------------
#Import subdirectory modules
#It uses for module base import such as:
#import scientimate as sm
#from scientimate import plotting
#sm.plotting.plot2d(x,y,'line_confid','blue_red','large')
#----------------------------------------------------------
from . import colormap
from . import data_downloading
from . import data_manipulating
from . import data_reading
from . import hurricane
from . import mapping
from . import ocean_wave_data_analysis
from . import ocean_wave_directional_analysis
from . import ocean_wave_parametric_model
from . import ocean_wave_properties
from . import ocean_wave_spectrum
from . import plotting
from . import signal_processing
from . import statistics
from . import swan_wave_model
from . import wind_engineering

#----------------------------------------------------------
#Import all modules
#It uses for direct use of functions such as:
#import scientimate as sm
#sm.plot2d(x,y,'line_confid','blue_red','large')
#----------------------------------------------------------

#Colormap
#--------
from .colormap.gencolormap import gencolormap
from .colormap.seqcolormap import seqcolormap
from .colormap.topocolormap import topocolormap

#Data Downloading
#----------------
from .data_downloading.downloadndbcdata import downloadndbcdata
from .data_downloading.downloadtidecurrentdata import downloadtidecurrentdata

#Data Manipulating
#-----------------
from .data_manipulating.downsamplex import downsamplex
from .data_manipulating.downsamplexy import downsamplexy
from .data_manipulating.downsamplexyz import downsamplexyz
from .data_manipulating.interpgrid2xyz import interpgrid2xyz
from .data_manipulating.interpxyz2grid import interpxyz2grid
from .data_manipulating.interpxyz2xyz import interpxyz2xyz
from .data_manipulating.replacemissing1d import replacemissing1d
from .data_manipulating.replacemissing2d import replacemissing2d
from .data_manipulating.replaceoutlier import replaceoutlier
from .data_manipulating.replacespike3dps import replacespike3dps
from .data_manipulating.replacespikediff import replacespikediff
from .data_manipulating.replacespikeenvelope import replacespikeenvelope

#Data Reading
#------------
from .data_reading.readasciitable import readasciitable
from .data_reading.readdatafile import readdatafile
from .data_reading.readtimeseriestable import readtimeseriestable
from .data_reading.readxyzfile import readxyzfile

#Hurricane
#---------
from .hurricane.hurricanebackgroundwind import hurricanebackgroundwind
from .hurricane.hurricanedpcpt import hurricanedpcpt
from .hurricane.hurricanepressureh80 import hurricanepressureh80
from .hurricane.hurricanetranslationvel import hurricanetranslationvel
from .hurricane.hurricanewavecontourcem import hurricanewavecontourcem
from .hurricane.hurricanewavecontourh16 import hurricanewavecontourh16
from .hurricane.hurricanewavecontoury88 import hurricanewavecontoury88
from .hurricane.hurricanewavemax import hurricanewavemax
from .hurricane.hurricanewindh08 import hurricanewindh08
from .hurricane.hurricanewindh80 import hurricanewindh80
from .hurricane.hurricanewindinflowangle import hurricanewindinflowangle
from .hurricane.hurricanewindvel import hurricanewindvel
from .hurricane.hurricanewindvelmax import hurricanewindvelmax
from .hurricane.hurricanewindvelmaxh08 import hurricanewindvelmaxh08
from .hurricane.hurricanewindvelmaxh80 import hurricanewindvelmaxh80
from .hurricane.readnhchurricane import readnhchurricane
from .hurricane.stormsurge1d import stormsurge1d

#Mapping
#-------
from .mapping.convertdir import convertdir
from .mapping.distancecart import distancecart
from .mapping.distancegc import distancegc
from .mapping.endpointcart import endpointcart
from .mapping.globalrelief import globalrelief
from .mapping.gridgenerator import gridgenerator
from .mapping.intersectgc import intersectgc
from .mapping.intersectlineedge import intersectlineedge
from .mapping.pointscart import pointscart
from .mapping.reckongc import reckongc
from .mapping.waypointsgc import waypointsgc
from .mapping.windfetch import windfetch
from .mapping.zprofilepath import zprofilepath

#Ocean Wave Data Analysis
#------------------------
from .ocean_wave_data_analysis.diagnostictail import diagnostictail
from .ocean_wave_data_analysis.pressure2surfaceelevfft import pressure2surfaceelevfft
from .ocean_wave_data_analysis.pressure2surfaceelevzcross import pressure2surfaceelevzcross
from .ocean_wave_data_analysis.seaswell1d import seaswell1d
from .ocean_wave_data_analysis.velocity2surfaceelevfft import velocity2surfaceelevfft
from .ocean_wave_data_analysis.velocity2surfaceelevzcross import velocity2surfaceelevzcross
from .ocean_wave_data_analysis.wavefrompressurepsd import wavefrompressurepsd
from .ocean_wave_data_analysis.wavefrompressurezcross import wavefrompressurezcross
from .ocean_wave_data_analysis.wavefromsurfaceelevpsd import wavefromsurfaceelevpsd
from .ocean_wave_data_analysis.wavefromsurfaceelevzcross import wavefromsurfaceelevzcross
from .ocean_wave_data_analysis.wavefromvelocitypsd import wavefromvelocitypsd
from .ocean_wave_data_analysis.wavefromvelocityzcross import wavefromvelocityzcross
from .ocean_wave_data_analysis.wavepropfrompsd import wavepropfrompsd

#Ocean Wave Directional Analysis
#-------------------------------
from .ocean_wave_directional_analysis.directionalpsd import directionalpsd
from .ocean_wave_directional_analysis.directionalpsdetauv import directionalpsdetauv
from .ocean_wave_directional_analysis.directionalpsdpuv import directionalpsdpuv
from .ocean_wave_directional_analysis.enu2truenorth import enu2truenorth
from .ocean_wave_directional_analysis.wavediretauv import wavediretauv
from .ocean_wave_directional_analysis.wavedirpuv import wavedirpuv

#Ocean Wave Parametric Model
#---------------------------
from .ocean_wave_parametric_model.asymptlimit import asymptlimit
from .ocean_wave_parametric_model.equivalentfetchdeep import equivalentfetchdeep
from .ocean_wave_parametric_model.equivalentfetchshallow import equivalentfetchshallow
from .ocean_wave_parametric_model.fullydevwave import fullydevwave
from .ocean_wave_parametric_model.mindurationdeep import mindurationdeep
from .ocean_wave_parametric_model.mindurationshallow import mindurationshallow
from .ocean_wave_parametric_model.parametricwavedeep import parametricwavedeep
from .ocean_wave_parametric_model.parametricwaveshallow import parametricwaveshallow
from .ocean_wave_parametric_model.wavedim2dimless import wavedim2dimless
from .ocean_wave_parametric_model.wavedimless2dim import wavedimless2dim

#Ocean Wave Properties
#---------------------
from .ocean_wave_properties.incidentreflectedwave import incidentreflectedwave
from .ocean_wave_properties.linearwavegenerator import linearwavegenerator
from .ocean_wave_properties.linearwavesuperposition import linearwavesuperposition
from .ocean_wave_properties.pressureresponse import pressureresponse
from .ocean_wave_properties.stokeswavegenerator import stokeswavegenerator
from .ocean_wave_properties.stokeswavesuperposition import stokeswavesuperposition
from .ocean_wave_properties.wavebedstress import wavebedstress
from .ocean_wave_properties.wavedispersion import wavedispersion
from .ocean_wave_properties.wavedispersionds import wavedispersionds
from .ocean_wave_properties.waveorbitalvelocity import waveorbitalvelocity
from .ocean_wave_properties.wavepower import wavepower
from .ocean_wave_properties.wavepowerfrompsd import wavepowerfrompsd
from .ocean_wave_properties.wavespectrum2timeseries import wavespectrum2timeseries
from .ocean_wave_properties.wavevel2wlconvfactor import wavevel2wlconvfactor

#Ocean Wave Spectrum
#-------------------
from .ocean_wave_spectrum.bretpsd import bretpsd
from .ocean_wave_spectrum.donelanpsd import donelanpsd
from .ocean_wave_spectrum.jonswappsd import jonswappsd
from .ocean_wave_spectrum.pmpsd import pmpsd
from .ocean_wave_spectrum.tmapsd import tmapsd

#Plotting
#--------
from .plotting.plot2d import plot2d
from .plotting.plot2dsubplot import plot2dsubplot
from .plotting.plot2dtimeseries import plot2dtimeseries
from .plotting.plot3d import plot3d
from .plotting.plot3ddem import plot3ddem
from .plotting.plot3dhillshades import plot3dhillshades
from .plotting.plot3dtopo import plot3dtopo

#Signal Processing
#-----------------
from .signal_processing.bartlettpsd import bartlettpsd
from .signal_processing.fftfrequency import fftfrequency
from .signal_processing.filtertimeseries import filtertimeseries
from .signal_processing.periodogrampsd import periodogrampsd
from .signal_processing.psd2timeseries import psd2timeseries
from .signal_processing.smoothsignal import smoothsignal
from .signal_processing.spectrogrampsd import spectrogrampsd
from .signal_processing.welchpsd import welchpsd

#Statistics
#----------
from .statistics.curvefit2d import curvefit2d
from .statistics.curvefit3d import curvefit3d
from .statistics.dataoverview import dataoverview
from .statistics.findextremum import findextremum
from .statistics.findknn import findknn
from .statistics.fitgoodness import fitgoodness
from .statistics.levelcrossing import levelcrossing
from .statistics.movingwindow import movingwindow
from .statistics.probability1d import probability1d
from .statistics.probability2d import probability2d
from .statistics.similaritymeasure import similaritymeasure

#Swan Wave Model
#---------------
from .swan_wave_model.swandepthgrid import swandepthgrid
from .swan_wave_model.swanvectorvarspconst import swanvectorvarspconst
from .swan_wave_model.swanvectorvarspvariedgrid import swanvectorvarspvariedgrid
from .swan_wave_model.swanvectorvarspvariedsct import swanvectorvarspvariedsct
from .swan_wave_model.swanwaterlevelspconst import swanwaterlevelspconst
from .swan_wave_model.swanwaterlevelspvariedgrid import swanwaterlevelspvariedgrid
from .swan_wave_model.swanwaterlevelspvariedsct import swanwaterlevelspvariedsct
from .swan_wave_model.swanwindspconst import swanwindspconst
from .swan_wave_model.swanwindspvariedgrid import swanwindspvariedgrid
from .swan_wave_model.swanwindspvariedsct import swanwindspvariedsct

#Wind Engineering
#----------------
from .wind_engineering.directionavg import directionavg
from .wind_engineering.smoothwind import smoothwind
from .wind_engineering.surfaceroughness import surfaceroughness
from .wind_engineering.sustainedwindduration import sustainedwindduration
from .wind_engineering.windavg import windavg
from .wind_engineering.winddrag import winddrag
from .wind_engineering.windgustfactor import windgustfactor
from .wind_engineering.windspectrum import windspectrum
from .wind_engineering.windspectrum2timeseries import windspectrum2timeseries
from .wind_engineering.windvelz1toz2 import windvelz1toz2
