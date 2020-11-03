API
===

Here is a list of the ScientiMate functions:

Data Manipulating
-----------------
.. toctree::
    :maxdepth: 1

    python_functions/data_manipulating/downsamplex.rst
    python_functions/data_manipulating/interpxyz2grid.rst
    python_functions/data_manipulating/replacemissing1d.rst
    python_functions/data_manipulating/replaceoutlier.rst

.. #downsamplexy.rst
.. #downsamplexyz.rst
.. #interpgrid2xyz.rst
.. #interpxyz2xyz.rst
.. #readdatafile.rst
.. #replacemissing2d.rst
.. #replacespike3dps.rst
.. #replacespikediff.rst
.. #replacespikeenvelope.rst


Hurricane
---------
.. toctree::
    :maxdepth: 1
    
    python_functions/hurricane/stormsurge1d.rst
    python_functions/hurricane/hurricanewindh80.rst

.. from .hurricanebackgroundwind import hurricanebackgroundwind
.. from .hurricanedpcpt import hurricanedpcpt
.. from .hurricanepressureh80 import hurricanepressureh80
.. from .hurricanetranslationvel import hurricanetranslationvel
.. from .hurricanewavecontourcem import hurricanewavecontourcem
.. from .hurricanewavecontourh16 import hurricanewavecontourh16
.. from .hurricanewavecontoury88 import hurricanewavecontoury88
.. from .hurricanewavemax import hurricanewavemax
.. from .hurricanewindh08 import hurricanewindh08
.. from .hurricanewindinflowangle import hurricanewindinflowangle
.. from .hurricanewindvel import hurricanewindvel
.. from .hurricanewindvelmax import hurricanewindvelmax
.. from .hurricanewindvelmaxh08 import hurricanewindvelmaxh08
.. from .hurricanewindvelmaxh80 import hurricanewindvelmaxh80
.. from .readnhchurricane import readnhchurricane

Mapping
-------
.. toctree::
    :maxdepth: 1
    
    python_functions/mapping/convertdir.rst
    python_functions/mapping/distancecart.rst
    python_functions/mapping/distancegc.rst
    python_functions/mapping/endpointcart.rst
    python_functions/mapping/globalrelief.rst
    python_functions/mapping/gridgenerator.rst
    python_functions/mapping/readxyz.rst
    python_functions/mapping/reckongc.rst

.. from .intersectgc import intersectgc
.. from .intersectlineedge import intersectlineedge
.. from .pointscart import pointscart
.. from .waypointsgc import waypointsgc
.. from .windfetch import windfetch
.. from .zprofilepath import zprofilepath

OCEANLYZ
--------
| OCEANLYZ, Ocean Wave Analyzing Toolbox, is a toolbox for analyzing the wave time series data collected by sensors in open body of water such as ocean, sea, and lake or in a laboratory.
| The Python version of OCEANLYZ toolbox is a part of the ScintiMate package.
| For more information, visit https://oceanlyz.readthedocs.io

.. toctree::
    :maxdepth: 1

    python_functions/oceanlyz/oceanlyz.rst
    python_functions/oceanlyz/PcorFFTFun.rst
    python_functions/oceanlyz/PcorZerocrossingFun.rst
    python_functions/oceanlyz/SeaSwellFun.rst
    python_functions/oceanlyz/WaveSpectraFun.rst
    python_functions/oceanlyz/WaveZerocrossingFun.rst

OCEANLYZ toolbox consists of 1 class and 5 functions.
The main class in OCEANLYZ toolbox is the oceanlyz(). To run OCEANLYZ toolbox, only the oceanlyz() class is required to be run.
Based on parameters set by a user, oceanlyz() calls appropriate function(s) to analyze data.
Note that, any of the OCEANLYZ functions might be used separately as well.

For more general purpose functions for wave analysis, see `Water Wave Data Analysis`_ section below.

Plotting
--------
.. toctree::
    :maxdepth: 1

    python_functions/plotting/plot2d.rst
    python_functions/plotting/plot2dsubplot.rst
    python_functions/plotting/plot2dtimeseries.rst
    python_functions/plotting/plot3d.rst

.. from .plot3ddem import plot3ddem
.. from .plot3dhillshades import plot3dhillshades
.. from .plot3dtopo import plot3dtopo

Signal Processing
-----------------
.. toctree::
    :maxdepth: 1

    python_functions/signal_processing/filtertimeseries.rst
    python_functions/signal_processing/smoothsignal.rst

.. from .bartlettpsd import bartlettpsd
.. from .fftfrequency import fftfrequency
.. from .periodogrampsd import periodogrampsd
.. from .psd2timeseries import psd2timeseries
.. from .spectrogrampsd import spectrogrampsd
.. from .welchpsd import welchpsd

Statistics
----------
.. toctree::
    :maxdepth: 1

    python_functions/statistics/curvefit2d.rst
    python_functions/statistics/dataoverview.rst
    python_functions/statistics/findextremum.rst
    python_functions/statistics/findknn.rst
    python_functions/statistics/fitgoodness.rst
    python_functions/statistics/movingwindow.rst
    python_functions/statistics/probability1d.rst
    python_functions/statistics/similaritymeasure.rst

.. from .curvefit3d import curvefit3d
.. from .levelcrossing import levelcrossing
.. from .probability2d import probability2d

Swan
----
.. toctree::
    :maxdepth: 1

    python_functions/swan/swandepthgrid.rst
    python_functions/swan/swanvectorvarspconst.rst
    python_functions/swan/swanwaterlevelspconst.rst
    python_functions/swan/swanwindspconst.rst

.. from .swanvectorvarspvariedgrid import swanvectorvarspvariedgrid
.. from .swanvectorvarspvariedsct import swanvectorvarspvariedsct
.. from .swanwaterlevelspvariedgrid import swanwaterlevelspvariedgrid
.. from .swanwaterlevelspvariedsct import swanwaterlevelspvariedsct
.. from .swanwindspvariedgrid import swanwindspvariedgrid
.. from .swanwindspvariedsct import swanwindspvariedsct

Water Wave Data Analysis
------------------------
.. toctree::
    :maxdepth: 1

    python_functions/water_wave_data_analysis/diagnostictail.rst
    python_functions/water_wave_data_analysis/seaswell1d.rst
    python_functions/water_wave_data_analysis/wavefrompressurepsd.rst
    python_functions/water_wave_data_analysis/wavefrompressurezcross.rst
    python_functions/water_wave_data_analysis/wavefromsurfaceelevpsd.rst
    python_functions/water_wave_data_analysis/wavefromsurfaceelevzcross.rst
    python_functions/water_wave_data_analysis/wavepropfrompsd.rst

.. from .pressure2surfaceelevfft import pressure2surfaceelevfft
.. from .pressure2surfaceelevzcross import pressure2surfaceelevzcross
.. from .velocity2surfaceelevfft import velocity2surfaceelevfft
.. from .velocity2surfaceelevzcross import velocity2surfaceelevzcross
.. from .wavefromvelocitypsd import wavefromvelocitypsd
.. from .wavefromvelocityzcross import wavefromvelocityzcross

For Ocean Wave Analyzing Toolbox, see `OCEANLYZ`_ section above.

Water Wave Directional Analysis
-------------------------------
.. toctree::
    :maxdepth: 1

    python_functions/water_wave_directional_analysis/directionalpsd.rst
    python_functions/water_wave_directional_analysis/enu2truenorth.rst

.. from .directionalpsdetauv import directionalpsdetauv
.. from .directionalpsdpuv import directionalpsdpuv
.. from .wavediretauv import wavediretauv
.. from .wavedirpuv import wavedirpuv

Water Wave Parametric Model
---------------------------
.. toctree::
    :maxdepth: 1

    python_functions/water_wave_parametric_model/parametricwavedeep.rst
    python_functions/water_wave_parametric_model/parametricwaveshallow.rst

.. from .asymptlimit import asymptlimit
.. from .equivfetchdeep import equivfetchdeep
.. from .equivfetchshallow import equivfetchshallow
.. from .fullydevwave import fullydevwave
.. from .mindurationdeep import mindurationdeep
.. from .mindurationshallow import mindurationshallow
.. from .wavedim2dimless import wavedim2dimless
.. from .wavedimless2dim import wavedimless2dim

Water Wave Properties
---------------------
.. toctree::
    :maxdepth: 1

    python_functions/water_wave_properties/bretpsd.rst
    python_functions/water_wave_properties/donelanpsd.rst
    python_functions/water_wave_properties/jonswappsd.rst
    python_functions/water_wave_properties/pmpsd.rst
    python_functions/water_wave_properties/pressureresponse.rst
    python_functions/water_wave_properties/tmapsd.rst
    python_functions/water_wave_properties/wavedispersion.rst

.. from .incidentreflectedwave import incidentreflectedwave
.. from .linearwavegenerator import linearwavegenerator
.. from .linearwavesuperposition import linearwavesuperposition
.. from .stokeswavegenerator import stokeswavegenerator
.. from .stokeswavesuperposition import stokeswavesuperposition
.. from .wavebedstress import wavebedstress
.. from .wavedispersionds import wavedispersionds
.. from .waveorbitalvelocity import waveorbitalvelocity
.. from .wavepower import wavepower
.. from .wavepowerfrompsd import wavepowerfrompsd
.. from .wavespectrum2timeseries import wavespectrum2timeseries
.. from .wavevel2wlconvfactor import wavevel2wlconvfactor

Wind
----
.. toctree::
    :maxdepth: 1

    python_functions/wind/directionavg.rst
    python_functions/wind/surfaceroughness.rst
    python_functions/wind/sustainedwindduration.rst
    python_functions/wind/windavg.rst
    python_functions/wind/winddrag.rst
    python_functions/wind/windgustfactor.rst
    python_functions/wind/windspectrum
    python_functions/wind/windvelz1toz2.rst

.. from .windspectrum2timeseries import windspectrum2timeseries
