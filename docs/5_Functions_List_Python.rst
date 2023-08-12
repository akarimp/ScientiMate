Functions List (Python)
=======================

Here is a list of the ScientiMate functions (Python):

Colormap
--------

.. table::
   :widths: auto
   :align: left

   ====================================================================================== ===========
   Function                                                                               Description
   ====================================================================================== ===========
   :doc:`python_functions/colormap/gencolormap`                                           Generate a colormap from input colors
   :doc:`python_functions/colormap/seqcolormap`                                           Generate sequential colormap for drawing lines
   :doc:`python_functions/colormap/topocolormap`                                          Export a topographic colormap
   ====================================================================================== ===========

.. toctree::
    :maxdepth: 1
    :hidden:

    python_functions/colormap/gencolormap.rst
    python_functions/colormap/seqcolormap.rst
    python_functions/colormap/topocolormap.rst

Data Downloading
-----------------

.. table::
   :widths: auto
   :align: left

   ====================================================================================== ===========
   Function                                                                               Description
   ====================================================================================== ===========
   :doc:`python_functions/data_downloading/downloadndbcdata`                              Download meteorological data from NOAA National Data Buoy Center
   :doc:`python_functions/data_downloading/downloadtidecurrentdata`                       Download meteorological data from NOAA’s Center for Operational Oceanographic Products and Services (CO-OPS)
   ====================================================================================== ===========

.. toctree::
    :maxdepth: 1
    :hidden:

    python_functions/data_downloading/downloadndbcdata.rst
    python_functions/data_downloading/downloadtidecurrentdata.rst

Data Manipulating
-----------------

.. table::
   :widths: auto
   :align: left

   ====================================================================================== ===========
   Function                                                                               Description
   ====================================================================================== ===========
   :doc:`python_functions/data_manipulating/downsamplex`                                  Downsample x data and retain given ratio
   :doc:`python_functions/data_manipulating/downsamplexy`                                 Downsample x and y data and retain given ratio
   :doc:`python_functions/data_manipulating/downsamplexyz`                                Downsample x, y,  and z data and retain given ratio
   :doc:`python_functions/data_manipulating/interpgrid2xyz`                               Interpolate 2d gridded data on given scatter point(s) using nearest neighbor method
   :doc:`python_functions/data_manipulating/interpxyz2grid`                               Interpolate x (longitude), y (latitude) and z (elevation) data into a defined mesh
   :doc:`python_functions/data_manipulating/interpxyz2xyz`                                Interpolate 2d scattered data on given point(s) by down-sampling the input data
   :doc:`python_functions/data_manipulating/replacemissing1d`                             Replace missing data points in 1d data such as time series
   :doc:`python_functions/data_manipulating/replacemissing2d`                             Replace missing data points in 2d array
   :doc:`python_functions/data_manipulating/replaceoutlier`                               Remove outliers in the time series using moving z-score window
   :doc:`python_functions/data_manipulating/replacespike3dps`                             Remove spikes in the time series based on 3D phase space method by Goring and Nikora (2002)
   :doc:`python_functions/data_manipulating/replacespikediff`                             Remove spikes in the time series using a local difference of data respect to a moving average window
   :doc:`python_functions/data_manipulating/replacespikeenvelope`                         Remove spikes in the time series that are outside a defined envelope
   ====================================================================================== ===========

.. toctree::
    :maxdepth: 1
    :hidden:

    python_functions/data_manipulating/downsamplex.rst
    python_functions/data_manipulating/downsamplexy.rst
    python_functions/data_manipulating/downsamplexyz.rst
    python_functions/data_manipulating/interpgrid2xyz.rst
    python_functions/data_manipulating/interpxyz2grid.rst
    python_functions/data_manipulating/interpxyz2xyz.rst
    python_functions/data_manipulating/replacemissing1d.rst
    python_functions/data_manipulating/replacemissing2d.rst
    python_functions/data_manipulating/replaceoutlier.rst
    python_functions/data_manipulating/replacespike3dps.rst
    python_functions/data_manipulating/replacespikediff.rst
    python_functions/data_manipulating/replacespikeenvelope.rst

Data Reading
------------

.. table::
   :widths: auto
   :align: left

   ====================================================================================== ===========
   Function                                                                               Description
   ====================================================================================== ===========
   :doc:`python_functions/data_reading/readasciitable`                                    Read and extract data from ASCII, text, Comma Separated Values (CSV), and spreadsheet (xlsx, xls, ods) file
   :doc:`python_functions/data_reading/readdatafile`                                      Read and extract data from ASCII, text, Comma Separated Values (CSV), Matlab .mat file
   :doc:`python_functions/data_reading/readtimeseriesfile`                                Read and extract time-series data from ASCII, text, Comma Separated Values (CSV), and spreadsheet (xlsx, xls, ods) file
   :doc:`python_functions/data_reading/readxyzfile`                                       Read and extract x (longitude), y (latitude) and z (elevation) data from ASCII gridded (tabular) xyz file
   ====================================================================================== ===========

.. toctree::
    :maxdepth: 1
    :hidden:

    python_functions/data_reading/readasciitable.rst
    python_functions/data_reading/readdatafile.rst
    python_functions/data_reading/readtimeseriesfile.rst
    python_functions/data_reading/readxyzfile.rst

Hurricane
---------

.. table::
   :widths: auto
   :align: left

   ====================================================================================== ===========
   Function                                                                               Description
   ====================================================================================== ===========
   :doc:`python_functions/hurricane/hurricanebackgroundwind`                              Calculate and add background wind velocity due to hurricane front motion to hurricane rotational wind velocity
   :doc:`python_functions/hurricane/hurricanedpcpt`                                       Calculate hurricane central pressure (Pc) intensity change over time (dPc/dt)
   :doc:`python_functions/hurricane/hurricanepressureh80`                                 Generate hurricane pressure data on given (x,y) points using method from Holland (1980)
   :doc:`python_functions/hurricane/hurricanetranslationvel`                              Calculate hurricane center translational (forward motion) velocity
   :doc:`python_functions/hurricane/hurricanewavecontourcem`                              Calculates hurricane wave height field (contours) on given mesh using method from Shore Protection Manual (SPM), U.S. Army Corps of Engineers (1984), and Young (1988)
   :doc:`python_functions/hurricane/hurricanewavecontourh16`                              Calculates hurricane wave height and wave period field (contours) on given mesh using method from Hwang (2016) and Hwang & Walsh (2016)
   :doc:`python_functions/hurricane/hurricanewavecontoury88`                              Calculate hurricane wave height field (contours) on given mesh using method from Young (1988)
   :doc:`python_functions/hurricane/hurricanewavemax`                                     Calculate hurricane maximum wave height and wave period at a location of maximum wind
   :doc:`python_functions/hurricane/hurricanewindh08`                                     Generate hurricane wind and pressure data on (x,y) points using method from Holland (2008)
   :doc:`python_functions/hurricane/hurricanewindh80`                                     Generate hurricane wind and pressure data on given (x,y) points using method from Holland (1980)
   :doc:`python_functions/hurricane/hurricanewindinflowangle`                             Calculate hurricane velocity tangential and inflow angle and inflow velocity in x (East) and y (North) directions
   :doc:`python_functions/hurricane/hurricanewindvel`                                     Generate hurricane wind velocity data on given (x,y) points
   :doc:`python_functions/hurricane/hurricanewindvelmax`                                  Calculate hurricane maximum wind velocity at the surface level
   :doc:`python_functions/hurricane/hurricanewindvelmaxh08`                               Calculate hurricane maximum wind velocity at the gradient level using Holland (2008) method
   :doc:`python_functions/hurricane/hurricanewindvelmaxh80`                               Calculate hurricane maximum wind velocity at the gradient level using Holland (1980) method
   :doc:`python_functions/hurricane/readnhchurricane`                                     Read and extracts hurricane data from National Hurricane Center (NHC) HURDAT2 file
   :doc:`python_functions/hurricane/stormsurge1d`                                         Calculate one dimensional storm surge using Dean Dalrymple (1991) method
   ====================================================================================== ===========

.. toctree::
    :maxdepth: 1
    :hidden:
    
    python_functions/hurricane/hurricanebackgroundwind.rst
    python_functions/hurricane/hurricanedpcpt.rst
    python_functions/hurricane/hurricanepressureh80.rst
    python_functions/hurricane/hurricanetranslationvel.rst
    python_functions/hurricane/hurricanewavecontourcem.rst
    python_functions/hurricane/hurricanewavecontourh16.rst
    python_functions/hurricane/hurricanewavecontoury88.rst
    python_functions/hurricane/hurricanewavemax.rst
    python_functions/hurricane/hurricanewindh08.rst
    python_functions/hurricane/hurricanewindh80.rst
    python_functions/hurricane/hurricanewindinflowangle.rst
    python_functions/hurricane/hurricanewindvel.rst
    python_functions/hurricane/hurricanewindvelmax.rst
    python_functions/hurricane/hurricanewindvelmaxh08.rst
    python_functions/hurricane/hurricanewindvelmaxh80.rst
    python_functions/hurricane/readnhchurricane.rst
    python_functions/hurricane/stormsurge1d.rst

Mapping
-------

.. table::
   :widths: auto
   :align: left

   ====================================================================================== ===========
   Function                                                                               Description
   ====================================================================================== ===========
   :doc:`python_functions/mapping/convertdir`                                             Convert direction from one system to another one
   :doc:`python_functions/mapping/distancecart`                                           Calculate distance from (x1,y1) to (x2,y2) on cartesian coordinate
   :doc:`python_functions/mapping/distancegc`                                             Calculate distance and azimuth (bearing) between (Latitude,Longitude) points using Great Circle
   :doc:`python_functions/mapping/endpointcart`                                           Find an end point of the straight line segment from its starting point (x1,y1) and its angle on cartesian coordinate
   :doc:`python_functions/mapping/globalrelief`                                           Return x (longitude), y (latitude) and z (elevation) data from ETOPO1 Global Relief Model (Amante & Eakins, 2009) interpolated on 0.08 degree grid
   :doc:`python_functions/mapping/gridgenerator`                                          Generate 2d x-y grid
   :doc:`python_functions/mapping/intersectgc`                                            Find intersection point between two line segments (line edges) on Great Circle
   :doc:`python_functions/mapping/intersectlineedge`                                      Find intersection point between two line segments (line edges)
   :doc:`python_functions/mapping/pointscart`                                             Generate points between point (x1,y1) and (x2,y2) on cartesian coordinate
   :doc:`python_functions/mapping/reckongc`                                               Calculate end point (Latitude,Longitude) from start point (Latitude,Longitude) and distance and azimuth (bearing) using Great Circle
   :doc:`python_functions/mapping/waypointsgc`                                            Generate (Latitude,Longitude) points between two (Latitude,Longitude) points using Great Circle
   :doc:`python_functions/mapping/windfetch`                                              Calculate a wind fecth and z (elevation) profile along a path over water for a given 2d x-y domain (map, image, ...)
   :doc:`python_functions/mapping/zprofilepath`                                           Calculate z (elevation, ...) profile along a path over a given 2d x-y domain (map, image, ...)
   ====================================================================================== ===========

.. toctree::
    :maxdepth: 1
    :hidden:
    
    python_functions/mapping/convertdir.rst
    python_functions/mapping/distancecart.rst
    python_functions/mapping/distancegc.rst
    python_functions/mapping/endpointcart.rst
    python_functions/mapping/globalrelief.rst
    python_functions/mapping/gridgenerator.rst
    python_functions/mapping/intersectgc.rst
    python_functions/mapping/intersectlineedge.rst
    python_functions/mapping/pointscart.rst
    python_functions/mapping/reckongc.rst
    python_functions/mapping/waypointsgc.rst
    python_functions/mapping/windfetch.rst
    python_functions/mapping/zprofilepath.rst

Ocean Wave Data Analysis
------------------------

.. table::
   :widths: auto
   :align: left

   ====================================================================================== ===========
   Function                                                                               Description
   ====================================================================================== ===========
   :doc:`python_functions/ocean_wave_data_analysis/diagnostictail`                        Replace a spectrum tail with JONSWAP (Hasselmann et al.,  1973) or TMA Spectrum (Bouws et al., 1985)
   :doc:`python_functions/ocean_wave_data_analysis/pressure2surfaceelevfft`               Calculate water surface elevation time series from water pressure time series by using Fast Fourier Transform
   :doc:`python_functions/ocean_wave_data_analysis/pressure2surfaceelevzcross`            Calculate water surface elevation time series from water pressure time series by using an upward zero crossing method
   :doc:`python_functions/ocean_wave_data_analysis/seaswell1d`                            Partition (separate) wind sea from swell in a power spectral density using an one dimensional method
   :doc:`python_functions/ocean_wave_data_analysis/velocity2surfaceelevfft`               Calculate water surface elevation time series from wave orbital velocity time series by using Fast Fourier Transform
   :doc:`python_functions/ocean_wave_data_analysis/velocity2surfaceelevzcross`            Calculate wave properties from wave orbital velocity by using an upward zero crossing method
   :doc:`python_functions/ocean_wave_data_analysis/wavefrompressurepsd`                   Calculate wave properties from water pressure by converting it to water surface elevation power spectral density
   :doc:`python_functions/ocean_wave_data_analysis/wavefrompressurezcross`                Calculate wave properties from water pressure by using an upward zero crossing method
   :doc:`python_functions/ocean_wave_data_analysis/wavefromsurfaceelevpsd`                Calculate wave properties from water surface elevation power spectral density
   :doc:`python_functions/ocean_wave_data_analysis/wavefromsurfaceelevzcross`             Calculate wave properties from water surface elevation by using an upward zero crossing method
   :doc:`python_functions/ocean_wave_data_analysis/wavefromvelocitypsd`                   Calculate wave properties from wave orbital velocity by converting it to water surface elevation power spectral density
   :doc:`python_functions/ocean_wave_data_analysis/wavefromvelocityzcross`                Calculate wave properties from wave orbital velocity by using an upward zero crossing method
   :doc:`python_functions/ocean_wave_data_analysis/wavepropfrompsd`                       Calculate wave properties from a power spectral density
   ====================================================================================== ===========

.. toctree::
    :maxdepth: 1
    :hidden:

    python_functions/ocean_wave_data_analysis/diagnostictail.rst
    python_functions/ocean_wave_data_analysis/pressure2surfaceelevfft.rst
    python_functions/ocean_wave_data_analysis/pressure2surfaceelevzcross.rst
    python_functions/ocean_wave_data_analysis/seaswell1d.rst
    python_functions/ocean_wave_data_analysis/velocity2surfaceelevfft.rst
    python_functions/ocean_wave_data_analysis/velocity2surfaceelevzcross.rst
    python_functions/ocean_wave_data_analysis/wavefrompressurepsd.rst
    python_functions/ocean_wave_data_analysis/wavefrompressurezcross.rst
    python_functions/ocean_wave_data_analysis/wavefromsurfaceelevpsd.rst
    python_functions/ocean_wave_data_analysis/wavefromsurfaceelevzcross.rst
    python_functions/ocean_wave_data_analysis/wavefromvelocitypsd.rst
    python_functions/ocean_wave_data_analysis/wavefromvelocityzcross.rst
    python_functions/ocean_wave_data_analysis/wavepropfrompsd.rst

| For ocean wave data analysis, you may use **OCEANLYZ** toolbox as well.
| OCEANLYZ, Ocean Wave Analyzing Toolbox, is a toolbox for analyzing the wave time series data collected by sensors in open body of water such as ocean, sea, and lake or in a laboratory.
| For more information, visit https://github.com/akarimp/Oceanlyz

Ocean Wave Directional Analysis
-------------------------------

.. table::
   :widths: auto
   :align: left

   ====================================================================================== ===========
   Function                                                                               Description
   ====================================================================================== ===========
   :doc:`python_functions/ocean_wave_directional_analysis/directionalpsd`                 Calculate wave directional spectrum using parametric directional spreading function
   :doc:`python_functions/ocean_wave_directional_analysis/directionalpsdetauv`            Calculate wave directional spectrum using water surface elevation and horizontal orbital velocity
   :doc:`python_functions/ocean_wave_directional_analysis/directionalpsdpuv`              Calculate wave directional spectrum using pressure and horizontal orbital velocity
   :doc:`python_functions/ocean_wave_directional_analysis/enu2truenorth`                  Convert mathematical direction (angle) from ENU (East North Up) coordinate system to compass direction with respect to true north
   :doc:`python_functions/ocean_wave_directional_analysis/wavediretauv`                   Calculate wave direction using water surface elevation and horizontal orbital velocity
   :doc:`python_functions/ocean_wave_directional_analysis/wavedirpuv`                     Calculate wave direction using pressure and horizontal orbital velocity
   ====================================================================================== ===========

.. toctree::
    :maxdepth: 1
    :hidden:

    python_functions/ocean_wave_directional_analysis/directionalpsd.rst
    python_functions/ocean_wave_directional_analysis/directionalpsdetauv.rst
    python_functions/ocean_wave_directional_analysis/directionalpsdpuv.rst
    python_functions/ocean_wave_directional_analysis/enu2truenorth.rst
    python_functions/ocean_wave_directional_analysis/wavediretauv.rst
    python_functions/ocean_wave_directional_analysis/wavedirpuv.rst

Ocean Wave Parametric Model
---------------------------

.. table::
   :widths: auto
   :align: left

   ====================================================================================== ===========
   Function                                                                               Description
   ====================================================================================== ===========
   :doc:`python_functions/ocean_wave_parametric_model/asymptlimit`                        Calculate dimensionless asymptotic limits of wind wave growth in shallow and intermediate water
   :doc:`python_functions/ocean_wave_parametric_model/equivalentfetchdeep`                Calculate an equivalent wind fetch for duration limited wave growth in deep water
   :doc:`python_functions/ocean_wave_parametric_model/equivalentfetchshallow`             Calculate an equivalent wind fetch for duration limited wave growth in shallow and intermediate water
   :doc:`python_functions/ocean_wave_parametric_model/fullydevwave`                       Calculate a fully developed condition for wind wave growth
   :doc:`python_functions/ocean_wave_parametric_model/mindurationdeep`                    Calculate a minimum required wind duration for wave to be fetch limited in deep water
   :doc:`python_functions/ocean_wave_parametric_model/mindurationshallow`                 Calculate a minimum required wind duration for wave to be fetch limited in shallow and intermediate water
   :doc:`python_functions/ocean_wave_parametric_model/parametricwavedeep`                 Calculate wave properties using parametric wave models in deep water
   :doc:`python_functions/ocean_wave_parametric_model/parametricwaveshallow`              Calculate wave properties using parametric wave models in shallow and intermediate water
   :doc:`python_functions/ocean_wave_parametric_model/wavedim2dimless`                    Calculate dimensionless numbers from dimensional numbers
   :doc:`python_functions/ocean_wave_parametric_model/wavedimless2dim`                    Calculate dimensional numbers from dimensionless numbers
   ====================================================================================== ===========

.. toctree::
    :maxdepth: 1
    :hidden:

    python_functions/ocean_wave_parametric_model/asymptlimit.rst
    python_functions/ocean_wave_parametric_model/equivalentfetchdeep.rst
    python_functions/ocean_wave_parametric_model/equivalentfetchshallow.rst
    python_functions/ocean_wave_parametric_model/fullydevwave.rst
    python_functions/ocean_wave_parametric_model/mindurationdeep.rst
    python_functions/ocean_wave_parametric_model/mindurationshallow.rst
    python_functions/ocean_wave_parametric_model/parametricwavedeep.rst
    python_functions/ocean_wave_parametric_model/parametricwaveshallow.rst
    python_functions/ocean_wave_parametric_model/wavedim2dimless.rst
    python_functions/ocean_wave_parametric_model/wavedimless2dim.rst

Ocean Wave Properties
---------------------

.. table::
   :widths: auto
   :align: left

   ====================================================================================== ===========
   Function                                                                               Description
   ====================================================================================== ===========
   :doc:`python_functions/ocean_wave_properties/incidentreflectedwave`                    Separate incident and reflected waves
   :doc:`python_functions/ocean_wave_properties/linearwavegenerator`                      Generate linear waves
   :doc:`python_functions/ocean_wave_properties/linearwavesuperposition`                  Superposition linear waves
   :doc:`python_functions/ocean_wave_properties/pressureresponse`                         Calculate a pressure response factor
   :doc:`python_functions/ocean_wave_properties/stokeswavegenerator`                      Generate second order stokes' waves
   :doc:`python_functions/ocean_wave_properties/stokeswavesuperposition`                  Superposition second order stokes' waves
   :doc:`python_functions/ocean_wave_properties/wavebedstress`                            Calculate the bottom shear velocity and shear stress from current velocity and wave
   :doc:`python_functions/ocean_wave_properties/wavedispersion`                           Solve water wave dispersion relation
   :doc:`python_functions/ocean_wave_properties/wavedispersionds`                         Solve water wave dispersion relation with presence of current (Doppler shift)
   :doc:`python_functions/ocean_wave_properties/waveorbitalvelocity`                      Calculate maximum wave orbital velocity and maximum wave orbital excursion using linear wave theory
   :doc:`python_functions/ocean_wave_properties/wavepower`                                Calculate wave power
   :doc:`python_functions/ocean_wave_properties/wavepowerfrompsd`                         Calculate wave energy and wave power from power spectral density
   :doc:`python_functions/ocean_wave_properties/wavespectrum2timeseries`                  Generate random water wave data from a given water wave spectrum using wave superposition
   :doc:`python_functions/ocean_wave_properties/wavevel2wlconvfactor`                     Calculate a water particle horizontal orbital velocity to the water surface elevation conversion factor
   ====================================================================================== ===========

.. toctree::
    :maxdepth: 1
    :hidden:

    python_functions/ocean_wave_properties/incidentreflectedwave.rst
    python_functions/ocean_wave_properties/linearwavegenerator.rst
    python_functions/ocean_wave_properties/linearwavesuperposition.rst
    python_functions/ocean_wave_properties/pressureresponse.rst
    python_functions/ocean_wave_properties/stokeswavegenerator.rst
    python_functions/ocean_wave_properties/stokeswavesuperposition.rst
    python_functions/ocean_wave_properties/wavebedstress.rst
    python_functions/ocean_wave_properties/wavedispersion.rst
    python_functions/ocean_wave_properties/wavedispersionds.rst
    python_functions/ocean_wave_properties/waveorbitalvelocity.rst
    python_functions/ocean_wave_properties/wavepower.rst
    python_functions/ocean_wave_properties/wavepowerfrompsd.rst
    python_functions/ocean_wave_properties/wavespectrum2timeseries.rst
    python_functions/ocean_wave_properties/wavevel2wlconvfactor.rst

Ocean Wave Spectrum
-------------------

.. table::
   :widths: auto
   :align: left

   ====================================================================================== ===========
   Function                                                                               Description
   ====================================================================================== ===========
   :doc:`python_functions/ocean_wave_spectrum/bretpsd`                                    Calculate Bretschneider spectrum (power spectral density), (Bretschneider, 1959), ITTC spectrum
   :doc:`python_functions/ocean_wave_spectrum/donelanpsd`                                 Calculate Donelan spectrum (power spectral density), (Donelan et al. 1985)
   :doc:`python_functions/ocean_wave_spectrum/jonswappsd`                                 Calculate JONSWAP spectrum (power spectral density), (Hasselmann et al. 1973)
   :doc:`python_functions/ocean_wave_spectrum/pmpsd`                                      Calculate Pierson-Moskowitz spectrum (power spectral density), (Pierson and Moskowitz 1964)
   :doc:`python_functions/ocean_wave_spectrum/tmapsd`                                     Calculate TMA spectrum (power spectral density), (Bouws et al. 1985)
   ====================================================================================== ===========

.. toctree::
    :maxdepth: 1
    :hidden:

    python_functions/ocean_wave_spectrum/bretpsd.rst
    python_functions/ocean_wave_spectrum/donelanpsd.rst
    python_functions/ocean_wave_spectrum/jonswappsd.rst
    python_functions/ocean_wave_spectrum/pmpsd.rst
    python_functions/ocean_wave_spectrum/tmapsd.rst

Plotting
--------

.. table::
   :widths: auto
   :align: left

   ====================================================================================== ===========
   Function                                                                               Description
   ====================================================================================== ===========
   :doc:`python_functions/plotting/plot2d`                                                Plot x and y data in 2-d plot
   :doc:`python_functions/plotting/plot2dsubplot`                                         Plot x and y data in 2-d subplots
   :doc:`python_functions/plotting/plot2dtimeseries`                                      Plot x data in 2-d timeseries
   :doc:`python_functions/plotting/plot3d`                                                Plot x , y, and z data in 2-d/3-d contour/surface plot
   :doc:`python_functions/plotting/plot3ddem`                                             Plot x (longitude), y (latitude) and z (elevation) data into a defined mesh
   :doc:`python_functions/plotting/plot3dhillshades`                                      Plot hillshades (shaded relief) of x (longitude), y (latitude) and z (elevation) data
   :doc:`python_functions/plotting/plot3dtopo`                                            Plot x (longitude), y (latitude) and z (elevation) data into a defined mesh
   ====================================================================================== ===========

.. toctree::
    :maxdepth: 1
    :hidden:

    python_functions/plotting/plot2d.rst
    python_functions/plotting/plot2dsubplot.rst
    python_functions/plotting/plot2dtimeseries.rst
    python_functions/plotting/plot3d.rst
    python_functions/plotting/plot3ddem.rst
    python_functions/plotting/plot3dhillshades.rst
    python_functions/plotting/plot3dtopo.rst

Signal Processing
-----------------

.. table::
   :widths: auto
   :align: left

   ====================================================================================== ===========
   Function                                                                               Description
   ====================================================================================== ===========
   :doc:`python_functions/signal_processing/bartlettpsd`                                  Calculate power spectral density using Bartlett's method
   :doc:`python_functions/signal_processing/fftfrequency`                                 Return frequencies for Fast Fourier Transform
   :doc:`python_functions/signal_processing/filtertimeseries`                             Filter time-series to retain signals with frequencies of fcL <= f <= fcH
   :doc:`python_functions/signal_processing/periodogrampsd`                               Calculate power spectral density using periodogram method
   :doc:`python_functions/signal_processing/psd2timeseries`                               Generate random wave data from a given spectrum
   :doc:`python_functions/signal_processing/smoothsignal`                                 Smooth input data using a window function
   :doc:`python_functions/signal_processing/spectrogrampsd`                               Calculate spectrogram following Welch's method without averaging
   :doc:`python_functions/signal_processing/welchpsd`                                     Calculate power spectral density using Welch's method
   ====================================================================================== ===========

.. toctree::
    :maxdepth: 1
    :hidden:

    python_functions/signal_processing/bartlettpsd.rst
    python_functions/signal_processing/fftfrequency.rst
    python_functions/signal_processing/filtertimeseries.rst
    python_functions/signal_processing/periodogrampsd.rst
    python_functions/signal_processing/psd2timeseries.rst
    python_functions/signal_processing/smoothsignal.rst
    python_functions/signal_processing/spectrogrampsd.rst
    python_functions/signal_processing/welchpsd.rst

Statistics
----------

.. table::
   :widths: auto
   :align: left

   ====================================================================================== ===========
   Function                                                                               Description
   ====================================================================================== ===========
   :doc:`python_functions/statistics/curvefit2d`                                          Fit curve to 2 dimensinal input dataset
   :doc:`python_functions/statistics/curvefit3d`                                          Fit curve to 3 dimensinal input dataset
   :doc:`python_functions/statistics/dataoverview`                                        Display an overview of the input data
   :doc:`python_functions/statistics/findextremum`                                        Find local extremum (minimum and maximum) in data
   :doc:`python_functions/statistics/findknn`                                             Find k-nearest neighbors using Euclidean distance
   :doc:`python_functions/statistics/fitgoodness`                                         Calculate goodness of fit parameters
   :doc:`python_functions/statistics/levelcrossing`                                       Calculate crossing point for a given level by using an upward zero crossing method
   :doc:`python_functions/statistics/movingwindow`                                        Calculate statistics of moving window through 1-d x data
   :doc:`python_functions/statistics/probability1d`                                       Calculate 1D probability density distribution for a given dataset
   :doc:`python_functions/statistics/probability2d`                                       Calculate 2D (joint) probability density distribution for two given datasets
   :doc:`python_functions/statistics/similaritymeasure`                                   Measure similarity between two arrays
   ====================================================================================== ===========

.. toctree::
    :maxdepth: 1
    :hidden:

    python_functions/statistics/curvefit2d.rst
    python_functions/statistics/curvefit3d.rst
    python_functions/statistics/dataoverview.rst
    python_functions/statistics/findextremum.rst
    python_functions/statistics/findknn.rst
    python_functions/statistics/fitgoodness.rst
    python_functions/statistics/levelcrossing.rst
    python_functions/statistics/movingwindow.rst
    python_functions/statistics/probability1d.rst
    python_functions/statistics/probability2d.rst
    python_functions/statistics/similaritymeasure.rst

SWAN Wave Model
---------------

.. table::
   :widths: auto
   :align: left

   ====================================================================================== ===========
   Function                                                                               Description
   ====================================================================================== ===========
   :doc:`python_functions/swan_wave_model/swandepthgrid`                                  Generate SWAN depth file and its associated x-y grid file
   :doc:`python_functions/swan_wave_model/swanvectorvarspconst`                           Generate SWAN file for spatially constant vector variable
   :doc:`python_functions/swan_wave_model/swanvectorvarspvariedgrid`                      Generate SWAN file for spatially varied vector variable from gridded input data
   :doc:`python_functions/swan_wave_model/swanvectorvarspvariedsct`                       Generate SWAN file for spatially varied vector variable from scattered input data
   :doc:`python_functions/swan_wave_model/swanwaterlevelspconst`                          Generate SWAN water level file for spatially constant water level
   :doc:`python_functions/swan_wave_model/swanwaterlevelspvariedgrid`                     Generate SWAN water level file for spatially varied water level from gridded input data
   :doc:`python_functions/swan_wave_model/swanwaterlevelspvariedsct`                      Generate SWAN water level file for spatially varied water level from scattered input data
   :doc:`python_functions/swan_wave_model/swanwindspconst`                                Generate SWAN wind file for spatially constant wind
   :doc:`python_functions/swan_wave_model/swanwindspvariedgrid`                           Generate SWAN wind file for spatially varied wind from gridded input data
   :doc:`python_functions/swan_wave_model/swanwindspvariedsct`                            Generate SWAN wind file for spatially varied wind from scattered input data
   ====================================================================================== ===========

.. toctree::
    :maxdepth: 1
    :hidden:

    python_functions/swan_wave_model/swandepthgrid.rst
    python_functions/swan_wave_model/swanvectorvarspconst.rst
    python_functions/swan_wave_model/swanvectorvarspvariedgrid.rst
    python_functions/swan_wave_model/swanvectorvarspvariedsct.rst
    python_functions/swan_wave_model/swanwaterlevelspconst.rst
    python_functions/swan_wave_model/swanwaterlevelspvariedgrid.rst
    python_functions/swan_wave_model/swanwaterlevelspvariedsct.rst
    python_functions/swan_wave_model/swanwindspconst.rst
    python_functions/swan_wave_model/swanwindspvariedgrid.rst
    python_functions/swan_wave_model/swanwindspvariedsct.rst

Wind Engineering
----------------

.. table::
   :widths: auto
   :align: left

   ====================================================================================== ===========
   Function                                                                               Description
   ====================================================================================== ===========
   :doc:`python_functions/wind_engineering/directionavg`                                  Average direction
   :doc:`python_functions/wind_engineering/smoothwind`                                    Smooth wind data using moving average window
   :doc:`python_functions/wind_engineering/surfaceroughness`                              Calculate shear velocity and surface roughness from a given velocity profile using Karimpour et al. (2012) method
   :doc:`python_functions/wind_engineering/sustainedwindduration`                         Calculate the sustained wind duration
   :doc:`python_functions/wind_engineering/windavg`                                       Average wind velocity and wind direction
   :doc:`python_functions/wind_engineering/winddrag`                                      Calculate wind drag coefficient, wind shear stress, and wind shear velocity
   :doc:`python_functions/wind_engineering/windgustfactor`                                Convert wind velocity of duration t0 to t
   :doc:`python_functions/wind_engineering/windspectrum`                                  Calculate wind spectrum with 513 frequencies
   :doc:`python_functions/wind_engineering/windspectrum2timeseries`                       Generate zero-mean wind velocity time series from a given spectrum
   :doc:`python_functions/wind_engineering/windvelz1toz2`                                 Convert wind velocity from first height, z1 (m), to second height, z2 (m), (e.g. 10 (m)) above surface
   ====================================================================================== ===========

.. toctree::
    :maxdepth: 1
    :hidden:

    python_functions/wind_engineering/directionavg.rst
    python_functions/wind_engineering/smoothwind.rst
    python_functions/wind_engineering/surfaceroughness.rst
    python_functions/wind_engineering/sustainedwindduration.rst
    python_functions/wind_engineering/windavg.rst
    python_functions/wind_engineering/winddrag.rst
    python_functions/wind_engineering/windgustfactor.rst
    python_functions/wind_engineering/windspectrum.rst
    python_functions/wind_engineering/windspectrum2timeseries.rst
    python_functions/wind_engineering/windvelz1toz2.rst
