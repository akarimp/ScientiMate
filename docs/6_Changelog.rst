Changelog
=========

Version 1.1
-----------

What is new in ver 1.1:

* MATLAB functions are added
* Global relief file ETOPO1_XYZ_0_08_Deg_Matlab.mat and ETOPO1_XYZ_0_08_Deg_Python.npz are added to GitHub
* OCEANLYZ (Python Version) is now a separate package and is removed from ScientiMate (https://oceanlyz.readthedocs.io)
* Release date: 2022-04-25

Version 1.0.5
-------------

What is new in ver 1.0.5:

* File MANIFEST.in added
* Global relief file ETOPO1_Z_0_125_Deg_Python.npz added
* oceanlyz.py bug fix, changed: warning -> warnings
* zero-crossing functions bug fix, changed: m!=0 -> m!=-1
* Bug fix (Python): removed -1 from index of if autofmaxpcorr=='on': in PcorFFTFun function (2021-08-07)
* Bug fix : Add 'if tailcorrection=='jonswap' or tailcorrection=='tma':' in WaveSpectraFun and SeaSwellFun functions (2021-10-26)
* Release date: 2021-05-26

Version 1.0.0
-------------

What is new in ver 1.0.0:

* The first version of ScientiMate is released
* Release date: 2020-11-03
