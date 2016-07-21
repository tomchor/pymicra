Pymicra - A Python tool for Micrometeorological Analyses
========================================================

This package was designed to make it easier to work with
micrometeorological data. Pymicra is currently fully written in python
and it's aimed towards aggregating all of functionality that is commonly
needed to work with micrometeorological data into one
productivity-enhancing tool.

The package is extensively (almost entirely) based on Pandas, mostly the
``pandas.DataFrame`` class. We use Pint for units control and
(generally) Numpy or Scipy for some numerical functions not contained in
Pandas.

Required Packages
-----------------

Most required packages already come with python. However, packages that
generally have to be manually installed beforehand are:

-  Pandas
-  Pint
-  Numpy
-  Scipy
-  setuptools (for installation only)

Main Features
-------------

Currently, this is a sum up of what Pymicra does:

-  Reading, separating and understanding micrometeorological data in
   virtually any column-separated ASCII format (thanks to pandas).
-  Quality control methods (max and min values check, spikes,
   reverse-arrangement test and etc..
-  Rotation of coordinates (2D).
-  Detrending of data in the most common ways (block averages, moving
   averages and polynomial detrending).
-  Correction of sensor drift.
-  Automatic calculation of most auxiliary variables based on actual
   measurements (air density, dry air density, etc.).
-  Calculation of spectra and cross-spectra.
-  Calculation fluxes and characteristic scales.
-  WPL correction.
-  Provide all the common constants generally used in micrometeorology.
-  Plus all native features of Pandas (interpolation, resampling,
   grouping, statistical tests, slicing, handling of missing data and
   etc.)


