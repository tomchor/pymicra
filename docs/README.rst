.. _Github page: https://github.com/tomchor/pymicra
.. _documentation: http://pymicra.readthedocs.io


Pymicra - A Python tool for Micrometeorological Analyses
========================================================

Pymicra is a Python package designed to make working with micrometeorological
datasets more efficient. It is aimed at improving productivity (by allowing us
to focus more on micrometeorology) while still being flexible enough to let us
program project-specific things.

Please check out the `Github page`_ and the documentation_.

Here's a quick (incomplete!) list of what Pymicra does:

-  Reading, separating and understanding micrometeorological data in
   virtually any column-separated ASCII format (thanks to pandas).
-  Quality control methods (max and min values check, spikes,
   reverse-arrangement test and etc..
-  Rotation of coordinates (2D).
-  Detrending of data in the most common ways (block averages, moving
   averages and polynomial detrending).
-  Correction of sensor drift.
-  Automatic calculation of most auxiliary variables based on
   measurements (air density, dry air density, etc.).
-  Calculation of spectra and cross-spectra.
-  Calculation fluxes and characteristic scales with or without WPL correction.
-  Provide common constants generally used in atmospheric sciences.
-  Plus all native features of Pandas (interpolation, resampling,
   grouping, statistical tests, slicing, handling of missing data and
   etc.).

The package is extensively (almost entirely) based on Pandas, mostly the
``pandas.DataFrame`` class. We use Pint for units control and (generally) Numpy
or Scipy for some numerical functions not contained in Pandas.

