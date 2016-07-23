Pymicra - A Python tool for Micrometeorological Analyses
========================================================

Pymicra is a Python package designed to make it easier to work with micrometeorological datasets. It is aimed at improving the productivity and allowing us to focus on the micrometeorological, rather than programming issues.

Please check out the |gitrepo|\_ and the |docs|\_.

Please check out the [Github page](https://github.com/tomchor/pymicra) and the [documentation](http://tomchor.github.io/pymicra/).

Here's a quick (incomplete!) list of what Pymicra does:

-   Reading, separating and understanding micrometeorological data in virtually any column-separated ASCII format (thanks to pandas).
-   Quality control methods (max and min values check, spikes, reverse-arrangement test and etc..
-   Rotation of coordinates (2D).
-   Detrending of data in the most common ways (block averages, moving averages and polynomial detrending).
-   Correction of sensor drift.
-   Automatic calculation of most auxiliary variables based on measurements (air density, dry air density, etc.).
-   Calculation of spectra and cross-spectra.
-   Calculation fluxes and characteristic scales with or without WPL correction.
-   Provide common constants generally used in atmospheric sciences.
-   Plus all native features of Pandas (interpolation, resampling, grouping, statistical tests, slicing, handling of missing data and etc.).

The package is extensively (almost entirely) based on Pandas, mostly the `pandas.DataFrame` class. We use Pint for units control and (generally) Numpy or Scipy for some numerical functions not contained in Pandas.
