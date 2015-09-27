# PyMicrA - Python Micrometeorology Analysis tool

This package was designed at Lemma, at the Federal University of Parana (UFPR), to make it easier to work with micrometeorological data. Currently the package is aimed towards agregating all of the functionality commonly needed to process, read, and extract fluxes and etc from micrometeorological data. This leaves Pandas in charge of the optimization. In a near future we will also focus on our own kinds of optimization.

*This package is still under construction. More will come soon.*

## Required Packages
* Pandas

## License
This package is licensed under GNU General Public License V3.0 (http://choosealicense.com/licenses/gpl-3.0/)

## Main Features
A small list of things this package is designed to do:

  - Automatically read and understand micrometeorological data in virtually any column-separated ASCII format (thanks to pandas).
  - Automatically rotate the coordinates given wind data, so that the v and w averages are null.
  - Interpolate and detrend data in the most common ways (block averages, moving averages and polynominal detrending).
  - Calculate fluxes and characteristic scales

