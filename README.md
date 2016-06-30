# Pymicra - A Python tool for Micrometeorological Analyses

This package was designed to make it easier to work with micrometeorological data. Pymicra is currently fully written in python and it's aimed towards aggregating all of the commonly needed functionality work with micrometeorological data into one productivity-enhancing tool.

The package is extensively (almost entirely) based on Pandas, which itself is based on Numpy.

## Required Packages
Most required packages already come with python. However, packages that generally have to be manually installed beforehand are:

* Pandas
* Pint
* setuptools (for installation only)

## License
This package is licensed under GNU General Public License V3.0 (http://choosealicense.com/licenses/gpl-3.0/)

## Main Features
Currently, this is a sum up of what Pymicra does:

  - Reading, separating and understanding micrometeorological data in virtually any column-separated ASCII format (thanks to pandas).
  - Quality control methods (max and min values check, spikes, reverse-arrangement test and etc..
  - Rotation of coordinates (2D).
  - Detrending of data in the most common ways (block averages, moving averages and polynomial detrending).
  - Correction of sensor drift.
  - Calculation of spectra and cross-spectra.
  - Calculation fluxes and characteristic scales.
  - Apply WPL correction.
  - Provide all the common constants generally used in micrometeorology.
  - Plus all native features of Pandas (interpolation, resampling, grouping, statistical tests, slicing, handling of missing data and etc.)

## Installation
To install Pymicra we recommend to install the `setuptools` package, which can be done with `sudo apt-get install python-setuptools` depending on your Linux distribution.

Download the package an unpack in somewhere. Then open a terminal and move to the directory created, whose name should be `pymicra`. Then run `sudo python setup.py install`. This should be enough to install the package.

To remove Pymicra, the easiest way is to use pip (`sudo apt-get install python-pip`) with the command `sudo pip uninstall pymicra`.

## Suggested notation (already default)
To work with Pandas (and hence Pymicra) each column of data preferably has to have its own name, so it is important that the names the user uses for columns match the ones that Pymicra uses (although a different-from-default notation can be passed to Pymicra). The next lines describe some general quantities, followed by the name that should be used inside straight brackets, followed by the units that compose it. The names provided are the default names and will be subjected to change by the user.
 - concentration [conc] - (mass of substance)/(mass of air)
 - density [rho] - (mass of substance)/(volume of air)
 - molar density [mrho] - (molar mass of substance)/(volume of air)
 - mixing ratio [r] - (mass of substance)/(mass of dry air)

We take humidity as the only exception to the rule above, following instead for the concentration of H2O:
 - specific humidity [q] - (mass of h2o)/(mass of air)

We also assume that the Reynolds decomposition of any variable "a" takes the form
    
    a = a_mean + a_fluctuation

where the suffixes "\_mean" and "\_fluctuation", as well as other suffixes, are written within the variables' names as
 - mean [\_mean] - indicates that the mean was taken on the variable whose name precedes this
 - fluctuation ['] - indicates that this is the fluctuation of the variable: a' = a - a\_mean
 - star/asterisk [\_star] - indicates a turbulence scale of the variable: a\_star = mean(u'\*a')

So that, in the Pymicra standard notation, the Reynolds decomposition is written as

    a = a_mean + a'

The standard notation for common meteorological variables is 
 - thermodynamic temperature: [theta] 
 - virtual temperature: [theta\_v]
 - pressure: [p]

The standard notation for micrometeorological variables currently is
 - Monin-Obukhov similarity variable: [zeta]
 - Monin-Obukhov length: [Lm]
 - momentum flux: [tau]
 - sensible heat flux: [H]
 - virtual sensible heat flux: [Hv]
 - water vapor flux: [E]
 - latent heat flux: [LE]
 - solute s flux: [F\_s]

Some examples of names of variables following the standard notation are (description - [variable name in Pymicra]):
 - u fluctuations: [u']
 - turbulence scale for virtual temperature: [theta\_v\_star]
 - fluctuations of water density on air: [rho\_h2o']
 - mean mixing ratio for carbon dioxide: [r\_co2\_mean]
