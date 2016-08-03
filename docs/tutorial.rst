Hands-on tutorial
=================

Here we give some examples of generally-used steps to obtain fluxes, spectra
and Ogives. This is just the commonly done procedures, but more can (and
sometimes should) be done dependnig on the data.

Calculating auxiliary variables and units
-----------------------------------------

After reading the data it's easy to rotate the coordinates. Currently, only the
2D rotation is implemented. But in the future more will come (you can
contribute with another one yourself!).  To rotate, use the ``pm.rotateCoor``
function:

.. ipython:: python

   fname = '../examples/ex_data/20131108-1300.csv'
   fconfig = pm.fileConfig('../examples/lake.config')
   data, units = pm.timeSeries(fname, fconfig, parse_dates=True)
   data = data.rotateCoor(how='2d')
   print(data.with_units(units).mean())

As you can see, the u, v, w components have been rotated and ``v`` and ``w``
mean are zero.   

Pymicra has a very useful function called ``preProcess`` which "expands" the
variables in the dataset by using the original variables to calculate new ones,
such as mixing ratios, moist and dry air densities etc. In the process some
units are also converted, such as temperature units, which are converted to
Kelvin, if they are in Celsius. The ``preProcess`` function has some options to
calculate some variables in determined ways and is very "verbose", so as to
show the user exactly what it is doing:

.. iython:: python

   data = pm.preProcess(data, units, expand_temperature=True, use_means=False, rho_air_from_theta_v=True, solutes=['co2'])
   print(data.with_units(units).mean())

Note that our dataset now has many other variables and that the temperatures
are now in Kelvin. Note also that you must, at this point, specify which
solutes you wish to consider. The only solute that Pymicra looks for by default
is water vapor.

Next we calculate the fluctuations. The way to do that is with the ``detrend`` function/method.

.. ipython:: python

   ddata = data.detrend(how='linear', units=units, ignore=['theta', 'p'])
   print(ddata.with_units(units).mean())

Note the ``ignore`` keyword with which you can pass which variables you don't
want fluctuations of. In our case the pressure fluctuations can't be measured
with the barometer used and the thermodynamic variable (measured with the
LI7500) are damped also because of sensor effects. This will prevent Pymicra
from inadventently using the fluctuations of these variables with these sensors
when calculating fluxes. ``theta`` fluctuations will instead be calculated with
virtual temperature fluctuations.

The ``how`` keyword supports also ``'movingmean'``, ``'movingmedian'``,
``'block'`` and ``'poly'``.  Each of those (with the exception of ``'block'``)
are passed to Pandas or Numpy (``pandas.rolling_mean``,
``pandas.rolling_median`` and ``numpy.polyfit``) and support their respective
keywords, such as ``data.detrend(how='movingmedian', units=units,
ignore=['theta', 'p']), min_periods=1)`` to determine the minimum number of
observations in window required to have a value.


Now we must expand our dataset with the fluctuations by `joining DataFrames
<http://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.join.html>`_:

.. ipython:: python

   data = data.join(ddata)
   print(data.with_units(units).mean())


Finally, we are able to calculte the fluxes and turbulent scales. For that we
use the ``eddyCovariance`` function along with the ``get_turbulent_scales``
keyword (which is true by default):

.. ipython:: python

   siteconf = pm.siteConfig('lake.site')
   results = pm.eddyCovariance(data, units, site_config=siteconf, get_turbulent_scales=True, wpl=True, solutes=['co2'])
   print(results.with_units(units).mean())

Note that once more we must list the solutes available.

Note also that the ``units`` dictionary is automatically updated at with every
function and method with the new variables created and their units! That way if
you're in doubt of which unit the outputs are comming, just check ``units``
directly or with the ``.with_units()`` method.

If you wish to change some units, Pymicra has some handy functions that convert between units using Pint.


Extracting fluxes
-----------------

Although you can extract the fluxes manually either using the DataFrame or extracting
the Numpy arrays, Pymicra has a couple of functions that come in handy.


Obtaining the spectra
---------------------

Using Numpy's fast Fourier transform implementation, Pymicra is also able to extract
spectra, co-spectra and quadratures.
