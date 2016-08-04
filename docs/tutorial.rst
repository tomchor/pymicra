Hands-on tutorial
=================

Here we give some examples of generally-used steps to obtain fluxes, spectra
and Ogives. This is just the commonly done procedures, but more can (and
sometimes should) be done depending on the data.

Calculating auxiliary variables and units
-----------------------------------------

After reading the data it's easy to rotate the coordinates. Currently, only the
2D rotation is implemented. But in the future more will come (you can
contribute with another one yourself!).  To rotate, use the ``pm.rotateCoor``
function:

.. ipython:: python
   :suppress:

   pm.notation = pm.Notation()

.. ipython:: python

   import pymicra as pm
   fname = '../examples/ex_data/20131108-1000.csv'
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

.. ipython:: python

   data = pm.preProcess(data, units, expand_temperature=True,
       use_means=False, rho_air_from_theta_v=True, solutes=['co2'])
   print(data.with_units(units).mean())

Note that our dataset now has many other variables and that the temperatures
are now in Kelvin. Note also that you must, at this point, specify which
solutes you wish to consider. The only solute that Pymicra looks for by default
is water vapor.

It is recommended that you carefully observe the output of ``preProcess`` to
see if it's doing what you think it's doing and to analyse its calculation
logic. By doing that you can have a better idea of what it does without having
to check the code, and how to prevent it from doing some calculations you don't
want. For example, if you're not happy with the way it calculates moist air
density, you can calculate it yourself by extracting Numpy arrays using the
``.values`` method:

.. code-block:: ipython

   T = data['theta'].values
   P = data['p'].values
   water_density = data['mrho_h2o'].values
   data['rho_air'] = my_rho_air_calculation(T, P, water_density)


Next we calculate the fluctuations. The way to do that is with the ``detrend``
function/method.

.. ipython:: python

   ddata = data.detrend(how='linear', units=units, ignore=['p', 'theta'])
   print(ddata.with_units(units).mean())

Note the ``ignore`` keyword with which you can pass which variables you don't
want fluctuations of. In our case both the pressure fluctuations and the
thermodynamic temperature (measured with the LI7500 and the datalogger's
internal thermometer) can't be properly measured with the sensors used. Passing
these variables as ignored will prevent Pymicra from calculating these
fluctuations, which later on prevents it from inadvertently using the
fluctuations of these variables with these sensors when calculating fluxes.
``theta`` fluctuations will instead be calculated with virtual temperature
fluctuations and pressure fluctuations are generally not used.

The ``how`` keyword supports ``'linear'``, ``'movingmean'``,
``'movingmedian'``, ``'block'`` and ``'poly'``.  Each of those (with the
exception of ``'block'`` and ``'linear'``) are passed to Pandas or Numpy
(``pandas.rolling_mean``, ``pandas.rolling_median`` and ``numpy.polyfit``) and
support their respective keywords, such as ``data.detrend(how='movingmedian',
units=units, ignore=['theta', 'p']), min_periods=1)`` to determine the minimum
number of observations in window required to have a value.


Now we must expand our dataset with the fluctuations by `joining DataFrames
<http://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.join.html>`_:

.. ipython:: python

   data = data.join(ddata)
   print(data.with_units(units).mean())



Creating a site configuration file
----------------------------------

In order o use some micrometeorological function, we need a site configuration
object, or site object. This object tells Pymicra how the instruments are organized
and how is experimental site (vegetation, roughness length etc.). This object is created
with the ``siteConfig()`` call and the easiest way is to use it with a site config file.

.. literalinclude:: ../examples/lake.site


Extracting fluxes
-----------------

Finally, we are able to calculate the fluxes and turbulent scales. For that we
use the ``eddyCovariance`` function along with the ``get_turbulent_scales``
keyword (which is true by default). Again, it is recommended that you carefully
check the output of this function before moving on:

.. ipython:: python

   siteconf = pm.siteConfig('../examples/lake.site')
   results = pm.eddyCovariance(data, units, site_config=siteconf, 
       get_turbulent_scales=True, wpl=True, solutes=['co2'])
   print(results.with_units(units).mean())

Note that once more we must list the solutes available.

Note also that the ``units`` dictionary is automatically updated at with every
function and method with the new variables created and their units! That way if
you're in doubt of which unit the outputs are coming, just check ``units``
directly or with the ``.with_units()`` method.

Check out the example to get fluxes of many files `here
<https://github.com/tomchor/pymicra/tree/master/examples/get_fluxes.py>`_ and
download the example data `here
<https://github.com/tomchor/pymicra/tree/master/examples/ex_data>`_.




Obtaining the spectra
---------------------

Using Numpy's fast Fourier transform implementation, Pymicra is also able to extract
spectra, co-spectra and quadratures.

