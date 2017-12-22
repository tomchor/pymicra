Hands-on tutorial
=================

Here we give some examples of generally-used steps to obtain fluxes, spectra
and Ogives. This is just the commonly done procedures, but more steps and
procedures can be done (and sometimes should, depending on the data).



.. include:: qcontrol.rst


Calculating auxiliary variables and units
-----------------------------------------

After reading the data it's easy to rotate the coordinates. Currently, only the
2D rotation is implemented. But in the future more will come (you can
contribute with another one yourself!).  To rotate, use the ``pm.rotateCoor``
function:

.. ipython:: python
   :suppress:

   import pymicra as pm
   pm.notation = pm.Notation()
   import pandas as pd
   pd.options.display.max_rows = 30



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

   data = pm.micro.preProcess(data, units, expand_temperature=True,
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
density, you can calculate it yourself by extracting Numpy_ arrays using the
``.values`` method:

.. code-block:: ipython

   T = data['theta'].values
   P = data['p'].values
   water_density = data['mrho_h2o'].values
   data['rho_air'] = my_rho_air_calculation(T, P, water_density)


Next we calculate the fluctuations. The way to do that is with the ``detrend()``
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

In order to use some micrometeorological function, we need a site configuration
object, or site object. This object tells Pymicra how the instruments are organized
and how is experimental site (vegetation, roughness length etc.). This object is created
with the ``siteConfig()`` call and the easiest way is to use it with a site config file. This
is like a ``fileConfig`` file, but simpler and for the micrometeorological aspects of the site.
Consider this example of a site config file, saved at ``../examples/lake.site``:

.. literalinclude:: ../examples/lake.site


As in the case of ``fileConfig``, the description is optional, but is handy for
organization purposes and printing on screen.

The ``measurement_height`` keyword is the general measurement height to be used
in the calculation, generally taken as the sonic anemometer height (in meters).

The ``canopy_height`` keyword is what it sounds like. Should be the mean canopy
height in meters. The ``displacement_height`` is the zero-plane displacement
height. If it is not given, it'll be calculated as 2/3 of the mean canopy
height.

The ``roughness_length`` is the roughness length in meters.

The creation of the site config object is done as

.. ipython:: python

   siteconf = pm.siteConfig('../examples/lake.site')
   print(siteconf)
 

Extracting fluxes
-----------------

Finally, we are able to calculate the fluxes and turbulent scales. For that we
use the ``eddyCovariance`` function along with the ``get_turbulent_scales``
keyword (which is true by default). Again, it is recommended that you carefully
check the output of this function before moving on:

.. ipython:: python

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

The output of this file is 

.. code-block:: python

                                                    tau                   H                  Hv                                 E                  LE  \
                       <kilogram / meter / second ** 2> <watt / meter ** 2> <watt / meter ** 2> <millimole / meter ** 2 / second> <watt / meter ** 2>   
   2013-11-08 10:00:00                         0.195731           50.837195           74.537609                          7.007486          306.392623   
   2013-11-08 12:00:00                         0.070182           -7.116972           14.737946                          6.647450          290.331014   
   2013-11-08 13:00:00                         0.059115           -3.325774           12.743194                          4.875838          212.906974   
   2013-11-08 16:00:00                         0.017075          -13.847040           -4.131275                          2.961474          128.982980   
   2013-11-08 17:00:00                         0.003284           -5.557620           -2.525949                          0.931213           40.576805   
   2013-11-08 18:00:00                         0.026564           -6.632330           -2.782190                          1.184226           51.654208   
   2013-11-08 19:00:00                         0.044550          -14.860873          -11.991508                          0.925637           40.476446   
   
                                 u_star theta_v_star theta_star            mrho_h2o_star          Lo            zeta          q_star  
                       <meter / second>     <kelvin>   <kelvin> <millimole / meter ** 3>     <meter> <dimensionless> <dimensionless>  
   2013-11-08 10:00:00         0.412886     0.156686   0.106865                16.971952  -83.429641       -0.032363        0.000266  
   2013-11-08 12:00:00         0.247697     0.051834  -0.025031                26.837060  -90.992476       -0.029673        0.000423  
   2013-11-08 13:00:00         0.227652     0.048903  -0.012763                21.417958  -81.621514       -0.033080        0.000338  
   2013-11-08 16:00:00         0.122980    -0.029651  -0.099383                24.080946   39.583979        0.068209        0.000384  
   2013-11-08 17:00:00         0.053951    -0.041353  -0.090986                17.260443    5.463569        0.494182        0.000276  
   2013-11-08 18:00:00         0.153328    -0.016003  -0.038149                 7.723494  113.854892        0.023714        0.000123  
   2013-11-08 19:00:00         0.198083    -0.053132  -0.065846                 4.672974   56.961515        0.047400        0.000074  
   

Note that the date parsing when reading the file comes is handy now, although
there are computationally faster ways to do it other than setting
``parse_dates=True`` in the ``timeSeries()`` call.

Obtaining the spectra
---------------------

Using Numpy's fast Fourier transform implementation, Pymicra is also able to
extract spectra, co-spectra and quadratures.

We begin normally with out data:

.. ipython:: python

   import pymicra as pm
   fname = '../examples/ex_data/20131108-1000.csv'
   fconfig = pm.fileConfig('../examples/lake.config')
   data, units = pm.timeSeries(fname, fconfig, parse_dates=True)
   data = data.rotateCoor(how='2d')
   data = pm.micro.preProcess(data, units, expand_temperature=True,
       use_means=False, rho_air_from_theta_v=True, solutes=['co2'])
   ddata = data.detrend(how='linear', units=units, ignore=['p', 'theta'])

   spectra = pm.spectra(ddata[["q'", "theta_v'"]], frequency=20, anti_aliasing=True)
   print(spectra)


We can plot it with the raw points, but it's hard to see anything

.. ipython:: python

   @savefig spectra.png
   spectra.plot(loglog=True, style='o')
   plt.show()

The best option it to apply a binning procedure before plotting it:

.. ipython:: python

   @savefig spectra_binned.png
   spectra.binned(bins_number=100).plot(loglog=True, style='o')
   plt.show()



We can also calculate the cross-spectra

.. ipython:: python

   crspectra = pm.crossSpectra(ddata[["q'", "theta_v'", "w'"]], frequency=20, anti_aliasing=True)
   print(crspectra)
   
We can then get the cospectra and plot it's binned version (the same can be
done with the ``.quadrature()`` method). It's important to note that, although
to pass from cross-spectra to cospectra one can merely do
``data.apply(np.real)``, it's recommended to use the ``.cospectra()`` method or
the ``pm.spectral.cospectra()`` function, since this way the notation on the
resulting DataFrame will be correctly passed on, as you can see by the legend
in the resulting plot.

.. ipython:: python

   cospectra = crspectra[[r"X_q'_w'", r"X_theta_v'_w'"]].cospectrum()
   @savefig cospectra_binned.png
   cospectra.binned(bins_number=100).plot(loglog=True, style='o')
   plt.show()
 

Now we can finally obtain an Ogive with the ``pm.spectral.Ogive()`` function.
Note that we plot ``Og/Og.sum()`` (i.e. we normalize the ogive) instead of
``Og`` merely for visualization purposes, as the scale of both ogives are too
different.

.. ipython:: python

   Og = pm.spectral.Ogive(cospectra)
   @savefig ogive.png
   (Og/Og.sum()).plot(logx=True)
   plt.show()


.. note::

 Include some examples for spectral corrections.
