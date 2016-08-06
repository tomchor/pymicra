.. include:: global.rst

Getting started
===============

This "Getting started" tutorial is a brief introduction to Pymicra. This is in
no way supposed to be a complete representation of everything that can be done
with Pymicra.

In this tutorial we use some example data and refer to some example python
scripts that can be download `here
<https://github.com/tomchor/pymicra/tree/master/examples>`_. These data and
scripts are from a measurement campaign in a very small island (about 20 meters
across) in a large artificial lake. At the time of these measurements the
island was almost completely immersed into about 5 cm of water. Please feel
free to explore both the example data and the example programs, as well as
modify the programs for your own learning process!


Notation
--------

Pymicra uses a specific notation to name each one of its columns. This notation
is extremely important, because it is by these labels that Pymicra knows which
variable is in each column. You can check the default notation with

.. ipython:: python

    %%capture
    import pymicra as pm
    print(pm.notation)


The output is too long to be reproduced here, but on the left you'll see the
full name of the variables (which corresponds to a notation
namespace/attribute) and on the right you'll see the default notation for that
variable.

We recommend to use the default notation for the sake of simplicity, however,
you can change Pymicra's notation at any time by altering the attributes of
``pm.notation``. For example, by default the notation for the mean is
``'%s_mean'``, and every variable follows this base notation:

.. ipython:: python

   pm.notation.mean
   pm.notation.mean_u
   pm.notation.mean_h2o_mass_concentration

To change this, you have to change the ``mean`` notation and then re-build the
whole notation with the ``build`` method:

.. ipython:: python

   pm.notation.mean = 'm_%s'
   pm.notation.build()
   pm.notation.mean_u
   pm.notation.mean_h2o_mass_concentration
   pm.notation.h2o='v'
   pm.notation.build()
   pm.notation.mean_h2o_mass_concentration

If you just want to change the notation of one variable, but not the full
notation, just don't re-build. For example:

.. ipython:: python

   pm.notation.mean_co2_mass_concentration = 'c_m'
   pm.notation.mean_co2_mass_concentration
   pm.notation.mean_h2o_mass_concentration

It is important to note that this changes the notation used throughout every
Pymicra function.  If, however, you want to use a different notation in a
specific part of the program (in one specific function for example) you can
create a Notation object and pass it to the function, such as 

.. ipython:: python
   :verbatim:

   mynotation = pm.Notation()
   mynotation.co2='c'
   mynotation.build()
   fluxes = pm.eddyCovariance(data, units, notation=mynotation) # For example

In the example above the default Pymicra notation is left untouched, and a
separate notation is defined which is then used in a Pymicra function
separately.


Creating file configurations file
---------------------------------

The easiest way to read data files is using a ``fileConfig`` object. This object holds
the configuration of the data files so you can just call this object when reading these files.
To make it easier, Pymicra prefers to read this configurations from a file. That way
you can write the configurations for some data files once, store it into a configuration file
and then use it from then on every time you want to read those data files. That is what
Pymicra calls a "file configuration file", or "config file" for short. From that
file, Pymicra can create a ``pymicra.fileConfig`` object. Consider, for example, the config
file below

.. literalinclude:: ../examples/lake.config

First of all, note that the ``.config`` file is written in Python syntax, so it
has to be able to actually be run on python. This has to be true for all
``.config`` files.

Furthermore, the extension of the file does not matter. We adopt the
``.config`` extension for clarity, but it could be anything else.

The previous config file describes the data files in the directory
``../examples/ex_data/``. Here's an example of one such file for comparison:

.. literalinclude:: ../examples/ex_data/20131108-1000.csv
   :lines: 1-10

Note that not all columns of this file are described. Columns that are not
described are also read but are discarded by default. You can change that using
``only_named_columns=False`` in the ``timeSeries`` function.

We obtain the config object with

.. ipython:: python

   fconfig = pm.fileConfig('../examples/lake.config')
   print(fconfig)

Each variable defined in this file works as a keyword, since it can also be
input manually when calling ``pymicra.fileConfig()``. Thus, for more
information, you can also use ``help(pymicra.fileConfig)``. Now we explain the
keywords one by one. In the next section we will explain how to use this object
for reading a data file.

description
...........

The description is optional. It's a string that serves only to better identify
the config file you're dealing with. It might useful for storage purposes and useful
when printing the config object.

variables
.........

The most important keyword is ``variables``. This is a python dictionary where
each key is a column and its corresponding value is the variable in that
column. Note that we are using here the default notation to indicate which
variable is in which column. If a different notation is to be used here, then
you will have to define a new notation in your program (refer back to
`Notation`_ for that).

.. note::

   From this point on, for simplicity,  we will assume that the default notation is used.


It is imperative that the columns be named accordingly. For example, measuring
|H2O| contents in mmol/m^3 is different from measuring it in g/m^3 or mg/g. The
first is a molar density (moles per volume), the second is a mass density (mass
per volume) and the third is a mass concentration (mass per mass). In the
default notation these are indicated by the names ``'mrho_h2o'``, ``'rho_h2o'``
and ``'conc_h2o'``, respectively, and Pymicra needs to know which one is which.


Columns that contain parts of the timestamp have to have their name matching
Python's `date format string directive
<https://docs.python.org/2/library/datetime.html#strftime-and-strptime-behavior>`_,
which themselves are the 1989 version default C standard format dates, which is
common in many platforms.

This is useful only in case you want to index your data by timestamp, which is
a huge advantage in some cases (check out what Pandas can do with
`timestamp-indexed data
<http://pandas.pydata.org/pandas-docs/stable/timeseries.html>`_) but Pymicra
can also work well without this. If you don't wish to work with timestamps and
want to work only by line number in each file, you can ignore these columns and
indicate that you don't want to parse dates. In fact, parsing of dates makes
Pymicra a lot slower. Reading a file parsing its dates is about 5.5 times
slower than reading the same file without parsing any dates!


units
.....

The ``units`` keyword is also very important. It tells Pymicra in which units
each variable is being measured. Units are handled by |pint|_, so for more
details on how to define the units please refer to their documentation. Suffices 
to say here that the format of the units are pretty intuitive. Some quick remarks
are

 - prefer to define units unambiguously (``'g/(m*(s**2))'`` is generally preferred to ``'g/m/s**2'``, although both will work).
 - to define that a unit is dimensionless, ``'1'`` will not work. Define it as ``'dimensionless'`` or ``'g/g'`` and so on.
 - if one variable does not have a unit (such as a sensor flag), you don't have to include that variable.
 - the keys of ``units`` **should exactly match** the values of ``variables``.


columns_separator
.................

The ``columns_separator`` keyword is what it sounds: what separates one column
from the other.  Generally it is one character, such as a comma. A special case
happens is if the columns are separated by whitespaces of varying length, or
tabs. In that case it should be ``"whitespace"``.


frequency
.........

The ``frequency`` keyword is the frequency of the data collection in Hertz.

header_lines
............

The keyword ``header_lines`` tells us which of the first lines are part of
the file header.  If there is no header then is should be ``None``. If there
are header lines than it should be a list or int. For example, if the first two
lines of the file are part of a header, it should be ``[0, 1]``. If it were the
4 first lines, ``[0, 1, 2, 3]`` (``range(4)`` would also be acceptable).

Header lines are not used by Pymicra and are therefore skipped.


filename_format
...............
The ``filename_format`` keyword tells Pymicra how the data files are named.


date_cols
.........

The ``date_cols`` keyword is optional. It is a list of integers that indicates
which of the columns are a part of the timestamp. If it's not provided, then
Pymicra will assume that columns whose names have the character "%" in them are
part of the date and will try to parse them. If the default notation is used,
this should always be true.


Reading data
------------

To read a data file or a list of data files we use the function ``timeSeries`` along with
a config file. Let us use the config file defined in the previous subsection with one of the data
file it describes:

.. ipython:: python
   :suppress:

   import pandas as pd
   pd.options.display.max_rows=20
   pd.options.display.width=100



.. ipython:: python

   fname = '../examples/ex_data/20131108-1000.csv'
   fconfig = pm.fileConfig('../examples/lake.config')
   data, units = pm.timeSeries(fname, fconfig, parse_dates=True)
   print(data)

Note that ``data`` is a ``pandas.DataFrame`` object which contains the whole
data available in the datafile with each column being a variable. Since we
indicated that we wanted to parse the dates with the option
``parse_dates=True``, each row has its respective timestamp.  If, otherwise, we
were to ignore the dates, the result would be a integer-indexed dataset:

.. ipython:: python

   data2, units = pm.timeSeries(fname, fconfig, parse_dates=False)
   print(data2)

And, as mentioned, the latter way is a lot faster:

.. ipython:: python

   %timeit pm.timeSeries(fname, fconfig, parse_dates=False)
   %timeit pm.timeSeries(fname, fconfig, parse_dates=True)


Viewing and manipulating data
-----------------------------

To view and manipulate data, mostly you have to follow Pandas's DataFrame rules. For
that we suggest that the user visit a Pandas tutorial. However, I'll explain
some main ideas here for the sake of completeness and introduce some few ideas specific
for Pymicra that don't exist for general Pandas DataFrames.

Printing and plotting
.....................

First, for viewing raw data on screen there's printing. Slicing and indexing
are supported by Pandas, but without support for units:

.. ipython:: python
   :suppress:

   import pandas as pd
   pd.options.display.max_rows = 12

.. ipython:: python

   print(data['theta_v'])
   print(data[['u', 'v', 'w']])
   print(data['20131108 10:15:00.000':'20131108 10:17:00.000'])

Note that Pandas "guesses" if the argument you pass (``'theta_v'`` or
``'2013-11-08 10:15:00'`` etc.) is a column indexer or a row indexer. To use
these unambiguously, use the ``.loc`` method as 

.. ipython:: python

    print(data.loc['2013-11-08 10:15:00':'2013-11-08 10:17:00', ['u','v','w']])

This method is actually preferred  and you can find more information on this topic `here
<http://pandas.pydata.org/pandas-docs/stable/indexing.html>`_.

To view these data with units, you can use the ``.with_units()`` method. 
The previous output would look like this using units:


.. ipython:: python

   print(data.with_units(units)['theta_v'])
   print(data.with_units(units)[['u', 'v', 'w']])
   print(data.with_units(units)['2013-11-08 10:15:00'])

.. warning::
 Note that, although this method returns
 a Pandas DataFrame, it is not meant for calculations. Currently the DataFrame it returns is meant for
 visualization purposes only! 


We can also plot the data on screen so we can view it interactively. This can be done directly
from the DataFrame with

.. ipython:: python

   from matplotlib import pyplot as plt
   @savefig uvw_plot_basics.png
   data[['u', 'v', 'w']].plot()
   plt.show()
 
Using the ``plt.show()`` command, the plot above would plot interactively. If
we had used ``plt.savefig('figure.png')`` instead, it would have saved the
figure as png. For more on plotting, you can checkout Pandas's `visualization
guide <http://pandas.pydata.org/pandas-docs/stable/visualization.html>`_ and
find out ways to make this plot look nicer, how to render it with LaTeX and
some more tricks.

Pymicra also has an ``.xplot`` method, which brings a little more options to
Pandas's ``.plot()`` method.

.. todo::

   give xplot examples


Converting units
................

You can manually convert between units using the contents from `Manipulating`
and the Pint package. But Pymicra has a very useful method to do this called ``.convert_cols`` (more exist,
but let's focus on this one).

Let's, for example, convert some units:

.. ipython:: python

   conversions = {'p':'pascal', 'mrho_h2o':'mole/m^3', 'theta_v':'kelvin'}
   print(data.convert_cols(conversions, units, inplace_units=False))

Note that the units dictionary is updated automatically if the
``inplace_units`` keyword is true. The default is false for safety reasons, but
passing this keyword as true is much simpler and compact:

.. ipython:: python

   conversions = {'theta':'kelvin', 'theta_v':'kelvin'}
   data = data.convert_cols(conversions, units, inplace_units=True)
   print(data.with_units(units))


Manipulating
............

Manipulating data is pretty intuitive with Pandas. For example

.. ipython:: python

   data['rho_air'] = data['p']/(287.058*data['theta_v'])
   print(data['rho_air'])

If, however, you're not familiar with Pandas and prefer to just stick with what
you know, you can get Numpy arrays from columns using the ``.values``
attribute:

.. ipython:: python

   P = data['p'].values
   Tv = data['theta_v'].values
   print(type(Tv))
   rho_air = P/(287.058*Tv)
   print(rho_air)
   print(type(rho_air))

Doing that you can step out of Pandas and do your own calculations using your
own Python or Numpy code. This is pretty advantageous if you have a lot of routines
that are already written in your own way.
