.. include:: global.rst

Getting started
===============

This "Getting started" tutorial is a brief introduction to Pymicra. This is in
no way supposed to be a complete representation of everything that can be done
with Pymicra.


Notation
--------

Pymicra uses a specific notation to name each one of its columns. This notation
is extremely important, because it is by these labels that Pymicra knows which
variable is in each column. You can check the default notation with

.. ipython:: python

    %%capture
    import pymicra as pm
    print(pm.notation)


On the left you see the full name of the variables (which corresponds to a
notation namespace) and on the right we see the default notation for that
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
you can write the configurations for some data files once, store it into a cofiguration file
and then use it from then on everytime you want to read those data files. That is what
Pymicra calls a "file configuration file", or "config file" for short. From that
file, Pymicra can create a ``pymicra.fileConfig`` object. Consider, for example, the config
file below

.. literalinclude:: ../examples/lake.config

First of all, note that the ``.config`` file is written in Python syntax, so it
has to be able to actually be run on python. This has to be true for all
``.config`` files.

Furthermore, the extension of the file does not matter. We adopt the
``.config`` extension for clarity, but it could be anything else.

The previous config file describes the data in the files in the directory
``../examples/ex_data/``. Here's an example of one such file for comparison

.. literalinclude:: ../examples/ex_data/20131108-1300.csv
   :lines: 1-10

We obtain the config object with

.. ipython:: python

   fconfig = pm.fileConfig('../examples/lake.config')
   print(fconfig)

Each variable defined in this file works as a keyword, since it can also be
input manually when calling ``pymicra.fileConfig()``. Thus, for more
information, you can also use ``help(pymicra.fileConfig``. Now we explain the
keywords one by one. In the next section we wil explain how to use this object
for reading a data file.

description
...........

The description is optional. It's a string that serves only to better identify
the config file you're dealing with. It might useful for storage purposes.

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

   fname = '../examples/ex_data/20131108-1300.csv'
   fconfig = pm.fileConfig('../examples/lake.config')
   data, units = pm.timeSeries(fname, fconfig, parse_dates=True)
   print(data)

Note that ``data`` is a ``pandas.DataFrame`` object which contains the whole data available in the
datafile with each column being a variable. Since we indicated that we wanted to parse the dates
with the option ``parse_dates=True``, each row has its respective timestamp.
If, otherwise, we were to ignore the dates, the result would be a integer-indexed dataset:

.. ipython:: python

   data2, units = pm.timeSeries(fname, fconfig, parse_dates=False)
   print(data2)

And, as mentioned, this way is a lot faster:

.. ipython:: python

   %timeit pm.timeSeries(fname, fconfig, parse_dates=False)
   %timeit pm.timeSeries(fname, fconfig, parse_dates=True)


Viewing and manipulating data
-----------------------------

To view and manipulate data, mostly you have to follow Pandas's DataFrame rules. For
that we suggest that the user visit a Pandas tutorial. However, I'll explain
some main ideas here for the sake of completeness and introduce some few ideas specific
for Pymicra that don't exist for general Pandas DataFrames.

First, for viewing raw data on screen, printing should be enough. Slicing and indexing
are supported:

.. ipython:: python

   print(data['u'])
   print(data['theta_v'])
   print(data[['u', 'v', 'w']])
   print(data['2013-11-08 13:15:00'])
