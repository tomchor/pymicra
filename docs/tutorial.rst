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

.. code-block:: python
    import pymicra
    print(pymicra.notation)

On the left you see the full name of the variables (which corresponds to a notation namespace) and on the right
we see the default notation for that variable.

You can change Pymicra's notation at any time by altering the attributes of
``pymicra.notation``.

.. todo::
   Extend this section with some examples


Reading datafiles
-----------------

The easiest way to read data files is using the ``timeSeries`` function with a
``fileConfig`` object such as 

.. code-block:: python
   config = pm.fileConfig('example.config')
   data = pm.timeSeries('1H_20130504.csv', config, parse_dates=False)

First,
let's explain what the ``fileConfig`` class does. This class holds
most of the configurations inherent to the datalogger that actually influence
the file output (such as columns separator), which variable is in each column,
variables units, frequency, file headers, etc..

Creating file configurations
............................

You can create a ``fileConfig`` object manually by passing each of
these informations as a keyword to the ``fileConfig`` function, but
the easiest method is to create a \texttt{.dlc} file (meaning datalogger
configuration), store it somewhere and create a ``fileConfig`` object
by reading it every time you work with files with that same configuration.

Consider the example \texttt{.dlc} file below, which will be explained next.

\lstinputlisting{ex_itaipu.dlc}

First of all, every .dlc file is writen in Python syntax, so it has to be able
to actually be run on python. The description is optional and the
\texttt{first\_time\_skip} and \texttt{date\_connector} keywords are generally
not necessary, so you'll rarely have to use them. They can be omitted from the
file.

The most important keyword is \texttt{varNames}. Since Pymicra works with labels for its
data, its best if all the names of the variables are properly writen, preferably
following the default Pymicra notation. Let look at how to write parts of the
timestamp first.

The columns that contain parts of the date have to have their name matching
Python's
\href{https://docs.python.org/2/library/datetime.html#strftime-and-strptime-behavior}{date
format string directive}. This is useful in case you want to index your data by
timestamp, which is a huge advantage in some cases (check out what Pandas can do
with
\href{http://pandas.pydata.org/pandas-docs/stable/timeseries.html}{timestamp-indexed
data}). If you don't wish to work with timestamps and want to work only by line
number in each file, you can ignore there columns.

As for the physical quantities, it is strongly advised to follow Pymicra's
notation, which is described in Pymicra's README.md file and explained in
Chapter \ref{chap:notation}. In the above example, which follows the Pymicra
notation, u, v and w are the three wind components, \texttt{theta\_v} stands for
the virtual temperature, and \texttt{mrho\_h2o} stands for the molar density of
\H2O. If it were the mass density the name would have been \texttt{rho\_h2o},
according to the notation.

The \texttt{date\_cols}

The \texttt{frequency}, in Hertz, is the frequency of the data collection.
Useful mostly when taking the spectrum or cospectrum.

The \texttt{columns} separator is what separates one column from the other.
Generaly it is one character, such as a comma or a whitespace. A special case
happens is if the columns are separated by whitespaces of varying length. In
that case just put "whitespace".

The keyword \texttt{header\_lines} is an int, or list of ints saying which rows
are headers, starting from zero. So if the first two rows of the file are a
header, the \texttt{header\_lines} in this case should be \texttt{[0, 1]}.

The {filename\_format}

The \tt{date\_connector}

The \tt{first\_time\_skip}

With this file ready, it is easy to read any datafile. Consider the example
below

EXAMPLE

In it we passed the path of a file to read and the \tt{fileConfig} object
containing the configuration of the file. The function \tt{timeSeries} reads the
file according to the options provided and returns a DataFrame that is put into
the \tt{data} variable. Printing \tt{data}, in this case, would print the
following on screen.

DATA PRINT

Notice that every variable defined in our datalogger configuration file appears in the
data variable.


Viewing and manipulating data
-----------------------------

To view and manipulate data you have to follow Pandas's DataFrame rules. For
that we suggest that the user visit a Pandas tutorial. However, I'll explain
some main ideas here for the sake of completeness.


Converting between different units
----------------------------------

Pymicra has some handy functions that convert between units using Pint.

Extracting fluxes
-----------------

Although you can extract the fluxes manually either using the DataFrame or extracting
the Numpy arrays, Pymicra has a couple of functions that come in handy.


Obtaining the spectra
---------------------

Using Numpy's fast Fourier transform implementation, Pymicra is also able to extract
spectra, cospectra and quadratures.
