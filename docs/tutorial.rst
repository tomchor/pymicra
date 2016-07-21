Getting started
===============


Reading datafiles
-----------------

The easiest way to read data files is using a Pymicra class called
\texttt{dataloggerConfig} along with the function \texttt{timeSeries}. First,
let's explain what the \texttt{dataloggerConfig} class does. This class holds
most of the configurations inherent to the datalogger that actually influence
the file output (such as columns separator), which variable is in each column,
variables units, frequency, file headers, etc..

You can create a \texttt{dataloggerConfig} object manually by passing each of
these informations as a keyword to the \texttt{dataloggerConfig} function, but
the easiest method is to create a \texttt{.dlc} file (meaning datalogger
configuration), store it somewhere and create a \texttt{dataloggerConfig} object
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
that case just put ``whitespace''.

The keyword \texttt{header\_lines} is an int, or list of ints saying which rows
are headers, starting from zero. So if the first two rows of the file are a
header, the \texttt{header\_lines} in this case should be \texttt{[0, 1]}.

The \tt{filename\_format}

The \tt{date\_connector}

The \tt{first\_time\_skip}

With this file ready, it is easy to read any datafile. Consider the example
below

EXAMPLE

In it we passed the path of a file to read and the \tt{dataloggerConfig} object
containing the configuration of the file. The function \tt{timeSeries} reads the
file according to the options provided and returns a DataFrame that is put into
the \tt{data} variable. Printing \tt{data}, in this case, would print the
following on screen.

DATA PRINT

Notice that every variable defined in our datalogger configuration file appears in the
data variable.



\chapter{Viewing and manipulating data}

To view and manipulate data you have to follow Pandas's DataFrame rules. For
that we suggest that the user visit a Pandas tutorial. However, I'll explain
some main ideas here for the sake of completeness.



\chapter{Extracting fluxes}


