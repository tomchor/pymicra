Introduction
============

\section{Presentation}

This material is meant to be a quick introduction to the Python tool for
Micrometeorological Analyses (Pymicra). Since the package is currently in an
experimental state, it is impossible for me to include 100\% of the
functions contained in Pymicra in this guide because new functions are
constantly being created. However, the most important classes, functions and
capabilities are described here and the user is encouraged to interactively use
the Pymicra and refer to the docstrings in the code to learn more.

Pymicra is a Python package that was created to condense many of the knowledge
of the micrometeorological community in one single fast-to-implement software
that is freely available to anyone. Because of that, I made an effort to include
a detailed docstring along with every function and class, so if you import the
package and run \texttt{help(pymicra.timeSeries)}, for example, you will see a
detailed description of what the function \texttt{timeSeries} does and so on for
any other functions. It is a good idea to import Pymicra using the IPython
terminal and use its autocomplete function by hitting tab after typing
\texttt{pymicra.} to see and explore the available functions.

Since Pymicra is meant to be a community package, improvements, suggestions of
improvement, and any kind of feedback are highly appreciated. The code is
available at \href{https://github.com/tomchor/pymicra}{its github page} and any
contact can be made through there (possibly creating an issue) or via e-mail.


\section{Description and suggestions}

In order for the user to program fast and effectively (and to reduce the time it
takes me to write its code), Pymicra was written on top of the
\href{http://pandas.pydata.org/}{Pandas package}, so that it is faster to run
the same code using Pymicra than it is running pure Python. As a consequence,
Pymicra makes extensive use of Pandas' DataFrame class, which is a very useful
2-D data structure optimized for performance and for timestamp-indexed data.

It is possible to use Pymicra without having to be familiar with Pandas, but
because Pymicra depends on Pandas, I suggest at the very least that the user
take a quick look at a Pandas tutorial
(\href{http://pandas.pydata.org/pandas-docs/stable/10min.html}{this one for
example}) so that one can be familiarized with the many functionalities that
Pandas offers in order to take full advantage of Pymicra. 


