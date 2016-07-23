Introduction
============

Pymicra is a Python package that was created to condense many of the knowledge
of the micrometeorological community in one single fast-to-implement software
that is freely available to anyone. Because of that, I made an effort to
include a detailed docstring along with every function and class. By doing that
I hope to have made every function, class etc. pretty self-explanatory, both
within the auto-generated docs (thanks to Sphinx) and by importing the package
and running ``help(pymicra.timeSeries)``, for example. 


Since Pymicra is meant to be a community package, improvements, suggestions of
improvement, and any kind of feedback are highly appreciated. The code is
available at its |gitrepo|_ and any contact can be made through there
(possibly creating an issue) or via e-mail.

In order for the user to program fast and effectively (and to reduce the time it
takes me to write its code), Pymicra was written on top of the
`Pandas package <http://pandas.pydata.org/>`_, so that it is faster to run
the same code using Pymicra than it is running pure Python. As a consequence,
Pymicra makes extensive use of Pandas' DataFrame class, which is a very useful
2-D data structure optimized for performance and for timestamp-indexed data.

It is possible to use Pymicra without having to be familiar with Pandas, but
because Pymicra depends on Pandas, I suggest at the very least that the user
take a quick look at a Pandas tutorial
`this one for example <http://pandas.pydata.org/pandas-docs/stable/10min.html>`_
so that one can be familiarized with the many functionalities that
Pandas offers in order to take full advantage of Pymicra. 


