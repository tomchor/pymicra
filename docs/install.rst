Installation
============

.. warning::
    The commands written here assume you are running a Ubuntu-based distribution of 
    Linux. Although the basic steps should be similar for all Linux distributions, you 
    should adapt the specific commands to your system in case you are using any other OS.

.. include:: global.rst

Most of the required packages already come with python. However, packages that
generally have to be manually installed beforehand are:

-  Pandas (recommended 0.17.1)
-  Pint (0.7.2)
-  Numpy
-  Scipy
-  setuptools (for installation only)


In order to install Pymicra the ``setuptools`` Python package should be
installed. If you don't have it installed already you can install it with
``sudo apt install python-setuptools`` or ``sudo pip install setuptools``.

Download the package from the |gitrepo|_ and unpack it somewhere. Then open a
terminal and move to the directory created, whose name should be ``pymicra``.
Then run ``sudo python setup.py install``. This should be enough to install the
package.


This should successfully install Pymicra. The only
requirements are Pandas and Pint. Version 0.17 of Pandas is suggested, but it
should work fine will any distribution from 0.13 up to 0.17.1 (0.18 is not
supported because of changed in the rolling functions API).

Although fairly general, I have tested the setup program on a limited number of
computers so far, so it is possible that an error occurs depending on the
version of some auxiliary packages you have installed. If that happens and the
installation fails for some reason, please contact me through email or creating
a github issue detailing your problem and I will improve the setup file.

In case errors arise, please try to manually install versions 1.11.0 of Numpy
and 0.17.0 of Scipy and then running ``sudo python setup.py install`` again. If
still installation fails, please create a |gitissue|_ stating your OS, Python
version and any relevant information.

To remove Pymicra, the easiest way is to use pip
(``sudo apt-get install python-pip``) with the command
``sudo pip uninstall pymicra``.

