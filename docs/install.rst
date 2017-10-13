Installation
============

.. warning::
    The commands written here assume you are running a Ubuntu-based distribution of 
    Linux. Although the basic steps should be similar for all Linux distributions, you 
    should adapt the specific commands to your system in case you are using any other OS.

.. include:: global.rst

Pymicra should work with both Python 2.7 and Python 3. Most of the required
packages already come with python, such as ``datetime`` or ``os``. The packages
that may not come as default are:

-  Pandas_
-  Pint_
-  Numpy_
-  Scipy_
-  setuptools (for installation only)

.. note::
   Version 0.17.1 of Pandas is suggested, but it should work fine with any
   distribution from 0.13 up to 0.17.1. However, version 0.18 upwards is not
   currently supported because of a change in the rolling functions API.


In order to install Pymicra the ``setuptools`` Python package should be
installed beforehand. If you don't have it installed already you can install it
with ``sudo apt install python-setuptools`` or ``sudo pip install setuptools``.

Once ``setuptools`` is installed, download the package from the `Github page`_,
unpack it somewhere, then run ``setup.py`` on a terminal with the ``install``
directive. Assuming the file is downloaded into the Downloads directory:

.. code-block:: bash

    cd ~/Downloads
    unzip pymicra-master.zip
    cp pymicra-master
    sudo python setup.py install


This should successfully install Pymicra. Note that the ``pymicra-master`` may
be different depending on whether you download directly from the master branch,
the dev branch or a specific release on Github.

Although fairly general, I have tested the setup program on a limited number of
computers so far, so it is possible that an error occurs depending on the
version of some auxiliary packages you have installed. If that happens, please
contact me through email or creating a `Github issue`_ detailing your problem and I
will try to improve the setup file accordingly. Alternatively, you may also try
to manually install versions 1.11.0 of Numpy_ and 0.17.0 of Scipy_ and then
running ``sudo python setup.py install`` again.

To remove Pymicra, the easiest way is to use pip (``sudo apt install
python-pip``) with the command ``sudo pip uninstall pymicra``.

