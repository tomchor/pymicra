#!/usr/bin/env python
from __future__ import print_function
import os
vfile = os.path.join(os.path.dirname(os.path.realpath(__file__)),'pymicra/version')
__version__ = open(vfile, 'rt').read().strip()

"""
Python Micrometeorological Analysis tool - Pymicra
"""

dependencies=['Pandas==0.17.1', 'Numpy', 'pint==0.7.2']

#-----------
# Be ready to work with or without setuptools
try:
    from setuptools import setup, find_packages
    print("Setuptools installed. Proceeding with installation...\n")
except ImportError:
    raise ImportError("No setuptools installed!\nPlease install setuptools in order to proceed with installation. Try: sudo apt-get install python-setuptools")
#-----------

#-----------
# Being sure to find every submodule within pymicra and set up the kwargs
packages = find_packages()
#-----------

extra_kwargs={'install_requires' : dependencies}

setup(name='pymicra',
      version = __version__,
      description='A Python tool for Micrometeorological Analyses',
      long_description=open('README.md').read(),
      url='https://github.com/tomchor/pymicra.git',
      author='Tomas Chor',
      author_email='tomaschor@gmail.com',
      license='GNU GPL V3.0',
      packages=packages,
      package_data={'pymicra':['pymicra.pint', 'version']},
      **extra_kwargs)

