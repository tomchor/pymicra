#!/usr/bin/env python
"""
Python Micrometeorological Analysis tool - Pymicra
"""
from __future__ import print_function
import os

#-----------
# Import version from one file used by the whole module
vfile = os.path.join(os.path.dirname(os.path.realpath(__file__)),'pymicra/version')
__version__ = open(vfile, 'rt').read().strip()
#-----------

dependencies=['Pandas', 'Numpy', 'pint']

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

extra_kwargs={'install_requires' : dependencies,
            'keywords':['meteorology', 'micrometeorology', 'atmosphere', 'science']}

setup(name='pymicra',
      version = __version__,
      description='A Python tool for Micrometeorological Analyses',
      long_description=open('README.rst').read(),
      url='https://github.com/tomchor/pymicra.git',
      author='Tomas Chor',
      author_email='tomaschor@gmail.com',
      license='GNU GPL V3.0',
      packages=packages,
      package_data={'pymicra':['pymicra.pint', 'version']},
      **extra_kwargs)

