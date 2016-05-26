#!/usr/bin/env python
from __future__ import print_function
from pymicra.__init__ import __version__
"""
Python Micrometeorological Analysis tool - Pymicra
"""

dependencies=['Pandas>=0.13', 'Numpy']

#-----------
# Be ready to work with or without setuptools
try:
    from setuptools import setup, find_packages
    print("Setuptools installed. Proceeding with installation...\n")
except ImportError:
    raise ImportError("No setuptools installed!\nPlease install setuptools in order to proceed with installation.")
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
      **extra_kwargs)

