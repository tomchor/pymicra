#!/usr/bin/env python
from __future__ import print_function
"""
Python Micrometeorological Analysis tool - Pymicra
"""

dependencies=['Pandas>=0.13', 'Numpy']

#-----------
# Be ready to work with or without setuptools
try:
    from setuptools import setup, find_packages
    print("Setuptools installed...\n")
    _has_setuptools = True
except ImportError:
    print("No setuptools installed.\nWill try to install with distutils...\n")
    from distutils.core import setup, find_packages
    _has_setuptools = False
#-----------

#-----------
# Being sure to find every submodule within pymicra and set up the kwargs
packages = find_packages()
#-----------

#-----------
if _has_setuptools:
    extra_kwargs={'install_requires' : dependencies}
else:
    import pkg_resources
    pkg_resources.require(dependencies)
    extra_kwargs={}
#-----------


setup(name='pymicra',
      version='0.1.3',
      description='A Python Micrometeorology Analysis tool',
      long_description=open('README.md').read(),
      url='https://github.com/tomchor/pymicra.git',
      author='Tomas Chor',
      author_email='tomaschor@gmail.com',
      license='GNU GPL V3.0',
      packages=packages,
      **extra_kwargs)
