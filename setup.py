#!/usr/bin/env python

"""
Parts of this file were taken from Pandas under the
Apache License, Version 2.0
"""

dependencies=['Pandas', 'Numpy']

try:
	from setuptools import setup
	_has_setuptools = True
except ImportError:
	# no setuptools installed
	from distutils.core import setup
	_has_setuptools = False

if _has_setuptools:
	extra_kwargs={'install_requires' : dependencies }
else:
	import pkg_resources
	pkg_resources.require(dependencies)
	extra_kwargs={}



setup(name='pymicra',
      version='0.1.0',
      description='A Python Micrometeorology Analysis tool',
      long_description=open('README.md').read(),
      url='https://github.com/tomchor/pymicra.git',
      author='Tomas Chor',
      author_email='t.chor0@gmail.com',
      license='GNU GPL V3.0',
      packages=['pymicra'],
      zip_safe=False,
      **extra_kwargs)
