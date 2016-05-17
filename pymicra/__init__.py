"""
 PyMicrA - Python Micrometeorology Analysis tool
 --------------------------------
 Author: Tomas Chor
 Date of start: 2015-08-15
 --------------------------------

 This package is a collection of functions aimed at analysing 
 micrometeorology data. It mainly uses Pandas and Numpy for its
 analysis so these are required packages.

 It is recommended that pint is installed for better handling the
 units inside the functions.
"""
from io import timeSeries, read_dlc, read_site, toUnitsCsv, readUnitsCsv
from util import qcontrol, separateFiles, correctDrift
from micro import *
from data import *
from core import *

import io
import physics
import util
import util2
import constants
from micro import spectral
import micro
import algs
__version__ = '0.1.4'

try:
    from pint.unit import UnitRegistry
    ureg = UnitRegistry()
except ImportError:
    pass
