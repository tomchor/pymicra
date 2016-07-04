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

#---------
# Here we create a "global" unitregistry for pymicra
try:
    from pint.unit import UnitRegistry
    import os
    abspath = os.path.dirname(os.path.abspath(__file__))
    ureg = UnitRegistry()
    ureg.load_definitions(os.path.join(abspath,'pymicra.pint'))
    Q_ = ureg.Quantity
except ImportError:
    pass
#---------
import decorators

from io import timeSeries, read_dlc, read_site, readUnitsCsv
from util import qcontrol, separateFiles, correctDrift
from micro import *
from data import *
from core import *

import io
import physics
import util
import constants
from micro import spectral
import micro
import algs
import methods
__version__ = '0.2.0'

notation = Notation()
