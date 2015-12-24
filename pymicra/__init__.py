"""
 PyMicrA - Python Micrometeorology Analysis tool
 --------------------------------
 Author: Tomas Chor
 Date: 2015-08-15
 --------------------------------

 This package is a collection of functions aimed at analysing 
 micrometeorology data. It mainly uses Pandas and Numpy for its
 analysis so these are required packages.
"""
from io import dataloggerConf, timeSeries, read_dlc, toUnitsCsv, readUnitsCsv
from micro import *
from data import *

import io
import physics
import algs
import util
import constants
import numeric
import notation
__version__='0.1.0'

notation.set_defaults()

global _notation

class ab:
    pass
_notation=ab()
_notation.pressure='hj'
