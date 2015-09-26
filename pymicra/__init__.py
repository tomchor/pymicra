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
from io import dataloggerConf, timeSeries, read_dlc
from micro import *
from data import *
import io
import physics
import algs
import util
import constants
__version__='0.1.0'
