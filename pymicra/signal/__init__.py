from __future__ import absolute_import
from .pysignal import *
try:
    from .csignal import *
    _numba=True
except ImportError:
    #print('Numba not found')
    _numba=False
