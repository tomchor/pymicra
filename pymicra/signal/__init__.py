from __future__ import absolute_import, print_function, division
from .pysignal import *
try:
    import numba
    _numba=True
except ImportError:
    _numba=False

if _numba:
    from .csignal import *
