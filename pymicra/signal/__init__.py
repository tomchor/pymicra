from pysignal import *
try:
    from csignal import *
except ImportError:
    print('Numba not found')
