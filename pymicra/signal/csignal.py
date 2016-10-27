#!/home/tomas/anaconda2/bin/python
"""
Defines functions useful to generic signal data
"""
from __future__ import print_function
from numba.pycc import CC
from numba import jit
import numpy as np
cc = CC('csignal')

def bulk_corr(x0, x1):
    """Bulk correlation coefficient according

    Bulk correlation coefficient according to 
    Cancelli, Dias, Chamecki. Dimensionless criteria for the production of...
    doi:10.1029/2012WR012127

    Parameters
    ----------
    x0, x1: np.array
        arrays of data for which ti get the bulk correlation coefficient
    """

    cov = np.cov(x0, x1)
    cov0 = cov[0, 0]
    cov1 = cov[1, 1]
    r = cov[0, 1]
    r = r / (np.sqrt(cov0)*np.sqrt(cov1))
    return r

@cc.export('reverse_arrangements', 'i4(f8[:], i4)')
def reverse_arrangements(array, points_number=0):
    """
    Performs the reverse arrangement test
    according to Bendat and Piersol - Random Data - 4th edition, page 96

    Parameters
    ----------
    array: np.array, list, tuple, generator
        array which to test for the reverse arrangement test
    points_number: integer
        number of chunks to consider to the test. Maximum is the length of the array.
        If it is less, then the number of points will be reduced by application of a mean

    WARNING! This fuction approximates table A.6 from Bendat&Piersol as a normal distribution.
    This may no be true, since they do not express which distribution they use to construct
    their table. However, in the range 9<N<101, this approximation is as good as 5% at N=10
    and 0.1% at N=100.

    Still not adapted for dataframes
    """

    #-----------
    # If number of points N is provided, we turn the run into a N-length array
    if points_number==0:
        points_number=len(array)

    if points_number==len(array):
        xarray=array
    else:
        chunklen = len(array)//points_number
        xarray = np.zeros(points_number, dtype=np.float64)
        for j in range(0, points_number):
            xarray[j] = np.mean(array[ (j*chunklen) : ((j+1)*chunklen) ])
    #-----------

    #-----------
    # We calculate the reverse arrangements
    A = 0
    for i in range(points_number - 1):
        hh = 0
        for j in range(i,len(xarray)):
            if (xarray[i] > xarray[j]):
                hh += 1
        A += hh
    #-----------

    return A


@cc.export('cvariogram', 'f8[:](f8[:], f8[:], i4)')
def cvariogram(a, b, endpoint=0):
    """
    Calculates the structure function for two arrays.

    From Stull, page 300, Eq. (8.3.1a)
    """
    N = len(a)
    if endpoint==0:
        endpoint = N//5

    if N != len(b):
        raise ValueError

    vario = np.zeros(endpoint)
    for j in range(0, endpoint):
        diffs = np.zeros(N-j)
        for k in range(0, N-j):
            diff = (a[k+j] - a[k])*(b[k+j] - b[k])
            diffs[k] = diff
        vario[j] = np.mean(diffs)
    return vario

if __name__ == '__main__':
    cc.compile()
