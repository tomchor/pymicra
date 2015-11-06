#!/usr/bin/python
"""
Author: Tomas Chor
Date: 2015-08-07
-------------------------

Modifications:
CHECKLIST
-INCLUDE QUALITY CONTROL FEATURE
-INCLUDE MAYBE DECODIFICAION OF DATA

"""

def check_spikes(dfs, maxcount=3):
    '''
    Applies spikes-check accorgin to Vickers and Mart (1997)

    Parameters:
    -----------

    dfs: list, tuple
        sequence of pandas.DataFrame objects

    Returns:
    --------
    fou: pandas.DataFrame
        the dataframe with interpolated spikes
    valid_cols: numpy.array
        array with the percentage of valid points in each column        
    '''
    valid_cols=np.zeros(len(dfs[0].columns))
    aux=pd.DataFrame()
    for i in range(len(dfs)):
        chunk=dfs[i].copy()
        chunk[ np.abs(chunk-chunk.mean()) > cut_coef*chunk.std() ] = float('nan')
        valid_cols+=np.array(chunk.count())
        chunk=fancyInterp(chunk, maxcount=maxcount)
        chunks[i]=chunk.copy()
    fou=pd.concat(chunks)
    fou_points=len(fou.index)
    valid_cols=valid_cols/fou_points
    return fou, valid_cols




def reverse_arrangement(array, points_number=None, alpha=0.05):
    '''
    Performs the reverse arrangement test
    according to Bendat and Piersol - Random Data - 4th edition, page 96

    Parameters
    ----------

    Array: np.array, list, tuple, generator
    array which to test for the reverse arrangement test

    points_number: integer
    number of chunks to consider to the test. Maximum is the length of the array.
    If it is less, then the number of points will be reduced by application of a mean

    alpha: float
    Significance level for which to apply the test

    WARNING! This fuction approximates table A.6 from Bendat&Piersol as a normal distribution.
    This may no be true, since they do not express which distribution they use to construct
    their table. However, in the range 9<N<101, this approximation is as good as 5% at N=10
    and 0.1% at N=100.
    '''
    if points_number==None:
        points_number=len(array)
        mean_Matrix=array
    elif points_number==len(array):
        mean_Matrix=array
    else:
        pts = len(array)//points_number
        mean_Matrix = []
        for j in range(0,points_number):
            mean_Matrix.append( mean(array[(j*pts):((j+1)*pts)]) ) # Calculo a media de cada um dos 50 intervalos
    A = []
    for i in range(len(mean_Matrix)):
        h = []
        for j in range(i,len(mean_Matrix)):
            if(mean_Matrix[i] > mean_Matrix[j]):
                h.append(1)
        A.append(sum(h))
    Atot = sum(A)
    N=len(mean_Matrix)
    mu=N*(N-1)/4
    variance=N*(2*N+5)*(N-1)/72
    #
    # USE BELOW FUNCTION
    def mu_var(N):
        mu=N*(N-1.)/4.
        variance=N*(2.*N+5.)*(N-1.)/72.
        return mu, variance

    f=inverse_normal_cdf(mu, math.sqrt(variance))
    phi1=1-alpha/2.
    phi2=alpha/2.
    A1=f(phi1)
    A2=f(phi2)
    if A1 < Atot < A2:
        return True
    else:
        return False

def printUnit(string, mode='L', trim=True, greek=True):
    """
    string: string
        string (unambiguous) that represents the unit
    mode: string
        {'L' : 'latex', 'P' : 'pretty', 'H' : 'html' }
    """
    try:
        from pint.unit import UnitRegistry
    except ImportError:
        raise ImportError('You must install python-pint in order to use this function')
    ur=UnitRegistry()
    ur.default_format=mode
    u=ur[string]
    u='{:~L}'.format(u)
    if trim:
        u=u[3:].strip()
    return u
