#!/usr/bin/python
from __future__ import print_function
"""
Author: Tomas Chor
Date: 2015-08-07
-------------------------

This module works with micrometeorological data using pandas,
numpy, datetime and several other packages.

"""

def rotCoor(data, notation_defs=None):
    """
    Rotates the coordinates of wind data
    ----------------------

    Parameters
    ----------
    data: pandas DataFrame 
        the dataFrame to be rotated
    notation_defs: notations object
        a notation object to know which are the wind variables
    """
    from math import atan2, sqrt
    import numpy as np

    #-------
    # Getting the names for u, v, w
    if notation_defs==None:
        from core import get_notation
        defs=get_notation()
    else:
        defs=notation_defs
    wind_vars = [ defs.u, defs.v, defs.w ]
    #-------

    #-------
    # Definition of the coefficients used to created the rotation matrix
    wind_vec = data[wind_vars].mean().values
    m_u, m_v, m_w = wind_vec
    alpha = atan2(m_v,m_u)
    beta =-atan2(m_w, sqrt((m_u**2.)+(m_v**2.)))
    #-------

    #-------
    # Definition of rotation matrix
    DC= np.zeros((3,3))
    DC[0,0],DC[0,1],DC[0,2] = np.cos(alpha)*np.cos(beta), np.cos(beta)*np.sin(alpha),-np.sin(beta)
    DC[1,0],DC[1,1],DC[1,2] =-np.sin(alpha)          , np.cos(alpha)          , 0.
    DC[2,0],DC[2,1],DC[2,2] = np.cos(alpha)*np.sin(beta), np.sin(alpha)*np.sin(beta), np.cos(beta)
    #-------

    #-------
    # Application of rotation as a matrix product
    data[wind_vars] = np.dot(DC, data[wind_vars].values.T).T
    #-------

    return data


def trend(data, mode='linear', rule=None, window=None, how='mean', **kwargs):
    """
    Wrapper to return the trend given data. Can be achieved using a moving avg, block avg or polynomial fitting

    Parameters
    ----------
    data: pandas DataFrame 
        the dataFrame to be trended
    mode: string
        mode of average to apply. Currently {'moving', 'block', 'linear'}.
    rule: string
        pandas offset string to define the block in the block average. Default is "10min".
    window: pandas date offset string
        if moving average is chosen, this tells us the window size
    how: str, function
        how to resample in block type. Default is mean but it can be any numpy function
        that returns a float. E.g, median.
    kwargs:
        degree: int
            degree of polynomial fit to be passed to algs.fitByDate
    """
    import pandas as pd
    import algs
    import numpy as np

    mode=mode.lower().replace('_',' ').replace('-','').replace(' ','')
    if any(w==mode for w in ['moving', 'movingaverage']):
        #-------
        # performs moving average on the data with window equivalent to the rule
        if window==None:
            print('Warning: approximating window size for moving average!\nShould provide window keyword.\n')
            window = int(len(data)/len(data.resample(rule)))
        return pd.rolling_mean(data, window=window, **kwargs)
        #-------

    elif any(w==mode for w in ['block']):
        #-------
        # performs block average on the data with the window being the rule. Assumes that frequency is constant
        if rule==None:
            return data.apply(lambda x: [np.mean(x)]*len(x), axis=0)
        else:
            freq=data.index.inferred_freq
            aux = data.resample(rule, how=how, **kwargs)
            aux.loc[ data.index[-1] ] = np.nan
            return aux.resample(freq, fill_method='pad')
        #-------

    elif any(w in mode for w in ['linear', 'polynomial']):
        #-------
        # performs a polynomial fit on the data in blocks of "rule"
        return algs.fitByDate(data, rule=rule, **kwargs)
        #-------

    else:
        #-------
        # if no mode can be identified
        raise KeyError('Mode defined is not correct. Options are "moving", "linear", "polynomial" and "block".')
        #-------

def detrend(data, mode='moving average', rule=None, suffix="'", **kwargs):
    """
    Returns the detrended fluctuations of a given dataset

    Parameters
    ----------
    data: pandas.DataFrame, pandas.Series
        dataset to be detrended
    mode: string
        what method to use in order to identify the trend
    rule: pandas offset string
        the blocks for which the trends should be calculated
    suffix: string
        suffix to add to variable names after fluctuation is extracted
    """
    from scipy import signal
    import algs

    mode=algs.stripDown(mode.lower(), args='-_')
    df=data.copy()

    #-----------
    # If possible, try to use scipy's optimized functions
    if mode=='linear' and rule==None:
        df=df.apply(signal.detrend, axis=0, type=mode)

    elif mode=='block' and rule==None:
        df=df.apply(signal.detrend, axis=0, type='constant')
    #-----------

    #-----------
    # If not, use our trending function
    else:
        df=df-trend(df, mode=mode, rule=rule, **kwargs)
    #-----------

    return df.add_suffix(suffix)



def spectrum(data, frequency=10, out_index='frequency', anti_aliasing=False, outname=None):
    """
    Calculates the spectrum for a set of data

    Parameters
    ----------
    data: pandas.DataFrame
        dataframe with one (will return the spectrum) or two (will return to cross-spectrum) columns
    frequency: float
        frequency of measurement of signal to pass to numpy.fft.rfftfreq
    out_index: str
        only accepting 'frequency' as keyword
    anti_aliasing: bool
        not working at the moment
    outname: str
        name of the output column

    Returns:
    --------
    spectrum: dataframe
        whose column is the spectrum or coespectrum of the input dataframe
    """
    import numpy as np
    import pandas as pd

    N = len(data)
    co = False
    if type(data)==pd.DataFrame:
        if len(data.columns)==1:
            pass
        elif len(data.columns)==2:
            co=True
        else:
            raise Exception('Too many columns of data. Chose one (spectrum) or two (cross-spectrum)')
    else:
        raise Exception('Input has to be pandas.DataFrame')

    cols=list(data.columns)
    #---------
    # This works if there's two columns or just one column
    spec = np.fft.rfft(data.iloc[:,0])
    spec = np.conj(spec)*np.fft.rfft(data.iloc[:,-1])
    #---------

    #---------
    # Now we normalize the spectrum and calculate their frequency
    spec*= 2./(frequency*N)
    freq = np.fft.rfftfreq(len(data), d=1./frequency)
    #---------

    #---------
    # names it cross-spectrum for cross-spectrum or spectrum
    if co:
        varname='cross-spectrum_{}_{}'.format(*cols)
    else:
        varname='spectrum_{}'.format(cols[0])
    #---------

    #---------
    if varname != None:
        varname = outname
    #---------

    if out_index=='frequency':
        aux=pd.DataFrame( data={ varname : spec }, index=freq )
        aux.index.name='frequencies'
    else:
        raise NameError

    if co:
        return aux
    else:
        return aux.apply(np.real)


def bulkCorr(data):
    """
    Bulk correlation coefficient according to
    Cancelli, Dias, Chamecki. Dimensionless criteria for the production of...
    doi:10.1029/2012WR012127

    Parameters:
    -----------
    data: pandas.dataframe
        a two-columns dataframe
    """
    import numpy as np

    a, b = data.columns[0], data.columns[-1]
    cov = data.cov()
    r = cov.loc[a,b]
    r = r / (np.sqrt(cov.loc[a,a])*np.sqrt(cov.loc[b,b]))
    return r

#----------
# Definition of bulk_correlation according to
# Cancelli, Dias, Chamecki. Dimensionless criteria for the production of...
# doi:10.1029/2012WR012127
def _bulk_corr(self):
    import numpy as np
    df = self.copy()
    cov = df.cov()
    out = cov.copy()
    for c in out.columns:
        out.loc[:, c] = out.loc[:, c]/np.sqrt(cov.loc[c, c])
    for idx in out.index:
        out.loc[idx, :] = out.loc[idx, :]/np.sqrt(cov.loc[idx, idx])
    return out
import pandas as pd
pd.DataFrame.bulk_corr = _bulk_corr
del pd      # prevents pandas to being pre-loaded in data.py
#----------


def reverse_arrangement(array, points_number=None, alpha=0.05):
    '''
    Performs the reverse arrangement test
    according to Bendat and Piersol - Random Data - 4th edition, page 96

    Parameters
    ----------
    array: np.array, list, tuple, generator
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

    Still not adapted for dataframes
    '''
    import algs
    import numpy as np

    #-----------
    # Definition of function that determines the mean and variance
    def mu_var(N):
        """
        from Bendat&Piersol
        """
        mu=N*(N-1.)/4.
        variance=N*(2.*N + 5.)*(N - 1.)/72.
        return mu, variance
    #-----------

    if points_number==None:
        points_number=len(array)
        xarray=array
    elif points_number==len(array):
        xarray=array
    else:
        pts = len(array)//points_number
        xarray = []
        for j in range(0,points_number):
            xarray.append( np.mean(array[(j*pts):((j+1)*pts)]) ) # Calculo a media de cada um dos 50 intervalos

    A = []
    for i in range(len(xarray)-1):
        h = []
        for j in range(i,len(xarray)):
            if(xarray[i] > xarray[j]):
                h.append(1)
        A.append(sum(h))
    Atot = sum(A)
    N=len(xarray)
    mu,variance=mu_var(N)
    f=algs.inverse_normal_cdf(mu, np.sqrt(variance))
    phi1=1.-alpha/2.
    phi2=alpha/2.
    A1=f(phi1)
    A2=f(phi2)
    print(A1, Atot, A2)
    if A1 < Atot < A2:
        return True
    else:
        return False

