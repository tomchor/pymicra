"""
Defines functions useful to generic signal data
"""
from __future__ import absolute_import, print_function, division
from .. import decorators as _decors


def mean(data, units, notation=None, inplace_units=True):
    """
    Calculates the mean of the DataFrame using Pymicra's notation

    Parameters
    ----------
    data: pandas.DataFrame
        dataset
    units: dict
        units dict
    notation: pymicra.Notation
        notation to be used
    inplace_units: bool
        whether to update units inplace or return a separate units dictionary
    """
    from .. import algs

    defs = algs.get_notation(notation)
    tomean = lambda x: defs.mean % x

    outunits = { tomean(key) : val for key, val in units.items() }
    if inplace_units:
        units.update(outunits)
        return data.rename(columns=tomean).mean()
    else:
        return data.rename(columns=tomean).mean(), outunits


def std(data, units, notation=None, inplace_units=True):
    """
    Calculates the std of the DataFrame using Pymicra's notation

    Parameters
    ----------
    data: pandas.DataFrame
        dataset
    units: dict
        units dict
    notation: pymicra.Notation
        notation to be used
    inplace_units: bool
        whether to update units inplace or return a separate units dictionary
    """
    from .. import algs

    defs = algs.get_notation(notation)
    tostd = lambda x: defs.std % x

    outunits = { tostd(key) : val for key, val in units.items() }
    if inplace_units:
        units.update(outunits)
        return data.rename(columns=tostd).std()
    else:
        return data.rename(columns=tostd).std(), outunits



def rotate2D(data, notation=None):
    """Rotates the coordinates of wind data

    Parameters
    ----------
    data: pandas DataFrame 
        the dataFrame to be rotated
    notation: notation object
        a notation object to know which are the wind variables

    Returns
    -------
    pandas.DataFrame
        the complete data input with the wind components rotated
    """
    from math import atan2, sqrt
    import numpy as np
    from .. import algs

    #-------
    # Getting the names for u, v, w
    defs = algs.get_notation(notation)
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


@_decors.pdgeneral(convert_out=True)
def trend(data, how='linear', rule=None, window=1200, block_func='mean', 
            center=True, units=None, ignore=None, **kwargs):
    """
    Wrapper to return the trend given data. Can be achieved using a moving avg, block avg or polynomial fitting

    Parameters
    ----------
    data: pandas.DataFrame or pandas.Series
        the data whose trend wee seek.
    how: string
        how of average to apply. Currently {'movingmean', 'movingmedian', 'block', 'linear'}.
    rule: string
        pandas offset string to define the block in the block average. Default is "10min".
    window: pandas date offset string or int
        if moving mean/median is chosen, this tells us the window size to pass to pandas. If int,
        this is the number of points used in the window. If string we will to guess the number of
        points from the index.
        Small windows (equivalent to 1min approx) work better when using rollingmedian.
    block_func: str, function
        how to resample in block type. Default is mean but it can be any numpy function
        that returns a float. E.g, median.
    degree: int
        degree of polynomial fit (only if how=='linear' or how=='polynomial')
    ignore: list
        not yet implemented
    units: units_dict
        not yet implemented

    Returns
    -------
    pandas.DataFrame or pandas.Series
        trends of data input
    """
    import pandas as pd
    from .. import algs
    import numpy as np

    how=algs.stripDown(how.lower(), args='-_')

    if ('moving' in how) or ('rolling' in how):
        if isinstance(window, str):
            window = int(len(data)/len(data.resample(window)))
        elif isinstance(window, int):
            pass
        else:
            raise TypeError('Window of moving function should be either an int or a pandas datetime offset string.')
            
        #-------
        # Performs moving average on the data with window
        if ('mean' in how) or ('average' in how):
            return pd.rolling_mean(data, window=window, center=center, **kwargs)
        #-------

        #-------
        # Performs moving median on the data with window
        elif 'median' in how:
            return pd.rolling_median(data, window=window, center=center, **kwargs)
        #-------

        #-------
        # In case how not found.
        else:
            raise KeyError('Method of trending not found. Check how keyword options with help(trend).')
        #-------

    elif any(w==how for w in ['block', 'blockaverage', 'blockmean']):
        #-------
        # performs block average on the data with the window being the rule. Assumes that frequency is constant
        if rule==None:
            if isinstance(data, pd.DataFrame):
                return data.apply(lambda x: [np.mean(x)]*len(x), axis=0)
            elif isinstance(data, pd.Series):
                return data.apply(lambda x: np.mean(x))
        else:
            freq=data.index.inferred_freq
            aux = data.resample(rule, how=block_func, **kwargs)
            aux.loc[ data.index[-1] ] = np.nan
            return aux.resample(freq, fill_method='pad')
        #-------

    elif any(w in how for w in ['linear', 'polynomial', 'poly']):
        #-------
        # performs a polynomial fit on the data in blocks of "rule"
        return data.polyfit(rule=rule, **kwargs)
        #-------

    else:
        #-------
        # if no how can be identified
        raise KeyError('Method of trending not found. Check how keyword options with help(trend).')
        #-------

@_decors.pdgeneral(convert_out=True)
def detrend(data, how='linear', rule=None, notation=None, suffix=None, units=None, inplace=True, ignore=[], **kwargs):
    """
    Returns the detrended fluctuations of a given dataset

    Parameters
    ----------
    data: pandas.DataFrame, pandas.Series
        dataset to be detrended
    how: string
        how of average to apply. Currently {'movingmean', 'movingmedian', 'block', 'linear', 'poly'}.
    rule: pandas offset string
        the blocks for which the trends should be calculated in the block and linear type
    window: pandas date offset string or int
        if moving mean/median is chosen, this tells us the window size to pass to pandas. If int,
        this is the number of points used in the window. If string we will to guess the number of
        points from the index.
        Small windows (equivalent to 1min approx) work better when using rollingmedian.
    block_func: str, function
        how to resample in block type. Default is mean but it can be any numpy function
        that returns a float. E.g, median.
    degree: int
        degree of polynomial fit (only if how=='linear' or how=='polynomial')

    Returns
    -------
    pandas.DataFrame or pandas.Series
        fluctuations of the input data
    """
    from scipy import signal
    from .. import algs
    import pandas as pd

    how=algs.stripDown(how.lower(), args='-_')
    df=data.copy()
    if ignore:
        df = df.drop(ignore, axis=1)
    defs = algs.get_notation(notation)

    #-----------
    # We can only use scipy's detrend function safely if index in not datetime
    if isinstance(df.index, pd.DatetimeIndex):
        df = df - trend(df, how=how, rule=rule, **kwargs)
    else:
        #-----------
        # If possible, try to use scipy's optimized functions
        if how=='linear' and rule==None:
            try:
                df = df.apply(signal.detrend, axis=0, type=how)
            except:
                df = df - trend(df, how=how, rule=rule, **kwargs)

        elif (any(w==how for w in ['block', 'blockaverage', 'blockmean'])) and (rule==None):
            try:
                df = df.apply(signal.detrend, axis=0, type='constant')
            except:
                df=df-trend(df, how=how, rule=rule, **kwargs)
        #-----------

        #-----------
        # If not, use our trending function
        else:
            df = df - trend(df, how=how, rule=rule, **kwargs)
        #-----------
    #-----------

    #-----------
    # We rename the columns names to indicate that they are fluctuations
    if units:
        newunits = { defs.fluctuations % el : units[el] if el in units.keys() else None for el in df.columns }
    #-----------
    
    #-----------
    # Rename the columns
    if suffix != '':
        df = df.rename(columns = lambda x: defs.fluctuations % x)
    #-----------

    #-----------
    # If units are not be changed in place, we copy them
    if units:
        if not inplace:
            units = units.copy()
        units.update(newunits)
    #-----------

    if inplace:
        return df
    else:
        return df, units


def crossSpectra(data, frequency=10, notation=None, anti_aliasing=True):
    """
    Calculates the cross-spectra for a set of data

    Parameters
    ----------
    data: pandas.DataFrame or pandas.Series
        dataframe with one (will return the spectrum) or two (will return to cross-spectrum) columns
    frequency: float
        frequency of measurement of signal to pass to numpy.fft.rfftfreq
    anti_aliasing: bool
        whether or not to apply anti-aliasing according to Gobbi, Chamecki & Dias, 2006 (doi:10.1029/2005WR004374)
    notation: notation object
        notation to be used

    Returns
    --------
    spectrum: pandas.DataFrame
        whose column is the spectrum or coespectrum of the input dataframe
    """
    from .. import notation
    from .. import algs
    import numpy as np
    import pandas as pd
    from itertools import combinations
    notation=algs.get_notation(notation)

    if len(data.columns) < 2:
        raise TypeError('DataFrame has to have more than 1 column.')

    N = len(data)
    combs = list(combinations(data.columns, 2))
    names = [ notation.cross_spectrum % (a, b) for a, b in combs ]

    specs = pd.DataFrame(columns = names)

    #---------
    # Calculate spectra here
    for (a, b) in combs: 
        spec_name = notation.cross_spectrum % (a, b)
        spec = np.fft.rfft(data.loc[:,a])
        specs.loc[:, spec_name ] = np.conj(spec) * np.fft.rfft(data.loc[:,b])
    #---------

    #---------
    # Anti-aliasing is done here
    if anti_aliasing:
        RA = np.array([ 1. + np.cos(np.pi*k/N) for k in range(N//2+1) ])/2.
        specs = specs.multiply(RA**2., axis='rows')
    #---------

    #---------
    # Now we normalize the spectrum and calculate their frequency
    specs *= 2./(frequency*N)
    freq = np.fft.rfftfreq(len(data), d=1./frequency)
    specs.index = freq
    specs.index.name='Frequency'
    #---------

    return specs




@_decors.pdgeneral(convert_out=True)
def spectra(data, frequency=10, notation=None, anti_aliasing=True):
    """
    Calculates the power spectra for a set of data

    Parameters
    ----------
    data: pandas.DataFrame or pandas.Series
        dataframe whose columns you want the power spectra of
    frequency: float
        frequency of measurement of signal to pass to numpy.fft.rfftfreq
    anti_aliasing: bool
        whether or not to apply anti-aliasing according to Gobbi, Chamecki & Dias, 2006 (doi:10.1029/2005WR004374)

    Returns
    -------
    spectra: pandas.DataFrame
        whose column is the spectrum or coespectrum of the input dataframe
    """
    from .. import notation
    from .. import algs
    import numpy as np
    import pandas as pd
    notation = algs.get_notation(notation)

    N = len(data)

    names = [ notation.spectrum % a for a in data.columns ]
    specs = pd.DataFrame(columns = names)

    #---------
    # Calculate cross-spectra here
    for name, col in zip(names, data.columns):
        spec = np.fft.rfft(data.loc[:, col ])
        specs.loc[:, name ] = np.conj(spec) * spec
    #---------

    #---------
    # Since it's the spectra, we can ignore the imaginary part
    specs = specs.apply(np.real)
    #---------

    #---------
    # Anti-aliasing is done here
    if anti_aliasing:
        RA = np.array([ 1. + np.cos(np.pi*k/N) for k in range(N//2+1) ])/2.
        specs = specs.multiply(RA**2., axis='rows')
    #---------

    #---------
    # Now we normalize the spectrum and calculate their frequency
    specs *= 2./(frequency*N)
    freq = np.fft.rfftfreq(len(data), d=1./frequency)
    specs.index = freq
    specs.index.name='Frequency'
    #---------

    return specs


def bulkCorr(data):
    """Bulk correlation coefficient according

    Bulk correlation coefficient according to 
    Cancelli, Dias, Chamecki. Dimensionless criteria for the production of...
    doi:10.1029/2012WR012127

    Parameters
    ----------
    data: pandas.dataframe
        a two-columns dataframe
    """
    import numpy as np

    a, b = data.columns[0], data.columns[-1]
    cov = data.cov()
    r = cov.loc[a,b]
    r = r / (np.sqrt(cov.loc[a,a])*np.sqrt(cov.loc[b,b]))
    return r


def test_reverse_arrangement(array, points_number=None, alpha=0.05, verbose=False):
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
    alpha: float
        Significance level for which to apply the test

    WARNING! This fuction approximates table A.6 from Bendat&Piersol as a normal distribution.
    This may no be true, since they do not express which distribution they use to construct
    their table. However, in the range 9<N<101, this approximation is as good as 5% at N=10
    and 0.1% at N=100.

    Still not adapted for dataframes
    """

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

    #-----------
    # Get reverse arrangements
    from .csignal import reverse_arrangements
    Atot = reverse_arrangements(array, points_number=points_number)
    #-----------

    #-----------
    # Define N
    if points_number==None:
        N = len(array)
    else:
        N = points_number
    #-----------

    mu, variance = mu_var(N)
    try:
        from scipy.stats import norm
        invnorm = norm(scale=np.sqrt(variance), loc=mu).isf
    except:
        invnorm = algs.inverse_normal_cdf(mu, np.sqrt(variance))

    phi1=1.-alpha/2.
    phi2=alpha/2.
    A1=invnorm(phi1)
    A2=invnorm(phi2)

    if verbose: print(A1,'<', Atot,'<', A2)

    if A1 < Atot < A2:
        return True
    else:
        return False

