from __future__ import print_function
"""
This module contains functions that test certain conditions on pandas.dataframes.
They all return True for the columns that pass the test and False for the columns
that fail the test.
"""
from . import algs

def check_replaced(replaced, max_count=180):
    '''
    Sums and checks if the number of replaced points is larger than the
    maximum accepted
    '''
    valid = replaced < max_count

    return valid


def check_nans(data, max_percent=0.1, replace_with='interpolation'):
    '''
    Checks data for NaN values
    ''' 
    from . import data as pmdata
    df = data.copy() 
    max_count = int(len(df)*max_percent/100.)

    #-----------
    # This counts the number of NaNs
    nan_count = df.isnull().sum()
    valid = nan_count < max_count
    #-----------

    #------------
    # Replace data with either its trend or by interpolating
    if replace_with=='trend':
        trend = pmdata.trend(df, how='linear')
        df = df.fillna(trend)
    elif replace_with=='interpolation':
        df = df.interpolate(method='index', limit_direction='both')
    #------------

    return valid, nan_count


def check_maxdif(data, tables, detrend=True, detrend_kw={'how':'movingmean', 'window':900}):
    '''
    Check the maximum and minimum differences between the fluctuations of a run.
    '''
    from . import data as pmdata
    from matplotlib import pyplot as plt

    detrended = pmdata.detrend(data, suffix='', **detrend_kw)
    maxdif = (detrended.max() - detrended.min()).abs()
    valid = tables.loc['dif_limits'] - maxdif
    valid = ~(valid < 0)

    #print(valid)
    #if (~valid).any():
        #print(~valid)
        #detrended['theta_v' ].plot()
        #plt.show()

    return valid



def check_stationarity(data, tables, detrend=False,
            detrend_kw={'how':'movingmean', 'window':900}, 
            trend=True, trend_kw={'how':'movingmedian', 'window':'1min'}):
    '''
    Check difference between the maximum and minimum values of the run trend agaisnt an upper-limit.
    This aims to flag nonstationary runs
    '''
    from . import data as pmdata

    #------------
    # If detrend==True, work with the fluctuations
    if detrend:
        df = pmdata.detrend(data, suffix='', **detrend_kw)
    else:
        df = data.copy()
    #------------

    #------------
    # If trend==True, work with the trend of df (df being either absolute values or the fluctuation)
    if trend:
        trend = pmdata.trend(df, **trend_kw)
    else:
        trend = df.copy()
    #------------

    maxdif = (trend.max() - trend.min()).abs()
    valid = tables.loc['dif_limits'] - maxdif
    valid = ~(valid < 0)

    #print(valid)
    return valid

 

def check_RA(data, detrend=True, detrend_kw={'how':'linear'},
            RAT_vars=None, RAT_points=50, RAT_significance=0.05):
    '''
    Performs the Reverse Arrangement Test in each column of data

    Parameters:
    -----------
    data: pandas.DataFrame
        to apply RAT to each column
    detrend_kw: dict
        keywords to pass to pymicra.detrend
    RAT_vars:
        list of variables to which apply the RAT
    RAT_points: int
        if it's an int N, then reduce each column to N points by averaging. If None,
        then the whole columns are used
    RAT_significance: float
        significance with which to apply the RAT

    Returns:
    --------
    valid: pd.Series
        True or False for each column. If True, column passed the test
    '''
    import data as pmdata

    #-----------
    # Detrend the data
    if detrend:
        df = pmdata.detrend(data, suffix='', **detrend_kw)
    else:
        df = data.copy()
    #-----------

    #-----------
    # If RAT_vars is given, apply reverse arrangement only on these variables
    if RAT_vars:
        valid = df[RAT_vars].apply(pmdata.reverse_arrangement, axis=0, points_number=RAT_points, alpha=RAT_significance)
    elif RAT_vars==None:
        valid = df.apply(pmdata.reverse_arrangement, axis=0, points_number=RAT_points, alpha=RAT_significance)
    else:
        raise TypeError('Check RAT_vars keyword')
    #-----------

    return valid



def check_std(data, tables, detrend=False, detrend_kw={'how':'linear'}, chunk_size='2min', falseverbose=False):
    '''
    Checks dataframe for columns with too small of a standard deviation

    Parameters:
    -----------
    data: pandas.DataFrame
        dataset whose standard deviation to check 
    tables: pandas.DataFrame
        dataset containing the standard deviation limits for each column
    detrend: bool
        whether to work with the absolute series and the fluctuations
    detrend_kw: dict
        keywords to pass to pymicra.detrend with detrend==True
    chunk_size: str
        pandas datetime offset string

    Returns:
    --------
    valid: pandas.Series
        contatining True of False for each column. True means passed the test.
    '''
    import data as pmdata
    import numpy as np
    import pandas as pd
    from . import algs

    #-----------
    # Detrend the data or not
    if detrend:
        df = pmdata.detrend(data, suffix='', **detrend_kw)
    else:
        df = data.copy()
    #-----------

    #-----------
    # Separate into smaller bits or just get the full standard deviation
    if chunk_size:
        #stds_list = df.resample(chunk_size, np.std).dropna()
        stds_list = algs.resample(df, chunk_size, how=np.std).dropna()
    else:
        stds_list = pd.DataFrame(index=[df.index[0]], columns = df.columns)
        stds_list.iloc[0, :] = df.apply(np.std)
    #-----------

    #-----------
    # Check each chunk separately
    validcols = ~( stds_list < tables.loc['std_limits'] )
    falseverbose=True
    if falseverbose and (False in validcols.values):
        falsecols = [ el for el in df.columns if False in validcols.loc[:, el].values ]
        print('STD test: failed variables and times are\n{0}\n'.format(validcols.loc[:, falsecols]))
    #-----------

    valid = validcols.all(axis=0)
    return valid
 

def check_limits(data, tables, max_percent=1., replace_with='interpolation'):
    '''
    Checks dataframe for lower and upper limits. If found, they are substituted by 
    the linear trend of the run. The number of faulty points is also checked for each
    column against the maximum percentage of accepted faults max_percent

    Parameters:
    -----------
    data: pandas dataframe
        dataframe to be checked
    tables: pandas.dataframe
        dataframe with the lower and upper limits for variables
    max_percent: float
        number from 0 to 100 that represents the maximum percentage of faulty
        runs accepted by this test.

    Return:
    -------
    df: pandas.DataFrame
        input data but with the faulty points substituted by the linear trend of the run.
    valid: pandas.Series
        True for the columns that passed this test, False for the columns that didn't.
    '''
    from . import trend as pmtrend
    import numpy as np
    import algs
    import pandas as pd

    df = data.copy()
    max_count = int(len(df)*max_percent/100.)
    low_count = pd.Series(0, index=tables.columns)
    upp_count = pd.Series(0, index=tables.columns)
    fault_count = pd.Series(0, index=tables.columns)

    #-----------
    # First we check the lower values
    if 'lower_limits' in tables.index.values:
        faulty = df < tables.loc['lower_limits']
        low_count = df[ faulty ].count() 
        df[ faulty ] = np.nan
    #-------------------------------
    
    #-------------------------------
    # Now we check the upper values
    if 'upper_limits' in tables.index.values:
        faulty = df > tables.loc['upper_limits']
        upp_count = df[ faulty ].count() 
        df[ faulty ] = np.nan
    #-------------------------------

    fault_count = low_count + upp_count
    valid = fault_count < max_count

    #------------
    # Replace data with either its trend or by interpolating
    if replace_with=='trend':
        trend = pmdata.trend(df, how='linear')
        df = df.fillna(trend)
    elif replace_with=='interpolation':
        df = df.interpolate(method='index', limit_direction='both')
    #------------

    #-------------------------------
    # Substitute faulty points by the linear trend
    #trend = data.polyfit()
    #df = df.fillna(trend)
    #-------------------------------

    return df, valid, fault_count
 

def check_spikes(data, chunk_size='2min',
                 detrend=True,
                 detrend_kw={'how':'linear'},
                 visualize=False, vis_col=1, max_consec_spikes=3,
                 cut_func = lambda x: (abs(x - x.mean()) > 5.*x.std()),
                 replace_with='interpolation',
                 max_percent=1.):
    '''
    Applies spikes-check according to Vickers and Mahrt (1997)

    Parameters:
    -----------
    data: pandas.dataframe
        data to de-spike
    chunk_size: str, int
        size of chunks to consider. If str should be pandas offset string. If int, number of lines.
    detrend: bool
        whether to detrend the data and work  with the fluctuations or to work with the absolute series.
    detrend_kw: dict
        dict of keywords to pass to pymicra.trend in order to detrend data (if detrend==True).
    visualize: bool
        whether of not to visualize the interpolation ocurring
    vis_col: str, int or list
        the column(s) to visualize when seeing the interpolation (only effective if visualize==True)
    max_consec_spikes: int
        maximum number of consecutive spikes to actually be considered spikes and substituted
    cut_func: function
        function used to define spikes
    replace_with: str
        method to use when replacing spikes. Options are 'interpolation' or 'trend'.
    max_percent: float
        maximum percentage of spikes to allow.
    '''
    import pandas as pd
    import algs
    import data as pmdata

    #------------
    if replace_with=='trend':
        def replace_nans(dframe):
            trend = pmdata.trend(dframe, how='linear')
            return dframe.fillna(trend)
    elif replace_with=='interpolation':
        def replace_nans(dframe):
            return dframe.interpolate(method='index', limit_direction='both')
    #------------

    original = data.copy()

    #------------
    # If dentreded == True we save the trend for later and work with the detrended data
    if detrend:
        origtrend = pmdata.trend(data, **detrend_kw)
        detrended = original - origtrend
        dfs = algs.splitData(detrended, rule=chunk_size)
    else:
        dfs = algs.splitData(original, rule=chunk_size)
    #------------

    max_count = int(len(original)*max_percent/100.)
    fault_count = pd.Series(len(original), index=dfs[0].columns)

    for i in range(len(dfs)):
        chunk=dfs[i].copy()

        #-------------------------------
        # This substitutes the spikes to NaNs so it can be replaced later
        if len(chunk)>max_consec_spikes:
            chunk=algs.limitedSubs(chunk, max_interp=max_consec_spikes, func=cut_func)
        fault_count = fault_count - chunk.count()
        #-------------------------------

        #-------------------------------
        # Substitution of spikes happens here
        #trend = pmdata.trend(chunk, how='linear')
        #chunk = chunk.fillna(trend)
        chunk = replace_nans(chunk)
        #-------------------------------

        #-------------------------------
        # We change the chunk in the original list of dfs to concatenate later
        dfs[i]=chunk.copy()
        #-------------------------------

    #---------------------
    # Now we put the chunks back together and maybe correct the trend
    despiked = pd.concat(dfs)
    if detrend:
        fou = despiked + origtrend
    else:
        fou = despiked
    valid = fault_count < max_count
    #---------------------

    #---------------------
    # Visualize what you're doing to see if it's correct
    if visualize:
        import matplotlib.pyplot as plt
        print('Plotting de-spiking...')
        original[vis_col].plot(style='g-', label='original')
        fou[vis_col].plot(style='b-', label='final')
        plt.title('Column: {}'.format(vis_col))
        plt.legend()
        plt.show()
        plt.close()
    #---------------------

    return fou, valid, fault_count


