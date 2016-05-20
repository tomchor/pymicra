from __future__ import print_function
"""
This module contains functions that test certain conditions on pandas.dataframes.
They all return True for the columns that pass the test and False for the columns
that fail the test.
"""


def check_maxdif(data, tables, detrend_kw={'how':'linear', 'window':900}):
    '''
    Check the maximum and minimum differences between the fluctuations of a run.
    '''
    detrended = data.detrend(fin, **detrend_kw)
    maxdif = (detrended.max() - detrended.min()).abs()
    valid = tables.loc['dif_limits'] - maxdif
    valid = ~(valid < 0)

    return valid



def check_stationarity(data, tables, trend_kw={'how':'movingmedian', 'window':'1min'}):
    '''
    Check difference between the maximum and minimum values of the run agaisnt an upper-limit.
    This aims to flag nonstationary runs
    '''
    trend = data.trend(fin, **trend_kw)
    maxdif = (trend.max() - trend.min()).abs()
    valid = tables.loc['stat_limits'] - maxdif
    valid = ~(valid < 0)

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
        stds_list = df.resample(chunk_size, np.std).dropna()
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
 

def check_limits(data, tables, max_percent=1.):
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
    if 'low_limits' in tables.index.values:
        faulty = df < tables.loc['low_limits']
        low_count = df[ faulty ].count() 
        df[ faulty ] = np.nan
    #-------------------------------
    
    #-------------------------------
    # Now we check the upper values
    if 'upp_limits' in tables.index.values:
        faulty = df > tables.loc['upp_limits']
        upp_count = df[ faulty ].count() 
        df[ faulty ] = np.nan
    #-------------------------------

    fault_count = low_count + upp_count
    valid = fault_count < max_count

    #-------------------------------
    # Substitute faulty points by the linear trend
    trend = data.polyfit()
    df = df.fillna(trend)
    #-------------------------------

    return df, valid
 

def check_spikes(data, chunk_size='2min',
                 spikes_detrend_kw={'how':'linear'},
                 visualize=False, vis_col=1, max_consec_spikes=3,
                 cut_func=lambda x: (abs(x - x.mean()) > 4.*x.std()),
                 max_percent=1.):
    '''
    Applies spikes-check according to Vickers and Mart (1997)

    Parameters:
    -----------
    data: pandas.dataframe
        data to de-spike
    chunk_size: str, int
        size of chunks to consider. If str should be pandas offset string. If int, number of lines.
    visualize: bool
        whether of not to visualize the interpolation ocurring
    vis_col: str, int or list
        the column(s) to visualize when seeing the interpolation (only effective if visualize==True)
    max_consec_spikes: int
        maximum number of consecutive spikes to actually be considered spikes and substituted
    cut_func: function
        function used to define spikes
    max_percent: float
        maximum percentage of spikes to allow.
    '''
    import pandas as pd
    import algs
    import matplotlib.pyplot as plt

    original = data.copy()
    max_count = int(len(original)*max_percent/100.)
    dfs = algs.splitData(original, chunk_size)

    fault_count = pd.Series(len(original), index=dfs[0].columns)

    for i in range(len(dfs)):
        #-------------------------------
        # We make a copy of the original df just in case
        chunk=dfs[i].copy()
        #-------------------------------

        #-------------------------------
        # This substitutes the spikes to NaNs so it can be interpolated later
        if len(chunk)>max_consec_spikes:
            chunk=algs.limitedSubs(chunk, max_interp=max_consec_spikes, func=cut_func)
        fault_count = fault_count - chunk.count()
        #-------------------------------

        #-------------------------------
        # Substitution of spikes happens here
        trend = chunk.polyfit()
        chunk = chunk.fillna(trend)
        #-------------------------------

        #-------------------------------
        # We change the chunk in the original list of dfs to concatenate later
        dfs[i]=chunk.copy()
        #-------------------------------

    fou=pd.concat(dfs)
    valid = fault_count < max_count

    #-------------------------------
    # Visualize what you're doing to see if it's correct
    if visualize:
        print('Plotting de-spiking...')
        original[vis_col].plot(style='g-', label='original')
        fou[vis_col].plot(style='b-', label='final')
        plt.title('Column: {}'.format(vis_col))
        plt.legend()
        plt.show()
        plt.close()
    #-------------------------------

    return fou, valid


