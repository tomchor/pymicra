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


