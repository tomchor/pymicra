#!/usr/bin/python
"""
"""


def to_UnitsCsv(data, units, filename, to_tex=False, **kwargs):
    """
    Writes s csv with the units of the variables as a second line

    Parameters:
    ----------

    data: pandas.DataFrame
        dataframe of data
    units: dict
        dictionary containing { nameOfVar : unit }
    filename: string
        path of output file
    to_tex: bool
        whether or not to convert the string of the unit to TeX format (useful for printing)
    """
    import pandas as pd

    if to_tex:
        from util import printUnit as pru
        units={ k : pru(v) for k,v in units.iteritems() }
    cols=data.columns
    unts=[ units[c] if c in units.keys() else ' ' for c in cols ]
    columns=pd.MultiIndex.from_tuples(zip(cols, unts))
    df=data.copy()
    df.columns=columns
    df.to_csv(filename, **kwargs)
    return


#---------------
# Creates a method to write to a unitsCsv
def _to_unitsCsv(self, units, filename, to_tex=False, **kwargs):
    """
    Wrapper around toUnitsCsv to create a method to print the contents of
    a dataframe plus its units into a unitsCsv file.
    
    Parameters:
    -----------
    self: dataframe
        dataframe to write
    units: dict
        dictionary with the names of each column and their unit
    filename: str
        path to which write the unitsCsv
    to_tex: bool
        whether to try and transform the units into TeX format
    kwargs:
        to be passed to pandas' method .to_csv
    """
    to_UnitsCsv(self, units, filename, to_tex=to_tex, **kwargs)
    return
import pandas as pd
pd.DataFrame.to_unitsCsv = _to_unitsCsv
#---------------


#---------------
def _as_dlc(self, outfile, dlc):
    """
    Should write a DataFrame in the exact format described by a dataloggerConfig object
    """
    df = self.copy()

    #---------------
    # Re-create date columns if they existed
    for datecol in dlc.date_col_names:
        df[ datecol ] = df.index.strftime(datecol)
    #---------------

    if isinstance(dlc.varNames, list):
        df = df[ dlc.varNames ]
    elif isinstance(dlc.varNames, dict):
        df = df[ dlc.varNames.values ]

    df.to_csv(outfile, header=dlc.header_lines, 
            sep=dlc.columns_separator, index=False, quoting=3, na_rep='nan')
#---------------


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
#----------


#---------------
def binwrapper(self, clean_index=True, **kwargs):
    """
    Method to return binned data from a dataframe using the function classbin
    """
    from . import algs
    import numpy as np
    import pandas as pd

    #----------
    # If input is a Series, make it a DataFrame
    if isinstance(self, pd.Series):
        series=True
        self = pd.DataFrame(self, index=self.index)
    else:
        series=False
    #----------

    x=np.array(self.index)#.astype(np.float64)
    out=pd.DataFrame(columns = self.columns)
    for c in self.columns:
        xsm, ysm = algs.classbin(x, self[c].astype(np.float64), **kwargs)
        out[c] = ysm
    out.index=xsm

    #----------
    # Remove rows where the index is NaN
    if clean_index:
        out = out[ np.isfinite(out.index) ]
    #----------

    if series==True:
        return out.iloc[:, 0]
    else:
        return out
pd.DataFrame.binned=binwrapper
pd.Series.binned=binwrapper
#---------------

del pd
