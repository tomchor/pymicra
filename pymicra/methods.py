#!/usr/bin/python
"""
"""
from . import decorators as _decors

#---------------
# Creates a method to write to a unitsCsv
def _to_unitsCsv(self, units, filename, **kwargs):
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
    kwargs:
        to be passed to pandas' method .to_csv
    """
    import pandas as pd

    data = self.copy()
    cols = data.columns
    unts = [ str(units[c]) if c in units.keys() else ' ' for c in cols ]
    columns = pd.MultiIndex.from_tuples(zip(cols, unts))
    data = data.copy()
    data.columns = columns
    data.to_csv(filename, **kwargs)
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

    df = df[ dlc.variables.values ]

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

#--------
# Define xplot method for pandas dataframes
def _xplot(self, xcol, reverse_x=False, return_ax=False, 
            fixed_cols=[], fcols_styles=[], latexify=False, **kwargs):
    """
    A smarter way to plot things with the x axis being one of the columns. Very useful for
    comparison of models and results

    Parameters:
    -----------
    self: pandas.DataFrame
        datframe to be plotted
    xcol: str
        the name of the column you want in the x axis
    reverse_x: bool
        whether to plot -xcol instead of xcol in the x-axis
    return_ax: bool
        whther to return pyplot's axis object for the plot
    fixed_cols: list of strings
        columns to plot in every subplot (only if you use subplot=True on keywords)
    fcols_styles: list of string
        styles to use for fixed_cols
    latexify: cool
        whether to attempt to transform names of columns into latex format
    **kwargs:
        kwargs to pass to pandas.plot method

    LATEXFY IS STILL BUGGY
    """
    from . import algs

    df = self.copy()

    #-----------
    # Try to display letters in the latex mathematical environment
    if latexfy:
	df.columns = algs.latexify(df.columns)
        xcol = algs.latexify([xcol])[0]
    #-----------

    #-----------
    # If you want ti plot against -xcol instead of xcol (useful for log plots)
    if reverse_x:
        newx = '-'+str(xcol)
        df[newx] = -df[xcol]
        df = df.drop(xcol, axis=1)
        xcol = newx
    #-----------

    subcols = df.columns.drop(fixed_cols)
    #--------
    # Checks for double xlim kwarg
    if kwargs.has_key('xlim'):
        xlim = kwargs['xlim']
        kwargs.pop('xlim')
    else:
        xlim=[df[xcol].min(), df[xcol].max()]
    #--------

    #--------
    # Here we actually plot the dataFrame
    try:
        axes = df[subcols].sort_values(xcol, axis=0).plot(x=xcol, xlim=xlim, **kwargs)
    except:
        axes = df[subcols].sort(xcol).plot(x=xcol, xlim=xlim, **kwargs)
    #--------

    #-------
    # Plot the fixed cols in each of the subplots
    if fixed_cols:
        #-------
        # The code won't work without style string
        if len(fcols_styles) == len(fixed_cols):
            pass
        else:
            fcols_styles = ['b-']*len(fixed_cols)
        #-------

        for fcol, fstyle in zip(fixed_cols, fcols_styles):
            df2 = df.sort_values(fcol, axis=0)
            for ax in axes[0]:
                ax.plot(df2[xcol], df2[fcol], fstyle, label=xcol)
    #-------

    if return_ax:
        return axes
    else:
        return
pd.DataFrame.xplot = _xplot
#--------



#---------
# Definition of dataframe method to fit
@_decors.pdgeneral(convert_out=True)
def _polyfit(self, degree=1, rule=None):
    """
    This method fits an n-degree polynomial to the dataset. The index can
    be a DateTimeIndex or not

    Parameters
    ----------
    data: pd.DataFrame, pd.Series
        dataframe whose columns have to be fitted
    degree: int
        degree of the polynomial. Default is 1.
    rule: str
        pandas offside string. Ex.: "10min".
    """
    import pandas as pd
    from .algs import fitWrap
    from . import algs

    data = self.copy()
    #if isinstance(self, pd.Series):
    #    data = pd.DataFrame(data)

    #-----------
    # If rule == None it should return a list of 1
    dflist = algs.splitData(data, rule=rule)
    #-----------

    out=pd.DataFrame()
    if isinstance(data.index, pd.DatetimeIndex):
        xx=data.index.to_julian_date()
        #for data in dflist:
        #    aux=data.apply(lambda x: fitWrap(xx, x, degree=degree), axis=0)
        #    out=out.append(aux)
    else:
        xx=data.index.values

    for data in dflist:
        aux=data.apply(lambda x: fitWrap(xx, x, degree=degree), axis=0)
        out=out.append(aux)
 
    #if isinstance(self, pd.Series):
    #    return out.iloc[:, 0]
    #else:
    #    return out
    return out
pd.DataFrame.polyfit = _polyfit
pd.Series.polyfit = _polyfit
#---------

#---------
# Detrend and trend methods for Series and DataFrames
from .data import detrend
pd.DataFrame.detrend = detrend
pd.Series.detrend  = detrend

from .data import trend
pd.DataFrame.trend = trend
pd.Series.trend  = trend
#---------


#---------
# Define convert_cols method here
from . import algs as _algs
pd.DataFrame.convert_cols = _algs.convert_cols
pd.Series.convert_indexes = _algs.convert_indexes
#---------

del pd
