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
    toUnitsCsv(self, units, filename, to_tex=to_tex, **kwargs)
    return
import pandas as pd
pd.DataFrame.to_unitsCsv = _to_unitsCsv
del pd
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


