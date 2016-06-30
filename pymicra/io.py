"""
Author: Tomas Chor
Date: 2015-08-07
-------------------------
"""
from __future__ import print_function

#-------------------------------------------
#-------------------------------------------
# INPUT OF DATA
#-------------------------------------------
#-------------------------------------------
def readDataFile(fname, varNames=None, only_named_cols=True, **kwargs):
    """
    Author: Tomas Chor

    Parameters
    ----------
    kwargs: dict
        dictionary with kwargs of pandas' read_csv function
        see http://pandas.pydata.org/pandas-docs/stable/generated/pandas.read_csv.html for more detail

    varNames: list or dict
        list or dictionary containing the names of each variable in the file (if dict, the keys must be ints)
        
    Returns
    ---------
    dataFrame: pandas.DataFrame object
    """
    import pandas as pd

    #------------
    # This makes it easier to read dates
    try:
        dtypes={ i : str for i,key in enumerate(varNames.values()) if r'%' in key }
    except:
        dtypes=None
    #------------

    #------------
    # If only_named_cols == True, read all columns in the file
    if not only_named_cols:
        usedcols = None
    else:
        usedcols = sorted(varNames.keys())
    #------------

    #------------
    # Should work, but just in case it doesn't
    try:
        data=pd.read_csv(fname, usecols=usedcols, dtype=dtypes, **kwargs)
    except ValueError:
        print('WARNING: Ignoring dtypes for date columns. This may cause problems parsing dates')
        data=pd.read_csv(fname, usecols=usedcols, **kwargs)
    #------------

    #------------
    # Renaming columns according to our variables
    data = data.rename(columns = varNames)
    #------------

    return data


def readDataFiles(flist, verbose=0, **kwargs):
    """
    Author: Tomas Chor
    Reads data from a list of files

    Parameters:
    -----------
    flist: sequence of strings
        files to be parsed
    verbose: bool
        whether to print
    **kwargs:
        readDataFile kwargs

    Returns:
    --------
    data: pandas.DataFrame
    """
    import pandas as pd

    if len(flist)==0:
        raise ValueError('Passed a list of files of zero length to be read.')
    dflist=[]
    for f in flist:
        if verbose==1:
            print('Reading',f)
        subdata=readDataFile(f, **kwargs)
        dflist.append(subdata)
    if verbose:
        print('Concatenating DataFrames...')
    data=pd.concat(dflist, ignore_index=True)
    if verbose:
        print('Done!')
    return data



def timeSeries(flist, datalogger, parse_dates=True, verbose=False,
        read_data_kw={}, parse_dates_kw={}, clean_dates=True, return_units=False, only_named_cols=True):
    """
    Creates a micrometeorological time series from a file or list of files.

    Parameters:
    -----------

    flist: list or string
        either list or names of files (dataFrame will be one concatenated dataframe) or the name of one file
    datalogger: pymicra.dataloggerConfig object
        configuration of the datalogger which is from where all the configurations of the file will be taken
    parse_date: bool
        whether or not to index the data by date. Note that if this is False many of the functionalities
        of pymicra will be lost.
        (i.d. there are repeated timestamps)
    verbose: int, bool
        verbose level
    """
    from . import algs

    #--------------
    # If datalogger is a string it should be the path to a .dlc file
    if isinstance(datalogger, str):
        datalogger = read_dlc(datalogger)
    #--------------
    
    #------------
    # We read the file(s)
    if isinstance(flist, str):
        flist=[flist]
    header_lines=datalogger.header_lines
    skiprows=datalogger.skiprows
    columns_separator=datalogger.columns_separator
    if columns_separator=='whitespace':
        timeseries=readDataFiles(flist, header=header_lines, skiprows=skiprows, delim_whitespace=True, 
            varNames=datalogger.varNames, only_named_cols=only_named_cols, **read_data_kw)
    else:
        timeseries=readDataFiles(flist, header=header_lines, skiprows=skiprows, sep=columns_separator, 
            varNames=datalogger.varNames, only_named_cols=only_named_cols, **read_data_kw)
    #------------

    #------------
    # We parse de dates
    if parse_dates:
        if verbose: print('Starting to parse the dates')
        timeseries=algs.parseDates(timeseries, dataloggerConfig=datalogger, **parse_dates_kw)
    #------------

    #------------
    # We clean the dates (if not cleaned already
    if clean_dates:
        if verbose: print('Cleaning the date columns')
        timeseries = timeseries[ [ col for col in timeseries.columns if col not in datalogger.date_col_names ] ]
    #------------

    if return_units:
        return timeseries, datalogger.units
    else:
        return timeseries


def read_dlc(dlcfile):
    """
    Reads datalogger configuration file

    WARNING! When defining the .dlc note that by default columns that are enclosed between doublequotes
    will appear without the doublequotes. So if your file is of the form :

    "2013-04-05 00:00:00", .345, .344, ...

    Then the .dlc should have: varNames=['%Y-%m-%d %H:%M:%S','u','v']. This is the default csv format of
    CampbellSci dataloggers. To disable this feature, you should parse the file with read_csv using the kw: quoting=3.
    """
    from core import dataloggerConfig

    globs={}
    dlcvars={}
    try:
        execfile(dlcfile, globs, dlcvars)
    except NameError:
        print('This version of python does not have an execfile function. This workaround should work but is yet to be fully tested')
        with open(dlcfile) as f:
            code=compile(f.read(), dlcfile, 'exec')
            exec(code, globs, dlcvars)
    return dataloggerConfig(**dlcvars)



def read_site(sitefile):
    """
    Reads .site configuration file, which holds siteConfig definitions

    The .site should have definitions as regular python syntax (in meters!):
        measurement_height  = 10
        canopy_height       = 5
        displacement_height = 3
        roughness_length    = 1.0

    sitedile: str
        path to .site file
    """
    from core import siteConfig, siteConfig

    globs={}
    sitevars={}
    try:
        execfile(sitefile, globs, sitevars)
    except NameError:
        print('This version of python does not have an execfile function. This workaround should work but is yet to be fully tested')
        with open(sitefile) as f:
            code=compile(f.read(), sitefile, 'exec')
            exec(code, globs, sitevars)

    #--------
    # First try new class, if not possible, try old one
    try:
        return siteConfig(**sitevars)
    except:
        return siteConfig(**sitevars)
    #--------


def readUnitsCsv(filename, names=0, units=1, **kwargs):
    """
    Reads a csv file in which the first line is the name of the variables
    and the second line contains the units

    Parameters:
    -----------
    filename: string
        path of the csv file to read
    names: int
        line number (starting from zero) that has the variables' names
    units: int
        line number (starting from zero) that has the variables' units

    Returns:
    --------
    df: pandas.DataFrame
        dataframe with the data
    unitsdic: dictionary
        dictionary with the variable names as keys and the units as values
    """
    import pandas as pd

    df=pd.read_csv(filename, header=[names, units], **kwargs)
    cols,units=zip(*df.columns)
    unitsdic={ k:v for k,v in zip(cols,units) }
    df.columns=cols
    return df, unitsdic


#-------------------------------------------
#-------------------------------------------
# OUTPUT OF DATA
#-------------------------------------------
#-------------------------------------------

def get_printable(data, units, to_tex_cols=True, to_tex_units=True):
    """
    Returns a csv that is pandas-printable. It does so changing the column names to add units to it.
    """
    if to_tex_cols==True:
        from constants import greek_alphabet
        columns=[ u'\\'+c if c in greek_alphabet.values() else c for c in data.columns ]
        units={ u'\\'+ c if c in greek_alphabet.values() else c : v for c,v in units.iteritems() }
    if to_tex_units==True:
        from util import printUnit as pru
        units={ k : pru(v) for k,v in units.iteritems() }
    columns=[ r'$\rm '+fl+r'\, \left({0}\right)$'.format(units[fl]) for fl in columns ]
    df=data.copy()
    df.columns=columns
    return df



