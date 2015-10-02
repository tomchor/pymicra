#!/usr/bin/python
"""
Author: Tomas Chor
Date: 2015-08-07
-------------------------

Modifications:

"""
import pandas as pd


#-------------------------------------------
#-------------------------------------------
# INPUT OF DATA
#-------------------------------------------
#-------------------------------------------
def readDataFile(fname, varNames=None, **kwargs):
    """
    Author: Tomas Chor

    Parameters
    ----------
    kwargs: dict
        dictionary with kwargs of pandas' read_csv function
        see http://pandas.pydata.org/pandas-docs/stable/generated/pandas.read_csv.html for more detail

    varNames: list
        list containing the names of each variable in the file.
        
    Returns
    ---------
    dataFrame: pandas.DataFrame object
    """
    data=pd.read_csv(fname, **kwargs)
    if varNames:
        data.columns=varNames + list(data.columns[len(varNames):])
    return data


def readDataFiles(flist, **kwargs):
    """
    Author: Tomas Chor
    ** needs to be tested! **

    -------------------------------------

    * kwargs are readDataFile kwargs

    * returns one pandas.DataFrame
    """
    if len(flist)==0:
        raise ValueError('Passed a list of files of zero length to be read.')
    data=pd.DataFrame()
    for f in flist:
        subdata=readDataFile(f, **kwargs)
        data=pd.concat( [data, subdata], ignore_index=True)
    return data

def parseDates(data, date_cols, connector='-', first_time_skip=0, clean=True):
    """
    Author: Tomas Chor
    date: 2015-08-10
    This routine parses the date from a pandas DataFrame when it is divided into several columns
    ----------------------------------

    data: pandas DataFrame

    date_cols: list of strings
    A list of the names of the columns in which the date is divided
    the naming of the date columns must be in accordance with the datetime directives,
    so if the first column is only the year, its name must be `%Y` and so forth.
    see https://docs.python.org/2/library/datetime.html#strftime-and-strptime-behavior
    """
    from datetime import timedelta,datetime
    #------------------------------------
    # joins the names of the columns, which must match the datetime directive (see __doc__)
    #------------------------------------
    date_format=connector.join(date_cols)
    auxformat='%Y-%m-%d %H:%M:%S.%f'
    #-------------------------------------
    # joins the appropriate pandas columns because pandas can read only one column into datetime
    #-------------------------------------
    try:
        aux=data[date_cols[0]].astype(int).astype(str)
    except ValueError:
        aux=data[date_cols[0]].astype(str)
    for col in date_cols[1:]:
        aux+=connector + data[col].astype(int).astype(str)
    dates=pd.to_datetime(aux, format=date_format)
    #-------------------------------------
    # The next steps are there to check if there are fractions that are not expressed in the datetime convention
    # and it assumes that the lowest time period expressed is the minute
    #-------------------------------------
    first_date=dates.unique()[1]
    n_fracs=len(dates[dates.values==first_date])
    if n_fracs>1:
        print 'Warninig! I identified that each timestamp entry contains', n_fracs, 'fractions.\n\
This generally means that the data is sampled at a frequency greater than the frequency of the timestamp.\
I will then proceed to guess the fractions based of the keyword "first_time_skip" and correct the index.'
    dates=[ date.strftime(auxformat) for date in dates ]
    aux=dates[0]
    cont=first_time_skip
    for i,date in enumerate(dates):
        if date==aux:
            pass
        else:
            cont=0
            aux=date
        dates[i]=datetime.strptime(date, auxformat) + timedelta(minutes=cont/float(n_fracs))
        cont+=1
    #-------------------------------------
    # setting new dates list as the index
    #-------------------------------------
    data=data.set_index([dates])
    #-------------------------------------
    # removing the columns used to generate the date
    #-------------------------------------
    if clean:
        data=data[ [col for col in data.columns if col not in date_cols] ]
    return data


#------------------------------------
#
#------------------------------------
class dataloggerConf(object):
    """
    This class defines a specific configuration of a datalogger output file
    --------------------------

    Parameters:
    ----------

    varNames: list of strings
    should be a list of strings with the names of the variables. If the variable
    is part if the date, then it should be provided as a datetime directive,
    so if the columns is only the year, its name must be `%Y` and so forth. While
    if it is the date in YYYY/MM/DD format, it should be `%Y/%m/%d`. For more info
    see https://docs.python.org/2/library/datetime.html#strftime-and-strptime-behavior

    date_cols: list of strings
    should be the subset of varNames that corresponds to the variables that compose
    the timestamp. If it is not provided the program will try to guess by getting
    all variable names that have a percentage sign.

    date_connector: string
    generally not really necessary. It is used to join and then parse the date_cols.

    columns_separator: string
    """
    def __init__(self, varNames,
             date_cols=None,
             frequency=None,
             date_connector='-', 
             columns_separator=',',
             header_lines=None,
             first_time_skip=False,
             units=None,
             description='generic datalogger configuration file'):
        #-------------------------
        # Makes sure that units is a dictionary type
        #-------------------------
        if units is not None:
            if not isinstance(units, dict):
                raise TypeError('units should be a dictionary. Ex.: {"u" : "m/s", "v" : "m/s", "theta" : "K" }')
        self.varNames=varNames
        if date_cols:
            self.date_cols=date_cols
        else:    #tries to guess the date columns by assuming that no other columns has a % sign on their name
            self.date_cols=[ el for el in varNames if '%' in el ]
        self.frequency=frequency
        self.date_connector=date_connector
        self.columns_separator=columns_separator
        self.header_lines=header_lines
        self.first_time_skip=first_time_skip
        self.units=units
        self.description=description


def timeSeries(flist, datalogger, index_by_date=True):
    """
    Creates a micrometeorological time series from a file or list of files.

    UNDER DEVELOPMENT
    It needs a dataloggerConf object.
    """
    
    if isinstance(flist, str):
        flist=[flist]
    header_lines=datalogger.header_lines
    columns_separator=datalogger.columns_separator
    date_cols=datalogger.date_cols
    date_connector=datalogger.date_connector
    series=readDataFiles(flist, header=header_lines, sep=columns_separator, varNames=datalogger.varNames)
    series=parseDates(series, date_cols, connector=date_connector, first_time_skip=datalogger.first_time_skip, clean=True)
    return series


def to_array(data):
    """
    Returns the contents of a timeSeries into an array type
    """
    vals=zip(*data.values)
    return [data.index.to_pydatetime]+vals


def read_dlc(dlcfile):
    """
    Reads datalogger configuration file

    WARNING! When defining the .dlc note that by default columns that are enclosed between doublequotes
    will appear without the doublequotes. So if your file is of the form :

    "2013-04-05 00:00:00", .345, .344, ...

    Then the .dlc should have: varNames=['%Y-%m-%d %H:%M:%S','u','v']. This is the default csv format of
    CampbellSci dataloggers. To disable this feature, you should parse the file with read_csv using the kw: quoting=3.
    """
    globs={}
    dlcvars={}
    try:
        execfile(dlcfile, globs, dlcvars)
    except NameError:
        print '''This version of python does not have an execfile function. This workaround should work but is yet to be fully tested'''
        with open(dlcfile) as f:
            code=compile(f.read(), dlcfile, 'exec')
            exec(code, globs, dlcvars)
    return dataloggerConf(**dlcvars)


def readUnitsCsv(filename, names=0, units=1):
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
    df=pd.read_csv(filename, header=[names, units], index_col=0, parse_dates=[0])
    cols,units=zip(*df.columns)
    unitsdic={ k:v for k,v in zip(cols,units) }
    df.columns=cols
    return df, unitsdic


#-------------------------------------------
#-------------------------------------------
# OUTPUT OF DATA
#-------------------------------------------
#-------------------------------------------

def toUnitsCsv(data, units, filename, totex=False, double_col=False):
    """
    Writes s csv with the units of the variables as a second line
    """
    if totex:
        from util import printUnit as pru
        units={ k : pru(v) for k,v in units.iteritems() }
    cols=data.columns
    unts=[ units[c] for c in cols ]
    columns=pd.MultiIndex.from_tuples(zip(cols, unts))
    df=data.copy()
    df.columns=columns
    df.to_csv(filename)
    return

def get_printable(data, units, totex_cols=True, totex_units=True):
    """
    Returns a csv that is pandas-printable. It does so changing the column names to add units to it.

    """
    if totex_cols==True:
        from constants import greek_alphabet
        columns=[ u'\\'+c if c in greek_alphabet.values() else c for c in data.columns ]
        units={ u'\\'+ c if c in greek_alphabet.values() else c : v for c,v in units.iteritems() }
    if totex_units==True:
        from util import printUnit as pru
        units={ k : pru(v) for k,v in units.iteritems() }
    columns=[ r'$\rm '+fl+r'\, \left({0}\right)$'.format(units[fl]) for fl in columns ]
    df=data.copy()
    df.columns=columns
    return df


