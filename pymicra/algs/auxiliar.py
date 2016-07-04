from __future__ import print_function
"""
"""
#import pandas as pd
#import numpy as np

#
#def splitData(data, rule='30min', return_index=False, **kwargs):
#    """
#    Splits a given pandas DataFrame into a series of "rule"-spaced DataFrames
#
#    Parameters
#    ----------
#    data: pandas dataframe
#        data to be split
#    rule: str or int 
#        If it is a string, it should be a pandas string offset.
#            Some possible values (that should be followed by an integer) are:
#            D   calendar day frequency
#            W   weekly frequency
#            M   month end frequency
#            MS  month start frequency
#            Q   quarter end frequency
#            BQ  business quarter endfrequency
#            QS  quarter start frequency
#            A   year end frequency
#            AS  year start frequency
#            H   hourly frequency
#            T   minutely frequency
#            Min minutely frequency
#            S   secondly frequency
#            L   milliseconds
#            U   microseconds
#
#        If it is a int, it should be the number of lines desired in each separated piece.
#
#        If it is None, then the dataframe isn't separated and a list containing only the
#        full dataframe is returned.
#        
#        check it complete at http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases
#        """
#    import pandas as pd
#
#    #---------
#    # Choose how to separate the data
#    if rule == None:
#        out = [ data ]
#    elif type(rule) == int:
#        out = [ data.iloc[ rule*i : rule*(i+1) ] for i in range(0, len(data)/rule) ]
#    else:
#        from itertools import izip_longest as izip
#        #---------
#        # We first create the index in which we base our separation
#        # THIS STEP CAN PROBABLY BE IMPROVED
#        res_dates = pd.Series(index=data.index).resample(rule, **kwargs).index
#        intervals = izip(res_dates, res_dates[1:], fillvalue=data.index[-1] + pd.DateOffset(microseconds=2))
#        #---------
#
#        out = [ data.loc[ bdate:edate - pd.DateOffset(microseconds=1) ] for bdate, edate in intervals ]
#    #---------
#
#    if return_index:
#        return out, res_dates
#    else:
#        return out
#
#
#def resample(df, rule, how=None, **kwargs):
#    '''
#    Extends pandas resample methods to index made of integers
#    '''
#    import pandas as pd
#    if how==None:
#        import numpy as np
#        how = np.mean
#
#    if isinstance(df.index, pd.DatetimeIndex) and isinstance(rule, str):
#        return df.resample(rule, how, **kwargs)
#    else:
#        s = (df.index.to_series() / rule).astype(int)
#        aux = df.groupby(s).apply(how).set_index(s.index[::rule])
#        return aux
#
#
#def fitWrap(x, y, degree=1):
#    """
#    A wrapper to numpy.polyfit and numpy.polyval that fits data given an x and y arrays.
#    This is specifically designed to be used with by pandas.DataFrame.apply method
#
#    Parameters
#    ----------
#    x: array, list
#    y: array, list
#    degree: int
#    """
#    import numpy as np
#
#    x = np.array(x)
#    y = np.array(y)
#
#    #----------
#    # Getting rid of NaNs and passing only "good" points to polyfit
#    idx = np.isfinite(x) & np.isfinite(y)
#    coefs = np.polyfit(x[idx], y[idx], degree)
#    #----------
#
#    yy = np.polyval(coefs, x)
#    return yy
#
#
#def fitByDate(data, degree=1, rule=None):
#    """
#    Given a pandas DataFrame with the index as datetime, this routine
#    fit a n-degree polynomial to the dataset
#
#    Parameters
#    ----------
#    data: pd.DataFrame, pd.Series
#        dataframe whose columns have to be fitted
#    degree: int
#        degree of the polynomial. Default is 1.
#    rule: str
#        pandas offside string. Ex.: "10min".
#    """
#    import pandas as pd
#
#    #-----------
#    # If rule == None it should return a list of 1
#    dflist=splitData(data, rule=rule)
#    #-----------
#
#    out=pd.DataFrame()
#    for data in dflist:
#        xx=data.index.to_julian_date()
#        aux=data.apply(lambda x: fitWrap(xx, x, degree=degree), axis=0)
#        out=out.append(aux)
#    return out
#

def stripDown(str, final='', args=['_', '-']):
    """
    Auxiliar function to strip down keywords from symbols
    """
    for arg in args:
        str=str.replace(arg,final)
    return str
        

def lenYear(year):
    """
    Calculates the length of a year in days
    Useful to figure out if a certain year is a leap year
    """
    import calendar
    feblen=calendar.monthrange(year,2)[1]
    otherlens=365-28
    return feblen + otherlens


#def limited_interpolation(data, maxcount=3):
#    '''
#    Interpolates linearly but only if gap is smaller of equal to maxcout
#
#    Parameters:
#    -----------
#
#    data: pandas.DataFrame
#        dataset to interpolate
#    maxcount: int
#        maximum number of consecutive NaNs to interpolate. If the number is smaller than that, nothing is done with the points.
#    '''
#    mask = data.copy()
#    grp = ((mask.notnull() != mask.shift().notnull()).cumsum())
#    grp['ones'] = 1
#    for i in data.columns:
#        mask[i] = (grp.groupby(i)['ones'].transform('count') <= maxcount) | data[i].notnull()
#    data=data.interpolate(mode='time').bfill()[mask]
#    return data
#

def completeHM(string):
    if string.isdigit():
        pass
    else:
        raise TypeError('String passed must contain only digits. Check the argument')
    if len(string)==3:
        string='0'+string
    return string

#
#def limitedSubs(data, max_interp=3, func=lambda x: abs(x) > abs(x.std()*4.) ):
#    """
#    Substitute elements for NaNs if a certain conditions given by fund is met at 
#    a maximum of max_interp times in a row.
#    If there are more than that number in a row, then they are not substituted.
#
#    Parameters:
#    -----------
#    data: pandas.dataframe
#        data to be interpolated
#    max_interp: int
#        number of maximum NaNs in a row to interpolate
#    func: function
#        function of x only that determines the which elements become NaNs. Should return
#        only True or False.
#
#    Returns:
#    --------
#    df: pandas.dataframe
#        dataframe with the elements substituted
#    """
#    import numpy as np
#
#    df=data.copy()
#    cond=func(df)
#    for c in df.columns:
#        grouper = (cond[c] != cond[c].shift(1)).cumsum() * cond[c]
#        fill = df.groupby(grouper)[c].transform(lambda x: x.size) <= max_interp
#        #fill = (df.groupby(grouper)[c].transform('size') <= max_interp)
#        df.loc[fill & cond[c], c] = np.nan
#    return df
#

def testValid(df_valid, testname='', falseverbose=True, trueverbose=True, filepath=None):
    '''
    Tests a boolean DataFrane obtained from the test and prints standard output

    Parameters:
    -----------
    df_valid: pandas.Series
        series contaning only True or False values for each of the variables, which should be the indexes
    testname: string
        the name of the test that generated the True/False values
    falseverbose: bool
        whether to return which variables caused a false result
    trueverbose: bool
        whether to print something successful cases

    Returns:
    --------
    result: bool
        True if the run passed the passed
    failed: list
        list of failed variables if result==False. None otherwise.
    '''
    from os.path import basename
    
    if False in df_valid.values:
        failed = df_valid[ df_valid==False ].index
        print(basename(filepath), ': !FAILED',testname,'test!\n')
        if falseverbose:
            print('Failed variable(s):', ', '.join(failed),'\n')
            print
        return False, failed
    else:
        if trueverbose: print(basename(filepath),'passed',testname,'test')
        return True, None

#
#def file_len(fname):
#    """
#    Returns length of a file through piping bash's function wc
#
#    Parameters:
#    -----------
#
#    fname: string
#        path of the file
#    """
#    import subprocess
#    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, 
#                                              stderr=subprocess.PIPE)
#    result, err = p.communicate()
#    if p.returncode != 0:
#        raise IOError(err)
#    return int(result.strip().split()[0])
#


def applyResult(result, failed, df, control=None, testname=None, filename=None, falseshow=False, index_n=None):
    """
    Auxiliar function to be used with util.qcontrol

    Parameters:
    -----------
    result: bool
        whether the test failed and succeeded
    failed: list
        list of failed variables. None object if the test was successful
    control: dictionary
        dictionary whose keys are the names of the tests and items are lists
    testname: string
        name of the test (has to match control dict)
    filename: string
        name or path or identifier of the file tested
    falseshow: bool
        whether to show the failed variables or not
    """
    import matplotlib.pyplot as plt
    if result==False:
        if falseshow:
            df[failed].plot()
            plt.show()
        if type(control)==dict:
            control[testname].append(filename)
        else:
            control.loc[ index_n, testname ] = filename
    return control


#def inverse_normal_cdf(mu, sigma):
#    """
#    Applied the inverse normal cumulative distribution
#
#    mu: mean
#    sigma: standard deviation
#    """
#    from scipy.special import erfinv
#    import numpy as np
#
#    def f(phi):
#        Z=np.sqrt(2.)*erfinv(-2.*phi+1.)
#        return sigma*Z + mu
#    return f
#
#
#def parseDates(data, dataloggerConfig=None,
#        date_col_names=None, connector='-', first_time_skip=0,
#        clean=True, correct_fracs=None, complete_zeroes=False, verbose=False):
#    """
#    Author: Tomas Chor
#    date: 2015-08-10
#    This routine parses the date from a pandas DataFrame when it is divided into several columns
#
#    Parameters:
#    -----------
#    data: pandas DataFrame
#        dataFrame whose dates have to be parsed
#    date_col_names: list
#        A list of the names of the columns in which the date is divided
#        the naming of the date columns must be in accordance with the datetime directives,
#        so if the first column is only the year, its name must be `%Y` and so forth.
#        see https://docs.python.org/2/library/datetime.html#strftime-and-strptime-behavior
#    connector: string
#        should be used only when the default connector causes some conflit
#    first_time_skip: int
#        the offset (mostly because of the bad converting done by LBA
#    clean: bool
#        remove date columns from data after it is introduced as index
#    correct_fracs: bool
#
#    complete_zeroes: list
#        list of columns that need to be padded with zeroes
#    """
#    from datetime import timedelta,datetime
#    import pandas as pd
#
#    #------------------------------------
#    # If dataloggerConfig object is provided, we take the keywords from it
#    if dataloggerConfig:
#        date_col_names = dataloggerConfig.date_col_names
#        connector = dataloggerConfig.date_connector
#        first_time_skip = dataloggerConfig.first_time_skip
#    elif date_col_names==None:
#        raise NameError('Must provide either dataloggerConfig or date_col_names')
#    #------------------------------------
#
#    #------------------------------------
#    # Joins the names of the columns, which must match the datetime directive (see __doc__)
#    if verbose: print('Using these columns: ', date_col_names)
#    date_format=connector.join(date_col_names)
#    auxformat='%Y-%m-%d %H:%M:%S.%f'
#    if complete_zeroes:
#        if type(complete_zeroes) == str:
#            complete_zeroes=[complete_zeroes]
#        for col in complete_zeroes:
#            data[col]=data[col].apply(completeHM)
#    #------------------------------------
#
#    #-------------------------------------
#    # joins the appropriate pandas columns because pandas can read only one column into datetime
#    try:
#        aux=data[date_col_names[0]].astype(str)
#    except ValueError:
#        aux=data[date_col_names[0]].astype(int).astype(str)
#    for col in date_col_names[1:]:
#        aux+=connector + data[col].astype(str)
#    dates=pd.to_datetime(aux, format=date_format)
#    #-------------------------------------
#
#    #-------------------------------------
#    # The next steps are there to check if there are fractions that are not expressed in the datetime convention
#    # and it assumes that the lowest time period expressed is the minute
#    first_date=dates.unique()[1]
#    n_fracs=len(dates[dates.values==first_date])
#    if n_fracs>1:
#        if correct_fracs == None:
#            print('Warning! I identified that there are', n_fracs, ' values (on average) for every timestamp.\n\
#This generally means that the data is sampled at a frequency greater than the frequency of the timestamp. \
#I will then proceed to guess the fractions based of the keyword "first_time_skip" and correct the index.')
#        if (correct_fracs==None) or (correct_fracs==True):
#            dates=[ date.strftime(auxformat) for date in dates ]
#            aux=dates[0]
#            cont=first_time_skip
#            for i,date in enumerate(dates):
#                if date==aux:
#                    pass
#                else:
#                    cont=0
#                    aux=date
#                dates[i]=datetime.strptime(date, auxformat) + timedelta(minutes=cont/float(n_fracs))
#                cont+=1
#        else:
#            print('\nWarning: fractions werent corrected. Check your timestamp data and the correct_fracs flag\n')
#    #-------------------------------------
#
#    #-------------------------------------
#    # setting new dates list as the index
#    data=data.set_index([dates])
#    data.index.name = date_format
#    #-------------------------------------
#
#    #-------------------------------------
#    # removing the columns used to generate the date
#    if clean:
#        data=data.drop(date_col_names, axis=1)
#    #-------------------------------------
#
#    return data
#
#
#def classbin(x, y, bins_number=100, function=np.mean, xfunction=np.mean, logscale=True):
#    '''
#    Separates x and y inputs into bins based on the x array.
#    x and y do not have to be ordered.
#
#    Parameters:
#    -----------
#    x: np.array
#        independent variable
#    y: np.array
#        dependent variable
#    bins_number: int
#        number of classes (or bins) desired
#    function: callable
#        funtion to be applied to both x and y-bins in order to smooth the data
#    logscale: boolean
#        whether or not to use a log-spaced scale to set the bins
#    '''
#    import warnings
#    import numpy as np
#
#    xmin=np.min(x)
#    xmax=np.max(x)
#    if logscale:
#        #-----------
#        # The following if statement gets rid of negative or zero values in the x array, since we are using log-scale
#        if (x<=0).any():
#            y=np.array([ yy for yy, xx in zip(y,x) if xx > 0 ])
#            x=np.array([ el for el in x if el > 0])
#            xmin=np.min(x)
#            xmax=np.max(x)
#        #-----------
#        bins=np.logspace(np.log(xmin), np.log(xmax), bins_number+1, base=np.e)
#    else:
#        bins=np.linspace(xmin, xmax, bins_number+1)
#    xsm = np.zeros(bins_number)
#    ysm = np.zeros(bins_number)
#    #-----------
#    # The following process is what actually bins the data using numpy
#    with warnings.catch_warnings():
#        warnings.simplefilter("ignore", category=RuntimeWarning)
#        for i in range(bins_number):
#            if i == bins_number - 1:
#                sel = bins[i] <= x
#            else:
#                sel = (bins[i] <= x) & (x < bins[i+1])
#            xsm[i] = xfunction(x[sel])
#            ysm[i] = function(y[sel])
#    #-----------
#    return xsm, ysm
#
#
#def get_index(x, y):
#    """
#    Just like the .index method of lists, except it works for multiple values
#
#    Parameter:
#    x: list or array
#        the main array
#    y: list of array
#        the subset of the main whose indexes are desired
#
#    Returns:
#    --------
#    indexes: np.array
#        array with the indexes of each element in y
#    """
#    import numpy as np
#
#    return np.nonzero([ col in y for col in x ])
#

def first_last(fname):
    """
    Returns first and last lines of a file
    """
    with open(fname, 'rb') as fin:
        first=fin.readline()
        for line in fin:
            pass
        last=line
    return first, last

#
#def name2date(filename, dlconfig):
#    """
#    Gets a date from a the name of the file according to a datalogger config object
#
#    Parameters:
#    -----------
#    filename: string
#        the (base) name of the file
#    dlconfig: pymicra.dataloggerConfig
#        configuration of the datalogger
#
#    Returns:
#    --------
#    cdate: datetime object
#
#    Warning:
#    Needs to be optimized in order to read question markers also after the date
#    """
#    from itertools import izip_longest
#    import datetime as dt
#    
#    filename_format=dlconfig.filename_format
#    f=''.join([ s for s,v in izip_longest(filename, filename_format) if v!='?' ])
#    fmt=filename_format.replace('?','')
#    cdate=dt.datetime.strptime(f, fmt)
#    return cdate
#
#
#
#def line2date(line, dlconfig):
#    """
#    Gets a date from a line of file according to dataloggerConfig object.
#
#    Parameters:
#    -----------
#    line: string
#        line of file with date inside
#    dlconfig: pymicra.dataloggerConfig
#        configuration of the datalogger
#
#    Returns:
#    --------
#    timestamp: datetime object
#    """
#    import datetime as dt
#    import numpy as np
#    import re
#
#    varnames=dlconfig.varNames
#    connector = dlconfig.date_connector
#    date_col_names = dlconfig.date_col_names
#    date_cols = dlconfig.date_cols
#    sep = dlconfig.columns_separator
#
#    datefmt=connector.join(date_col_names)
#
#    #-------
#    # Dealing with whitespace as delimiter using regex
#    if sep=='whitespace':
#        sep = r"\s*"
#    #-------
#
#    line = np.array( re.split(sep, line.strip()) )
#    s    = connector.join(line[ date_cols ])
#    return dt.datetime.strptime(s, datefmt)
#
#
#def diff_central(x, y):
#    """
#    Applies the central finite difference scheme
#
#    Parameters:
#    -----------
#    x: array
#        independent variable
#    y: array
#        dependent variable
#
#    Returns:
#    --------
#    dydx: array
#        the dependent variable differentiated
#    """
#    x0 = x[:-2]
#    x1 = x[1:-1]
#    x2 = x[2:]
#    y0 = y[:-2]
#    y1 = y[1:-1]
#    y2 = y[2:]
#    f = (x2 - x1)/(x2 - x0)
#    return (1-f)*(y2 - y1)/(x2 - x1) + f*(y1 - y0)/(x1 - x0)
#

#def find_nearest(array, value):
#    """
#    Smart and small function to find the index of the nearest value, in an array, of some other value
#
#    Parameters:
#    -----------
#    array: array
#        list or array
#    value: float
#        value to look for in the array
#    """
#    import numpy as np
#
#    idx = (np.abs(array-value)).argmin()
#    return idx
#
#
#def mad(data, axis=None):
#    import numpy as np
#    return np.median(np.absolute(data - np.median(data, axis)), axis)
#
#
#def parseUnits(unitstr):
#    '''
#    Gets unit from string, list of strings, or dict's values, using the UnitRegistry
#    defined in __init__.py
#    '''
#    try:
#        from .. import ureg, Q_
#    except ImportError:
#        print('You should have pint installed to use units!')
#        return unitstr
#
#    if isinstance(unitstr, str):
#        return ureg[unitstr]
#    elif isinstance(unitstr, list):
#        return [ ureg[el] for el in unitstr ]
#    elif isinstance(unitstr, dict):
#        return { key: ureg[el] for key, el in unitstr.items() }
#

#def convert_to(data, inunit, outunit, inplace=False, key=None):
#    '''
#    Converts data from one unit to the other
#
#    Parameters:
#    -----------
#    data: pandas.series
#        to be chanhed from one unit to the other
#    inunit: pint.quantity or dict
#        unit(s) that the data is in
#    outunit: str
#        convert to this unit
#    inplace: bool
#        if inunit is a dict, the dict is update in place. "key" keyword must be provided
#    key: str
#        if inunit is a dict, it is the name of the variable to be chamged
#    '''
#    from .. import Q_
#
#    if key:
#        Q = inunit[key].to(outunit)
#    else:
#        Q = inunit.to(outunit)
#
#    coef = Q.magnitude
#    outunit = Q_(1, Q.units)
#    if inplace:
#        inunit.update({key : outunit})
#        return data*coef
#    else:
#        return data*coef, outunit

#
#def get_notation(notation_def):
#    '''
#    Auxiliar function ro retrieve notation
#    '''
#    if notation_def != None:
#        return notation_def
#    else:
#        from .. import notation
#        return notation
#

#def latexify(variables, math_mode=True):
#    '''
#    '''
#    from ..constants import greek_alphabet
#
#    latex = []
#    for var in variables:
#        new = var
#        for letter in greek_alphabet.values():
#            new = new.replace(letter, r'\{}'.format(letter))
#        if math_mode:
#            latex.append('${}$'.format(new))
#        else:
#            latex.append(new)
#    return latex
