#!/usr/bin/python
"""

"""
import pandas as pd
import numpy as np
import scipy.stats as st


def combine(levels, order='Crescent'):
    """ Combines the elements of the levels in non-repeting pairs
    for a n length array, gives n*(n-1)/2 combinations

    Parameters
    ----------
    levels: list
    list whose elements should be combined
    """
    print order
    if order.lower()=='decrescent':
        levels=levels[::-1]
    combinations=[]
    for i in range(0, len(levels)):
        for j in range(i+1, len(levels)):
            if j>i:
                combinations.append([levels[j],levels[i]])
    return combinations


def splitData(data, rule='30min', return_index=False, **kwargs):
    """
    Splits a given pandas DataFrame into a series of "rule"-spaced DataFrames

    Parameters
    ----------
    data: pandas dataframe
        data to be split
    rule: pandas string offset
        Some possible values (that should be followed by an integer) are:
        D   calendar day frequency
        W   weekly frequency
        M   month end frequency
        MS  month start frequency
        Q   quarter end frequency
        BQ  business quarter endfrequency
        QS  quarter start frequency
        A   year end frequency
        AS  year start frequency
        H   hourly frequency
        T   minutely frequency
        Min minutely frequency
        S   secondly frequency
        L   milliseconds
        U   microseconds
        
        check it complete at http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases
        """
    res_index = pd.Series(index=data.index).resample(rule, **kwargs).index
    out=[]
    pdate=res_index[0]
    for date in res_index:
        aux=data.ix[pdate:date]#[:-1]
        print aux
        #aux=data[pdate:date][:-1]
        if len(aux.values)>0:
            out.append(aux)
        pdate=date
    out.append(data[pdate:])
    if return_index:
        return out, res_index
    else:
        return out

def fitWrap(x,y,degree=1):
    """
    A wrapper to numpy.polyfit and numpy.polyval that fits data given an x and y arrays.
    This is specifically designed to be used with by pandas.DataFrame.apply method

    Parameters
    ----------
    x: array, list
    y: array, list
    degree: int
    """
    coefs=np.polyfit(x,y,degree)
    yy=np.polyval(coefs,x)
    return yy

def fitByDate(data, degree=1, rule=None):
    """
    Given a pandas DataFrame with the index as datetime, this routine
    fit a n-degree polynomial to the dataset

    Parameters
    ----------
    data: pd.DataFrame, pd.Series
        dataframe whose columns have to be fitted
    """
    if rule==None:
        dflist=[data]
    else:
        dflist=splitData(data, frequency=rule)
    out=pd.DataFrame()
    for data in dflist:
        xx=data.index.to_julian_date()
        aux=data.apply(lambda x: fitWrap(xx, x, degree=degree), axis=0)
        out=out.append(aux)
    return out


def stripDown(str, final='', args=['_', '-']):
    """
    Auxiliar function to strip down kaywords from symbols
    """
    for arg in args:
        str=str.replace(arg,final)
    return str
        

def auxCov(data):
    """
    Auxiliar function to obtain covariances in a pandas.DataFrame

    data: pandas.DataFrame, pandas.Series
        the data whose covariance you want to obtain. Maximum number of columns is 2
    """
    colnum=len(data.columns)
    if colnum==1:
        return np.mean( data[data.columns[0]]*data[data.columns[0]] )
    elif colnum==2:
        return np.mean( data[data.columns[0]]*data[data.columns[1]] )
    else:
        raise TypeError('Datset with more than two columns')



def lenYear(year):
    """
    Calculates the length of a year in days
    Useful to figure out if a certain year is a leap year
    """
    import calendar
    feblen=calendar.monthrange(year,2)[1]
    otherlens=365-28
    return feblen + otherlens


def limited_interpolation(data, maxcount=3):
    '''
    Interpolates linearly but only if gap is smaller of equal to maxcout

    Parameters:
    -----------

    data: pandas.DataFrame
        dataset to interpolate
    maxcount: int
        maximum number of consecutive NaNs to interpolate. If the number is smaller than that, nothing is done with the points.
    '''
    mask = data.copy()
    grp = ((mask.notnull() != mask.shift().notnull()).cumsum())
    grp['ones'] = 1
    for i in data.columns:
        mask[i] = (grp.groupby(i)['ones'].transform('count') <= maxcount) | data[i].notnull()
    data=data.interpolate(mode='time').bfill()[mask]
    return data


def inverse_normal_cdf(mu, sigma):
    """
    Inverse normal CDF

    Parameters
    """
    from scipy.special import erfinv
    def f(phi):
        Z=math.sqrt(2)*erfinv(-2.*phi+1.)
        return sigma*Z + mu
    return f


def completeHM(string):
    if string.isdigit():
        pass
    else:
        raise TypeError('String passed must contain only digits. Check the argument')
    if len(string)==3:
        string='0'+string
    return string


def limitedSubs(data, max_interp=3, func=lambda x: abs(x) > abs(x.std()*4.) ):
    """
    """
    df=data.copy()
    cond=func(df)
    for c in df.columns:
        grouper = (cond[c] != cond[c].shift(1)).cumsum() * cond[c]
        fill = (df.groupby(grouper)[c].transform('size') <= max_interp)
        df.loc[fill, c] = np.nan
    return df


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
    '''
    
    if False in df_valid.values:
        failed=df_valid[df_valid==False].index
        if falseverbose:
            print filepath, 'failed',testname,'test'
            print 'Failed variable(s):', ', '.join(failed)
            print
        return False, failed
    else:
        if trueverbose: print filepath,'Passed',testname,'test'
        return True, None


def file_len(fname):
    """
    Returns length of a file through piping bash's function wc

    Parameters:
    -----------

    fname: string
        path of the file
    """
    import subprocess
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, 
                                              stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])

def check_numlines(fname, numlines=18000):
    """
    Checks length of file against a correct value.
    Returns False is length is wrong and True if length is right

    Parameters:
    ----------

    fname: string
        path of the file to check
    numlines: int
        correct number of lines that the file has to have
    """
    lines=file_len(fname)
    if lines==numlines:
        return True
    else:
        return False

def inverse_normal_cdf(mu, sigma):
    """
    """
    from scipy.special import erfinv
    def f(phi):
        Z=np.sqrt(2.)*erfinv(-2.*phi+1.)
        return sigma*Z + mu
    return f

def applyResult(result, failed, df, control=None, testname=None, filename=None, falseshow=False):
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
        if control:
            control[testname].append(filename)
    return control


def parseDates(data, date_cols, connector='-', first_time_skip=0,
  clean=True, correct_fracs=None, complete_zeroes=False):
    """
    Author: Tomas Chor
    date: 2015-08-10
    This routine parses the date from a pandas DataFrame when it is divided into several columns

    Parameters:
    -----------
    data: pandas DataFrame
        dataFrame whose dates have to be parsed
    date_cols: list of strings
        A list of the names of the columns in which the date is divided
        the naming of the date columns must be in accordance with the datetime directives,
        so if the first column is only the year, its name must be `%Y` and so forth.
        see https://docs.python.org/2/library/datetime.html#strftime-and-strptime-behavior
    connector: string
    first_time_skip: int
        the offset (mostly because of the bad converting done by LBA
    clean: bool
        remove date columns from data after it is introduced as index
    correct_fracs: bool

    complete_zeroes: list
        list of columns that need to be padded with zeroes
    """
    from datetime import timedelta,datetime
    from algs import completeHM
    #------------------------------------
    # joins the names of the columns, which must match the datetime directive (see __doc__)
    #------------------------------------
    date_format=connector.join(date_cols)
    auxformat='%Y-%m-%d %H:%M:%S.%f'
    if complete_zeroes:
        if type(complete_zeroes) == str:
            complete_zeroes=[complete_zeroes]
        for col in complete_zeroes:
            data[col]=data[col].apply(completeHM)
    #-------------------------------------
    # joins the appropriate pandas columns because pandas can read only one column into datetime
    #-------------------------------------
    try:
        aux=data[date_cols[0]].astype(str)
    except ValueError:
        aux=data[date_cols[0]].astype(int).astype(str)
    for col in date_cols[1:]:
        aux+=connector + data[col].astype(str)
    dates=pd.to_datetime(aux, format=date_format)
    #-------------------------------------
    # The next steps are there to check if there are fractions that are not expressed in the datetime convention
    # and it assumes that the lowest time period expressed is the minute
    #-------------------------------------
    first_date=dates.unique()[1]
    n_fracs=len(dates[dates.values==first_date])
    if n_fracs>1:
        if correct_fracs == None:
            print 'Warning! I identified that there are', n_fracs, ' values (on average) for every timestamp.\n\
This generally means that the data is sampled at a frequency greater than the frequency of the timestamp. \
I will then proceed to guess the fractions based of the keyword "first_time_skip" and correct the index.'
        if (correct_fracs==None) or (correct_fracs==True):
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
        else:
            print '\nWarning: fractions werent corrected. Check your timestamp data and the correct_fracs flag\n'
    #-------------------------------------
    # setting new dates list as the index
    #-------------------------------------
    data=data.set_index([dates])
    #-------------------------------------
    # removing the columns used to generate the date
    #-------------------------------------
    if clean:
        data=data.drop(date_cols, axis=1)
    return data


def classbin(x, y, bins_number=100, function=np.mean, log_scale=True):
    '''
    Separates x and y inputs into bins based on the x array.
    x and y do not have to be ordered.

    Parameters:
    -----------
    x: np.array
        independent variable
    y: np.array
        dependent variable
    bins_number: int
        number of classes (or bins) desired
    function: callable
        funtion to be applied to both x and y-bins in order to smooth the data
    log_scale: boolean
        whether or not to use a log-spaced scale to set the bins
    '''
    import warnings
    xmin=np.min(x)
    xmax=np.max(x)
    if log_scale:
        #-----------
        # The following if statement gets rid of negative or zero values in the x array, since we are using log-scale
        #-----------
        if (x<=0).any():
            print 'Warning: zero and/or negative values exist in x array about to be log-scaled. Will try to ignore but errors might arise.'
            y=np.array([ yy for yy, xx in zip(y,x) if xx > 0 ])
            x=np.array([ el for el in x if el > 0])
            xmin=np.min(x)
            xmax=np.max(x)
        bins=np.logspace(np.log(xmin), np.log(xmax), bins_number+1, base=np.e)
    else:
        bins=np.linspace(xmin, xmax, bins_number+1)
    xsm = np.zeros(bins_number)
    ysm = np.zeros(bins_number)
    #-----------
    # The following process probably should be optimized
    #-----------
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        for i in range(bins_number):
            if i == bins_number - 1:
                sel = bins[i] <= x
            else:
                sel = (bins[i] <= x) & (x < bins[i+1])
            xsm[i] = function(x[sel])
            ysm[i] = function(y[sel])
    return xsm, ysm


def binwrapper(self, **kwargs):
    x=self.index.astype(np.float64)
    out=pd.DataFrame(columns = self.columns)
    for c in self.columns:
        xsm, ysm = classbin(self[c].astype(np.float64))
        out[c] = ysm
    out.index=xsm
    return out
pd.DataFrame.binned=binwrapper


def classlogbin(maxcl, indx, x, y, pr_sign=+1.0, geometric_mean=True, function=None):
    '''
    Author: Nelson L. Dias
    Modified by: Tomas L. Chor (2015-11-23) to include arithmetic mean

    Class-averages a big array: the abscissas are divided into maxcl
    classes, then it returns the mean of the abscissas and the ordinates for
    each class. The classes are log-spaced so as to appear linear in logplot
    This function contemplates the possibility that some classes may be empty.

    Parameters:
    -----------
    maxcl: int
        the maximum # of classes desired
    indx: list, array
        the index ordering x from smallest to biggest
    x: list, array
        the big array of abscissas
    y: list, array
        the big array of ordinates
    pr_sign:
        prevailing sign 
    geometric_mean: bool
        whether or not to use geometric mean. If False, arithmetic mean is used
    function: function
        overwrites argument geometric_mean and uses this function to perform the binning
    NEEDS TO BE TESTED!

    Returns:
    --------
    nsm: int
        the number of non-empty classes
    npclsign:
        the number of points used for averaging per class
    xsm: list
        the smoothed array of abscissas
    ysm: list
        the smoothed array of ordinates
    '''
    if geometric_mean:
        cmean=st.gmean
    else:
        cmean=np.mean
    if function != None:
        assert hasattr(function, '__call__')
        cmean=function
    ntotal = len(x)
# ------------------------------------------------------------------------------
# number of points per class: will need two of them!
# ------------------------------------------------------------------------------
    npclsign = np.zeros(maxcl,int)
    npclass = np.zeros(maxcl,int)
# ------------------------------------------------------------------------------
# the smoothed arrays
# ------------------------------------------------------------------------------
    xsm = np.zeros(maxcl,float)
    ysm = np.zeros(maxcl,float)
# ------------------------------------------------------------------------------
# log increments
# ------------------------------------------------------------------------------
    logxmin = np.log(x[indx[0]])
    logxmax = np.log(x[indx[ntotal-1]])
# ------------------------------------------------------------------------------
# counts the number of points per class
# ------------------------------------------------------------------------------
    i = 0
    classe = 0
    while ( i < ntotal ):
        npclass[classe] = 0
# ------------------------------------------------------------------------------
# warning! truncation errors?
# ------------------------------------------------------------------------------
        alfa = float(classe+1)/float(maxcl)
        assert ( alfa <= 1.0 )
        logxsup = alfa * logxmax + (1.0 - alfa)*logxmin
        while ( (i < ntotal) and ( np.log(x[indx[i]]) <= logxsup) ):
            npclass[classe] += 1
            i += 1
        classe += 1
    assert ( classe == maxcl )
# ------------------------------------------------------------------------------
# begins averaging
# ------------------------------------------------------------------------------
    nsm = 0 
    first = 0 
    for classe in range(maxcl): 
        if ( npclass[classe] > 0 ) :
# ------------------------------------------------------------------------------
# needs a local array for storing contiguously the members of a class    
# ------------------------------------------------------------------------------
            xlocal = np.zeros( npclass[classe],float) 
            ylocal = np.zeros( npclass[classe],float)
# ------------------------------------------------------------------------------
# contiguous storage
# ------------------------------------------------------------------------------
            k = 0 
            for j in range(npclass[classe]):
                indlocal = indx[first+j]
# ------------------------------------------------------------------------------
# checks index consistency
# ------------------------------------------------------------------------------
                if ( indlocal > ntotal - 1 ):
                    print("stat-->classlogavg: incorrect index in class averaging\n")
                    exit(1)
# ------------------------------------------------------------------------------
# checks signal consistency
# ------------------------------------------------------------------------------
                yy = pr_sign * y[indlocal] 
                if ( yy > 0.0 ) :
                    xlocal[k] = x[indlocal] 
                    assert( xlocal[k] <= 10.0)
                    # INCLUDED .REAL HERE IN ORDER TO AVOID WARNING
                    ylocal[k] = yy.real
                    k += 1
# ------------------------------------------------------------------------------
# the VALID number of points per class is now limited to those points having the
# prevailing sign: calculates the x- and y- averages of each class for non-empty
# classes
# ------------------------------------------------------------------------------
            if ( k > 0 ) :
                xavg = cmean(xlocal[0:k])
                yavg = cmean(ylocal[0:k])
                xsm[nsm] = xavg 
                ysm[nsm] = pr_sign * yavg 
                npclsign[nsm] = k 

                nsm += 1
            pass
        pass

        first += npclass[classe] 
    return (nsm, npclsign, xsm, ysm)


def get_index(x, y):
    return np.nonzero([ col in y for col in x ])

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


def line2date(line, dlconfig):
    import datetime as dt
    varnames=dlconfig.varNames
    date_cols = dlconfig.date_cols
    sep = dlconfig.columns_separator

    datefmt=' '.join(date_cols)
    indexes = get_index(varnames, date_cols)
    line=np.array(line.split(sep))
    s=' '.join(line[indexes])
    return dt.datetime.strptime(s, datefmt)

def diff_central(x, y):
    x0 = x[:-2]
    x1 = x[1:-1]
    x2 = x[2:]
    y0 = y[:-2]
    y1 = y[1:-1]
    y2 = y[2:]
    f = (x2 - x1)/(x2 - x0)
    return (1-f)*(y2 - y1)/(x2 - x1) + f*(y1 - y0)/(x1 - x0)
