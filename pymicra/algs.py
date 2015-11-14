#!/usr/bin/python
"""

"""
import pandas as pd
import numpy as np


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


def splitData(data, frequency='30Min'):
    """
    Splits a given pandas DataFrame into a series of "frequency"-spaced DataFrames

    Parameters
    ----------

    data: pandas dataframe

    frequency: pandas string offset
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
    res_index=data[data.columns[0]].resample(frequency).index
    out=[]
    pdate=res_index[0]
    for date in res_index:
        aux=data[pdate:date][:-1]
        if len(aux.values)>0:
            out.append(aux)
        pdate=date
    out.append(data[pdate:])
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


def testValid(df_valid, testname='', falseverbose=True, trueverbose=True):
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
        if falseverbose:
            failed=df_valid[df_valid==False].index
            print 'Failed',testname,'test'
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
        Z=np.sqrt(2)*erfinv(-2.*phi+1.)
        return sigma*Z + mu
    return f

