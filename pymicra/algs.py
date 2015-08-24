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
    res_index=data.resample(frequency).index
    out=[]
    pdate=res_index[0]
    for date in res_index:
        aux=data[pdate:date][:-1]
        if len(aux.values)>0:
            out.append(aux)
        pdate=date
    out.append(data[pdate:])
    return out


def fitByDate(data, degree=1, rule=None):
    """
    Given a pandas DataFrame with the index as datetime, this routine
    fit a n-degree polynomial to the dataset

    SHOULD LATER ADD MULTIPLE COLUMNS FUNCTIONALITY
    """
    if rule==None:
        dflist=[data]
    else:
        dflist=splitData(data, frequency=rule)
    out=pd.DataFrame()
    for data in dflist:
        y=data.values
        x=data.index.to_julian_date()
        coefs=np.polyfit(x,y,degree)
        yy=np.polyval(coefs,x)
        aux=pd.DataFrame(yy, index=data.index)
        out=out.append(aux)
    return out

