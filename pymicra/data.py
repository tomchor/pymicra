#!/usr/bin/python
"""
Author: Tomas Chor
Date: 2015-08-07
-------------------------

This module works with micrometeorological data using pandas, numpy, datetime and several other packages

-------------------------

Modifications:


"""
import pandas as pd
import numpy as np
import physics

def rotCoor(data, wind_vars=['u','v','w'], verbose=0):
    """
    Rotates the coordinates of wind data
    ----------------------

    Parameters
    ----------
    data: pandas DataFrame 
        the dataFrame to be rotated
    wind_vars: list
        a list with the names of columns which define the wind speed, typically ['u','v','w']
    """
    from numpy import cos,sin,zeros,dot
    from math import atan2, sqrt
    DC= np.zeros((3,3))
    wind_vec = data[wind_vars].mean().values
    m_u, m_v, m_w = wind_vec
    alpha = atan2(m_v,m_u)
    beta =-atan2(m_w, sqrt((m_u**2.)+(m_v**2.)))
    if verbose:
        print "alpha: ",alpha,"\nbeta: ",beta
    # definition of rotation matrix
    DC[0,0],DC[0,1],DC[0,2] = np.cos(alpha)*np.cos(beta), np.cos(beta)*np.sin(alpha),-np.sin(beta)
    DC[1,0],DC[1,1],DC[1,2] =-np.sin(alpha)          , np.cos(alpha)          , 0.
    DC[2,0],DC[2,1],DC[2,2] = np.cos(alpha)*np.sin(beta), np.sin(alpha)*np.sin(beta), np.cos(beta)
    if verbose:
        print "Rotation matrix is:\n", DC
    # application of rotation as a matrix product
    data[wind_vars] = np.dot(DC, data[wind_vars].values.T).T
    return data


def trend(data, mode='moving average', rule=None, window=None, **kwargs):
    """
    Wrapper to return the trend given data. Can be achieved using a moving avg, block avg or polynomial fitting
    -------------

    Parameters
    ----------

    data: pandas DataFrame 
        the dataFrame to be rotated
    mode: string
        mode of average to apply. Currently {'moving', 'block'}.
    rule: string
        pandas offset string to define the block in the block average. Default is "10min".
    """
    mode=mode.lower().replace('_',' ').replace('-','').replace(' ','')
    if any(w==mode for w in ['moving', 'movingaverage']):
        # performs moving average on the data with window equivalent to the rule
        if window==None:
            print('warning: approximating window size for moving average')
            window=int(len(data)/len(data.resample(rule)))
        return pd.rolling_mean(data, window=window, **kwargs)
    elif any(w==mode for w in ['block', 'blockaverage']):
        # performs block average on the data with the window being the rule. Assumes that frequency is constant
        if rule==None:
            return data.apply(lambda x: [np.mean(x)]*len(x), axis=0)
        else:
            print 'Warning. Might be bugged. Check results.'
            freq=data.index.inferred_freq
            return data.resample(rule, how='mean', **kwargs).resample(freq, fill_method='pad')
    elif any(w in mode for w in ['linear', 'polynomial', 'fit']):
        # performs a polynomial fit on the data in blocks of "rule"
        from algs import fitByDate
        return fitByDate(data, rule=rule, **kwargs)
    else:
        # if no mode can be identified
        raise KeyError('Mode defined is not correct. Options are "moving" and "block".')

def detrend(data, mode='moving average', rule=None, suffix="'", **kwargs):
    """
    Returns the detrended fluctuations of a given dataset

    Parameters
    ----------

    data: pandas.DataFrame, pandas.Series

    mode: string
        what method to use in order to identify the trend
    rule: pandas offset string
        the blocks for which the trends should be calculated
    suffix: string
        suffix to add to variable names after fluctuation is extracted
    """
    from algs import stripDown
    mode=stripDown(mode.lower(), args='-_')
    data=data-trend(data, mode=mode, rule=rule, **kwargs)
    return data.add_suffix(suffix)


def spectrum(data, window='1min', frequency=10, absolute=True, T=30):
    """
    Calculates the spectrum for a set of data

    Parameters
    ----------

    data: pandas DataFrame/ timeSeries

    frequency: float
    frequency of measurement of signal to pass to numpy.fft.rfftfreq

    absolute: bool
    wether or not the results will be given in absolute value

    T: int, float
    period in minutes
    """
    T=T*60.
    sig2=None
    var1=''
    if type(data)==pd.Series:
        sig=data.values
    elif type(data)==pd.DataFrame:
        if len(data.columns)==1:
            var1=data.columns[0]
        elif len(data.columns)==2:
            var1,var2=data.columns
            sig2=data[var2].values
        else:
            raise Exception('Too many columns of data. Chose one (spectrum) or two (cross-spectrum')
        sig=data[var1].values
    else:
        if len(data)==1:
            sig=data
        elif len(data)==2:
            sig,sig2=data
        else:
            raise Exception('Too many columns of data. Chose one (spectrum) or two (cross-spectrum')
    spec= np.fft.rfft(sig)
    freq= np.fft.rfftfreq(len(sig), d=1./frequency)
    if sig2 != None:
        spec2=np.fft.rfft(sig2)
        spec=np.real((2./T)*(spec*spec2.conjugate()))
    elif absolute==True:
        spec=np.real((2./T)*(spec*spec.conjugate()))
    aux=pd.DataFrame( data={var1+' spectrum':spec}, index=freq )
    aux.index.name='frequencies'
    return aux


def spectrum_old(data, variable=None, frequency=10, absolute=True, T=30):
    """
    Author: Tomas Chor

    Calculates the spectrum for a set of data

    Parameters
    ----------
    SHOULD NOT BE USED - PROBABLY OBSOLETE

    data: pandas DataFrame/ timeSeries

    col: str
        the column for which to calculate the spectrum
    frequency: float
        frequency of measurement of signal to pass to numpy.fft.rfftfreq
    absolute: bool
        whether or not the results will be given in absolute value
    T: int, float
        period in minutes
    """
    T=T*60.
    if variable==None:
        variable=data.columns[0]
    sig=data[variable].values
    spec= np.fft.rfft(sig)
    freq= np.fft.rfftfreq(len(sig), d=1./frequency)
    if absolute==True:
        #specnp.real((2./T)*(spec*spec.conjugate()))
        spec=map(abs,spec)
    aux=pd.DataFrame( data={variable+' spectrum':spec}, index=freq )
    aux.index.name='frequencies'
    return aux



