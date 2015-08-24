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


def trend(data, mode='moving average', rule='10min', window=None, **kwargs):
    """
    Wrapper to return the moving/block average or polynomial fit of a given data
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
        # performs block average on the data with the window being the rule
        # assumes that frequency is constant
        freq=str((data.index[1]-data.index[0]).microseconds)
        return data.resample(rule, how='mean', **kwargs).resample(freq+'U', fill_method='pad')
    elif any(w in mode for w in ['linear', 'detrend', 'fit']):
        # performs a polynomial fit on the data in blocks of "rule"
        from algs import fitByDate
        return fitByDate(data, rule=rule, **kwargs)
    else:
        # if no mode can be identified
        raise KeyError('Mode defined is not correct. Options are "moving" and "block".')

def detrend(data, variables, mode='moving average', rule='10Min', **kwargs):
    """
    Returns the detrended data of variables
    """
    mode=mode.lower().replace('_',' ').replace('-','').replace(' ','')
    return False 


def spectrum(data, variable=None, frequency=10, absolute=True):
    """
    Author: Tomas Chor

    Calculates the spectrum for a set of data

    Parameters
    ----------

    data: pandas DataFrame/ timeSeries

    col: str
    the column for which to calculate the spectrum

    frequency: float
    frequency of measurement of signal to pass to numpy.fft.rfftfreq

    absolute: bool
    wether or not the results will be given in absolute value
    """
    if variable==None:
        variable=data.columns[0]
    sig=data[variable].values
    spec= np.fft.rfft(sig)
    freq= np.fft.rfftfreq(len(sig), d=1./frequency)
    if absolute==True:
        spec=map(abs,spec)
    aux=pd.DataFrame( data={variable+' spectrum':spec}, index=freq )
    aux.index.name='frequencies'
    return aux



def gradients(data, levels, order='Crescent'):
    """
    Calculates the gradients for data considering the levels provided.
    UNDER DEVELOPMENT
    
    Parameters
    ----------

    data: pandas DataFrame/ timeSeries
    the data for which the gradients should be calculated

    levels: list
    the columns considered to calculate the gradients

    """
    from algs import combine
    flux=pd.DataFrame(index=data.index)
    for pair in combine(levels, order=order):
        a,b=pair
        flux[str(a)+'-'+str(b)]=data[a]-data[b]
    return flux


def obukhovLen(theta_v_star, theta_v_mean, u_star, siteConst):
    """
    Calculates the Monin-Obukhov stability length
    """
    g=siteConst.constants.g
    kappa=siteConst.constants.g
    z=siteConst.variables_height
    d=siteConst.displacement_height
    return - ((z* u_star *u_star* theta_v_mean) / (kappa *g* theta_v_star)



def MonObuVar(theta_v_star, theta_v_mean, u_star, siteConst):
    """
    Calculates the Monin-Obukhov stability variable
    """
    L0=obukhovLen(theta_v_star, theta_v_mean, u_star, siteConst=siteConst)
    return (siteConst.variables_height['u']-siteConst.displacement_height)/L0



def calcLenghts(data,
  mode='linear fit',
  rule='10Min',
  varDict={'u':'u',
  'v':'v',
  'w':'u',
  'pressure':'p',
  'temperature':'T',
  'specific humidity':'q',
  'relative humidity':'rh'},
  **kwargs):
    """
    Calculates the usual characteristic lengths
    """
    from algs import splitData
    chunks=splitData(data, frequency=rule)
    for chunk in chunks:
        u_prime=detrend(data, varDict['u'], mode=mode, **kwargs)
        w_prime=detrend(data, varDict['w'], mode=mode, **kwargs)
        theta_prime=detrend(data, varDict['thermodynamic temp'], mode=mode, **kwargs)
        p_mean=data[varDict['pressure']].mean()
        T_mean=data[varDict['temperature']].mean()
        q_mean=data[varDict['specific humidity']].mean()
        rho_wet=physics.wetAirDens(p=p_mean, T=T_mean, q=q_mean)


#------------------------------------
#
#------------------------------------
class siteConstants(object):
    """
    Author: Tomas Chor

    Keeper of the characteristics and constants of an experiment.

    Attributes:
    -----------
        variables height: should be a dict with the keys being the timeSeries' variables
        canopy height: should be a float
        displacement_height: should be a float (will be estimated as 2/3*canopy_height if not given)
        variables_units: a dict with the units for every variable to which units are aplicable
        gravity: the acceleration of gravity in m/(s*s) (assumed to be 9.81 if not provided)
        R: Universal gas constant
        Rs: Specific gas constant for dry air
        Rv: Specific gas constant for water vapor
        mu: Rv/Rs
    """
    from physics import constants

    def __init__(self, variables_height, canopy_height,
             displacement_height=None,
             variables_units=None,
             gravity=9.81,
             Rs=287.04, 
             Rv=461.50, 
             cp=1003.5):
        # Checks if the variables_height is actually a Dict
        if not isinstance(variables_height, dict):
            raise TypeError('variables_height should be a dictionary. Ex.: {"u" : 10, "v" : 10, "theta" : 12 }')
        self.variables_height = variables_height #meters
        self.canopy_height = canopy_height         #meters
        if displacement_height==None:
            self.displacement_height = (2./3.)*self.canopy_height #meters
        else:
            self.displacement_height=displacement_height
        self.gravity = gravity        #meters/(s**2)
        self.Rs = Rs    #specific gas constant for dry air J/(kg.K)
        self.Rv = Rv    #J/(kg.K)
        self.mu = Rv/Rs
        self.cp = cp     #specific heat for constant pressure J/(kg.K)


