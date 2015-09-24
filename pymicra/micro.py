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
import algs

def obukhovLen(theta_v_star, theta_v_mean, u_star, siteConst):
    """
    Calculates the Monin-Obukhov stability length
    """
    g=siteConst.constants.gravity
    kappa=siteConst.constants.kappa
    z=siteConst.variables_height['u']
    d=siteConst.displacement_height
    zeta=MonObuVar(theta_v_star, theta_v_mean, u_star, g=g)
    return (z-d)/zeta


def MonObuVar(theta_v_star, theta_v_mean, u_star, g=9.81):
    """
    Calculates the Monin-Obukhov stability variable
    """
    kappa=.4
    return - (kappa *g* theta_v_star) / (u_star *u_star* theta_v_mean)


def get_scales(data, u="u'", w="w'", theta="theta'", theta_v="theta_v'", q="q'"):
    """
    Assumes
    u_*^2 = -mean(u' * w')
    theta_* = -mean(theta' * w') / u_*
    theta_v_* = -mean(theta_v' * w') / u_*
    q_* = -mean(q' * w') / u_*
    """
    return np.sqrt(-algs.auxCov(data[uw]))



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
    def __init__(self, variables_height, canopy_height,
             displacement_height=None,
             variables_units=None,
             cp=1003.5):
        # Checks if the variables_height is actually a Dict
        if not isinstance(variables_height, dict):
            raise TypeError('variables_height should be a dictionary. Ex.: {"u" : 10, "v" : 10, "theta" : 12 }')
        self.constants=physics.constants()
        self.variables_height = variables_height #meters
        self.canopy_height = canopy_height         #meters
        if displacement_height==None:
            self.displacement_height = (2./3.)*self.canopy_height #meters
        else:
            self.displacement_height=displacement_height
        self.cp = cp     #specific heat for constant pressure J/(kg.K)


