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

def MonObuSimVar(L_m, siteConst):
    """
    Calculates the Monin-Obukhov Similarity Variable
    according to:

    GARRAT, 
    zeta = z/L
    """
    z=siteConst.variables_height
    d=siteConst.displacement_height
    return (z-d)/L_m


def MonObuLen(theta_v_star, theta_v_mean, u_star, g=9.81, kappa=0.4):
    """
    Calculates the Monin-Obukhov Length
    according to:

    GARRAT, The atmospheric boundary layer, 1992 (eq. 1.11, p. 10)
    L = ( u_star^2 * theta_v ) / ( kappa * g * theta_v_star )

    KUNDU, Fluid mechanics, 1990 (eq 71, chap. 12, p. 462)
    L_M = - u_star^3 / ( kappa * alpha * g * cov(w,T') )

    ARYA, Introduction to micrometeorology (eq. 11.1, p. 214)
    L = - u_star^3 / (kappa * (g/T_0) * (H_0/(rho*c_p)) )

    STULL, An introduction to Boundary layer meteorology, 1988 (eq. 5.7b, p. 181)
    L = - ( theta_v * u_star^3 ) / ( kappa *g* cov(w',theta_v') )
    """
    return - ( (u_star**2) * theta_v_mean) / (kappa *g* theta_v_star)


def get_scales(data, siteConst,  
  varDict={'u':r"u'",
  'v':r"v'",
  'w':r"w'",
  'pressure':'p',
  'temperature':'theta_v',
  'temperature fluctuations':r"theta_v'",
  'specific humidity':r"q'",
  'relative humidity':'rh'},
  updt={}):

    """
    Calculates characteristics lengths for data
    Assumes
    u_*^2 = -mean(u' * w')
    theta_* = -mean(theta' * w') / u_*
    theta_v_* = -mean(theta_v' * w') / u_*
    q_* = -mean(q' * w') / u_*

    The names of the variables are retrived out the dictionary. You can update the dictionary
    and change the names by using the updt keyword with the keys
    u,v,w,pressure, temperature, temperature fluctuations, specific humidity, relative humidity

    CHECKLIST:
    NEEDS IMPROVEMENT IN ORDER TO GET CONSTANTS FROM SITECONST OBJECT
    ADD MIXED-LAYER CONVECTION SCALES FOR VELOCITY (w*) AND TEMPERATURE (t*)
    """
    varDict.update(updt)
    u=varDict['u']
    v=varDict['v']
    w=varDict['w']
    p=varDict['pressure']
    theta_v=varDict['temperature']
    theta_v_fluc=varDict['temperature fluctuations']
    q=varDict['specific humidity']
    
    u_star=np.sqrt(-algs.auxCov( data[[u,w]] ))
    theta_v_star=algs.auxCov( data[[theta_v_fluc,w]] )/u_star
    theta_v_mean=data[theta_v].mean()
    q_star=algs.auxCov( data[[q,w]] )/u_star

    L_m=MonObuLen(theta_v_star, theta_v_mean, u_star, g=siteConst.constants.gravity, kappa=siteConst.constants.kappa)
    zeta=MonObuSimVar(L_m, siteConst)
    return zeta, u_star, (L_m, theta_v_star, q_star)


def get_fluxes(u_star, q_star, theta_star, theta_v_star, c_star, rho_mean, cp=1):
    """
    Get fluxes according to char lengths
    
    TO-DO LIST:
    add more concentrations to the variables
    """
    tau=rho_mean* (u_star**2.)
    H=  rho_mean* cp* u_star* theta_star
    E=  rho_mean* cp* u_star* q_star
    Hv= rho_mean* cp* u_star* theta_v_star
    F=  rho_mean* u_star* c_star
    return, tau, H, E, Hv, F


def phi_H(zeta):
    """
    Currently using Businger-Dyer eqs.

    TO-DO LIST:
    include more types
    """
    if zeta<0:
        return np.sqrt(1. - 16.*zeta)
    if zeta>0:
        return 1. + 5.*zeta



#
#def calcLenghts(data,
#  mode='linear fit',
#  rule='10Min',
#  varDict={'u':'u',
#  'v':'v',
#  'w':'u',
#  'pressure':'p',
#  'temperature':'T',
#  'specific humidity':'q',
#  'relative humidity':'rh'},
#  **kwargs):
#    """
#    Calculates the usual characteristic lengths
#    """
#    from algs import splitData
#    chunks=splitData(data, frequency=rule)
#    for chunk in chunks:
#        u_prime=detrend(data, varDict['u'], mode=mode, **kwargs)
#        w_prime=detrend(data, varDict['w'], mode=mode, **kwargs)
#        theta_prime=detrend(data, varDict['thermodynamic temp'], mode=mode, **kwargs)
#        p_mean=data[varDict['pressure']].mean()
#        T_mean=data[varDict['temperature']].mean()
#        q_mean=data[varDict['specific humidity']].mean()
#        rho_wet=physics.wetAirDens(p=p_mean, T=T_mean, q=q_mean)
#




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

#        if not isinstance(variables_height, dict):
#            raise TypeError('variables_height should be a dictionary. Ex.: {"u" : 10, "v" : 10, "theta" : 12 }')
        import constants
        self.constants=constants
        self.variables_height = variables_height #meters
        self.canopy_height = canopy_height         #meters
        if displacement_height==None:
            self.displacement_height = (2./3.)*self.canopy_height #meters
        else:
            self.displacement_height=displacement_height
        self.cp = cp     #specific heat for constant pressure J/(kg.K)


