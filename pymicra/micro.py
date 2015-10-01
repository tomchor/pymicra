#!/usr/bin/python
"""
Author: Tomas Chor
Date: 2015-08-07
-------------------------

This module works with micrometeorological data using pandas, numpy, datetime and several other packages

-------------------------

TO-DO LIST

Maybe add all these variables that EddyPro calculates:
http://www.licor.com/env/help/EddyPro3/Content/Topics/Calculating_Micromet_Variables.htm

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
    Lm= - ( (u_star**2) * theta_v_mean) / (kappa *g* theta_v_star)
    return Lm


def get_scales(data, siteConst,
  varDict={
  'u'   :r"u'",
  'v'   :r"v'",
  'w'   :r"w'",
  'pressure'        :'p',
  'temperature'     :'theta',
  'virtual temperature' :'theta_v',
  'virtual temperature fluctuations'    :"theta_v'",
  'co2 fluctuations':    "co2'",
  'specific humidity'   :r"q'",
  'relative humidity'   :'rh'},
  updt={},
  output_as_df=True):


    """
    Calculates characteristics lengths for data

    Assumes:
    u_*^2 = -mean(u' * w')
    theta_* = -mean(theta' * w') / u_*
    theta_v_* = -mean(theta_v' * w') / u_*
    q_* = -mean(q' * w') / u_*

    The names of the variables are retrived out the dictionary. You can update the dictionary
    and change the names by using the updt keyword with the keys
    u,v,w,pressure, temperature, temperature fluctuations, specific humidity, relative humidity

    Parameters:
    -----------
    data: pandas.DataFrame
        dataset to be used. Must have a minimum of columns in order for it to work
    siteConst: pymicra.siteConstants object
        currently not used
    varDict: dictionary
        dictionary saying which names the variables have in the dataframe
    updt: dictionary
        use this dictionary to change a small amount of the names that the variables have in the dataframe
    output_as_df: boolean
        True if you want the output to be a one-line pandas.DataFrame

    CHECKLIST:
    NEEDS IMPROVEMENT IN ORDER TO GET CONSTANTS FROM SITECONST OBJECT
    ADD MIXED-LAYER CONVECTION SCALES FOR VELOCITY (w*) AND TEMPERATURE (t*)
    """
    varDict.update(updt)
    u=varDict['u']
    v=varDict['v']
    w=varDict['w']
    p=varDict['pressure']
    theta=varDict['temperature']
    theta_v=varDict['virtual temperature']
    theta_v_fluc=varDict['virtual temperature fluctuations']
    q=varDict['specific humidity']
    c=varDict["co2 fluctuations"]
    
    u_star=np.sqrt(-algs.auxCov( data[[u,w]] ))
    u_mean=data[u].mean()

    theta_v_star=algs.auxCov( data[[theta_v_fluc,w]] )/u_star
    theta_v_mean=data[theta_v].mean()

    q_star=algs.auxCov( data[[q,w]] )/u_star
    q_mean=data[q].mean()

    c_star=algs.auxCov( data[[c,w]] )/u_star
    c_mean=data[c].mean()

    theta_mean=data[theta].mean()
    theta_star=(theta_v_star - 0.61*theta_mean*q_star)/(1.+0.61*q_mean)

    Lm=MonObuLen(theta_v_star, theta_v_mean, u_star, g=siteConst.constants.gravity, kappa=siteConst.constants.kappa)
    zeta=MonObuSimVar(Lm, siteConst)
    if output_as_df:
        namespace=locals()
        columns=['zeta', 'Lm', 'u_mean', 'u_star', 'theta_v_mean', 'theta_v_star', 'theta_mean', 'theta_star', 'q_mean', 'q_star', 'c_mean', 'c_star']
        dic={ col : [namespace[col]] for col in columns }
        out=pd.DataFrame(dic, index=[data.index[0]])
        return out
    else:
        return zeta, Lm, (u_mean, u_star), (theta_v_mean, theta_v_star), (theta_mean, theta_star), (q_mean, q_star), (c_mean, c_star)


def get_fluxes_DF(data, cp=None, wpl=True):
    """
    Get fluxes according to char lengths
    
    add more concentrations to the variables

    Parameters:
    -----------
    data: pandas.DataFrame
        dataframe with the characteristic lengths calculated
    cp: float (optional)
        value for the specific heat capacity at constant pressure
    """
    if cp==None:
        from constants import cp_dry
    cp=cp_dry
    rho_mean=data['rho_mean']
    u_star=data['u_star']
    theta_star=data['theta_star']
    theta_v_star=data['theta_v_star']
    q_star=data['q_star']
    c_star=data['c_star']
    out=pd.DataFrame(index=data.index)
    out['tau']=rho_mean* (u_star**2.)
    out['H']=  rho_mean* cp* u_star* theta_star
    out['Hv']= rho_mean* cp* u_star* theta_v_star
    # APPLY WPL CORRECTION. PAGES 34-35 OF MICRABORDA
    if wpl:
        out['E']=  rho_mean* u_star* q_star
        out['F']=  rho_mean* u_star* c_star
    else:
        out['E']=   (1+mu*rv)*(cov(w,rho_v) + rho_v_mean*(cov(w*theta)/theta_mean))
        out['F']=   rho_c*(1 + mu*rv)*cov(w,theta)/theta_mean + mu*rc*cov(w,rho_v) + cov(w,rho_c)
    return out


def get_fluxes(u_star, q_star, theta_star, theta_v_star, c_star, rho_mean, cp=None):
    """
    Get fluxes according to char lengths
    
    TO-DO LIST:
    add more concentrations to the variables
    """
    if cp==None:
        from constants import cp_dry
    cp=cp_dry
    tau=rho_mean* (u_star**2.)
    H=  rho_mean* cp* u_star* theta_star
    Hv= rho_mean* cp* u_star* theta_v_star
    # APPLY WPL CORRECTION. PAGES 34-35 OF MICRABORDA
    E=  rho_mean* u_star* q_star
    F=  rho_mean* u_star* c_star
    return tau, H, E, Hv, F


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




def Cx(x, za, zb, d, z0, Lm, Psif=None):
    """
    Get Cx for gradient-flux method.

    Parameters:
    -----------

    x: string
        {'tau', 'H', 'E', 'c', 'F'}
    za: float
        height at which substance (water, c2o, temperature, etc) is measured
    zb: float
        height at which wind velocity is measured
    d: float
        displacement height
    z0: float
        roughness length
    Lm: float
        Monin-Obukhov Length
    Psif: function (optional)
        deviation function Psi=log(abs(zeta)) - Phi(zeta)
    """
    from constants import kappa
    if Psif==None:
        print 'Should get a repository for functions'
    elif hasattr(Psif, '__call__')==False:
        raise TypeError('Psi argument has to be a function')
    else:
        pass
    if x=='tau':
        cx=1./(auxLogMinPsi(zb,d,z0,Psif,Lm))**2.
    elif x=='H':
        pass
    elif x=='E':
        pass
    elif x=='F':
        pass
    else:
        raise NameError('x is not on the list. See help(Cx) for more information')
    cx*=kappa**2
    return cx



def auxLogMinPsi(za, d, z0, Psi, Lm):
    """
    Auxiliar function that is used in the denominator of the C_tau-like coefficients

    Taken from Dias, Micrometeorological Approaches to to estimate fluxes of greenhouse gases between the surface and the
    atmosphere, (Chapter 7, eqs. 94-97, p. 22.)

    CONSIDER TAKING THIS SOMEWHERE ELSE!
    """
    if hasattr(Psi, '__call__')==False:
        raise TypeError('Psi argument has to be a function')
    a=np.log((za-d)/z0)
    b=Psi((za-d)/Lm)
    return a-b

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


