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
import notation

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


def get_scales(data, siteConst, notation_defs=None,
  output_as_df=True, solutes=[]):
    """
    Calculates characteristic lengths for data

    The names of the variables are retrived out the dictionary. You can update the dictionary
    and change the names by using the notation_defs keyworkd, which is a notation object

    Parameters:
    -----------
    data: pandas.DataFrame
        dataset to be used. Must have a minimum of columns in order for it to work
    siteConst: pymicra.siteConstants object
        currently not used
    updt: dictionary
        use this dictionary to change a small amount of the names that the variables have in the dataframe
    output_as_df: boolean
        True if you want the output to be a one-line pandas.DataFrame

    CHECKLIST:
    ADD MIXED-LAYER CONVECTION SCALES FOR VELOCITY (w*) AND TEMPERATURE (t*)
    """
    if notation_defs==None:
        defs=notation.get_notation()
    else:
        defs=notation_defs
    flup, flus = defs.fluctuation_preffix, defs.fluctuation_suffix
    u       =   flup + defs.u + flus
    v       =   flup + defs.v + flus
    w       =   flup + defs.w + flus
    p       =   defs.pressure
    theta   =   defs.thermodyn_temp
    theta_v =   defs.virtual_temp
    theta_v_fluc= flup + defs.virtual_temp + flus
    q           = defs.specific_humidity
    qfluct      = flup + defs.specific_humidity + flus
    solutesf    = [ flup + el + flus for el in solutes ]
    
    u_star=np.sqrt(-algs.auxCov( data[[u,w]] ))
    u_std=data[u].std()

    theta_v_star=algs.auxCov( data[[theta_v_fluc,w]] )/u_star
    theta_v_mean=data[theta_v].mean()
    theta_v_std=data[theta_v].std()

    q_star=algs.auxCov( data[[qfluct,w]] ) / u_star
    q_mean=data[q].mean()
    q_std=data[qfluct].std()

    c_stars=[]
    c_stds =[]
    for c in solutesf:
        c_stars.append( algs.auxCov( data[[c,w]] ) / u_star )
        c_stds.append( data[c].std() )

    theta_mean=data[theta].mean()
    theta_std=data[theta].std()
    theta_star=(theta_v_star - 0.61*theta_mean*q_star)/(1.+0.61*q_mean)

    Lm=MonObuLen(theta_v_star, theta_v_mean, u_star, g=siteConst.constants.gravity, kappa=siteConst.constants.kappa)
    zeta=MonObuSimVar(Lm, siteConst)

    if output_as_df:
        namespace=locals()
        columns=['zeta', 'Lm', 'u_std', 'u_star', 'theta_v_std', 'theta_v_star', 'theta_std', 'theta_star', 'q_std', 'q_star']
        dic={ col : [namespace[col]] for col in columns }
        out=pd.DataFrame(dic, index=[data.index[0]])
        # We have to input the solutes separately
        for solute, c_star, c_std in zip(solutes, c_stars, c_stds):
                out[ '{}_std'.format(solute) ] = c_std
                out[ '{}_star'.format(solute)] = c_star
        return out
    else:
        return zeta, Lm, (u_std, u_star), (theta_v_std, theta_v_star), (theta_std, theta_star), (q_std, q_star), (c_std, c_star)




def ste(data, w_fluctuations="w'"):
    """
    Returns the Symmetric Transfer Efficiency in the time domain, ste
    according to Cancelli, Dias, Chamecki, Dimensionless criteria for the production-dissipation equilibrium
    of scalar fluctuations and their implications for scalar similarity, Water Resources Research, 2012

    NEEDS TO BE VALIDATED!
    """
    from data import bulkCorr
    wcol=w_fluctuations
    df=data.copy()
    w=df[wcol]
    df=df.drop(wcol, axis=1)
    a,b=df.columns
    a,b=df[a],df[b]
    rwa=bulkCorr([w,a])
    rwb=bulkCorr([w,b])
    rwa=np.abs(rwa)
    rwb=np.abs(rwb)
    return 1. - np.abs( rwa - rwb )/( rwa + rwb )
 

def eddyCov(data, wpl=True,
        notation_defs=None, solutes=[]):
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
    import constants
    cp=constants.cp_dry

    if notation_defs==None:
        defs=notation.get_notation()
    else:
        defs=notation_defs

    stap, stas = defs.star_preffix, defs.star_suffix

    p       =   defs.pressure
    theta_mean   =   defs.mean_preffix + defs.thermodyn_temp + defs.mean_suffix
    theta_v =   defs.virtual_temp
    theta_star= stap + defs.thermodyn_temp + stas
    theta_v_star= stap + defs.virtual_temp + stas
    u_star  =   stap + defs.u + stas
    q_star  = stap + defs.specific_humidity + stas

    mu=1./constants.mu
    rho_mean=data['rho_air_mean']
    rho_h2o_mean=data['rho_h2o_mean']
    rho_co2_mean=data['rho_co2_mean']
    rho_dry_mean=data['rho_dry_mean']
    c_stars = []
    for solute in solutes:
        c_stars.append( stap + solute + stas )

    out=pd.DataFrame(index=data.index)
    out['tau']=rho_mean* (u_star**2.)
    out['H']=  rho_mean* cp* u_star* theta_star
    out['Hv']= rho_mean* cp* u_star* theta_v_star
    out['E']=  rho_mean* u_star* q_star
    out['F']=  rho_mean* u_star* c_star
    #------------------------
    # APPLY WPL CORRECTION. PAGES 34-35 OF MICRABORDA
    if wpl:
        rv=rho_h2o_mean/rho_dry_mean
        rc=rho_co2_mean/rho_dry_mean
        out['E']=   (1. +mu*rv)*( out['E'] + rho_h2o_mean*( (theta_star*u_star)/theta_mean))
        out['F']=   out['F'] + rho_co2_mean*(1. + mu*rv)*(theta_star*u_star)/theta_mean + mu*rc*out['E'] 
    #------------------------
    return out


def phi(zeta, x=None):
    """
    Currently using Businger-Dyer eqs.
    """
    if zeta<0:
        if x=='tau':
            return (1. - 16.*zeta)**(1./4.)
        else:
            return np.sqrt(1. - 16.*zeta)
    if zeta>0:
        return 1. + 5.*zeta




def _Cx(x, za, zb, d, z0, Lm, Psif=None):
    """
    Get Cx for gradient-flux method.
    In this code, zb means tau and za means anything else and the bottom level
    is the ground (height zero)

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
    #--------
    # Checks for Psi and if it's a function
    if Psif==None:
        print 'Should get a repository for functions'
    elif hasattr(Psif, '__call__')==False:
        raise TypeError('Psi argument has to be a function')
    else:
        Psif = Psi
    #--------

    if x=='tau':
        cx=1./(auxLogMinPsi(zb, d, z0, Psif, Lm))**2.
    elif x=='H':
        pass
    elif x=='E':
        pass
    elif x=='F':
        cx = 1./( auxLogMinPsi(zb, d, z0, Psif, Lm) * auxLogMinPsi(za, d, z0, Psif, Lm) )
    else:
        raise NameError('x is not on the list. See help(Cx) for more information')
    cx*=kappa**2
    return cx



def _auxLogMinPsi(z_i, d, z0, Psi, Lm):
    """
    Auxiliar function that is used in the denominator of the C_tau-like coefficients

    Taken from Dias, Micrometeorological Approaches to to estimate fluxes of greenhouse gases between the surface and the
    atmosphere, (Chapter 7, eqs. 94-97, p. 22.)

    CONSIDER TAKING THIS SOMEWHERE ELSE!
    """
    if hasattr(Psi, '__call__')==False:
        raise TypeError('Psi argument has to be a function')
    a=np.log((z_i-d)/z0)
    b=Psi((z_i-d)/Lm)
    return a-b



def Psi(zeta, x='tau', zeta0=0.):
    """
    Integral Monin-Obukhov scale or deviation function, which is the deviation of a
    variable (x) in relation nto their logarithmic profiles due the stability zeta != 0

    Parameters:
    -----------
    zeta: float
        the stability variable
    x: string
        the variable. Options are 'tau', 'H', 'E', 'F'
    zeta0: float
        value of zeta_{0 x}. Only used for the stable case
    """
    #---------
    # For unstable conditions
    if zeta < 0:
        b = (1. - 16.*zeta)**(1./4.)
        if x=='tau':
            psi = np.log((b*b + 1.)/2.) + 2.*np.log((b+1.)/2.) - 2.*np.arctan(b) + np.pi/2.
        if any([ x == var for var in ['H', 'E', 'F'] ]):
            psi = 2.*np.log((b*b + 1.)/2.)
        return psi
    #---------

    #---------
    # For stable conditions
    else:
        return -5.*(zeta - zeta0)
    #---------





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
             variables_units=None):
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


