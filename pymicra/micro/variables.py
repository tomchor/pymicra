#!/usr/bin/python
"""
Author: Tomas Chor
Date: 2015-08-07
-------------------------

This module works with micrometeorological data using pandas, numpy, datetime and several other packages

-------------------------

TO-DO LIST
add direct EddyCov method (without stars)

Maybe add all these variables that EddyPro calculates:
http://www.licor.com/env/help/EddyPro3/Content/Topics/Calculating_Micromet_Variables.htm

"""
import pandas as pd
import numpy as np
from .. import notation
from .. import constants

def MonObuVar(L_m, siteConst):
    """
    Redirects to MonObuSimVar
    """
    return MonObuSimVar(L_m, siteConst)

def MonObuSimVar(L_m, siteConst):
    """
    Calculates the Monin-Obukhov Similarity Variable
    defined as

    zeta = (z-d)/Lm
    where d is the displacement height or zero-plane displacement
    and L_m is the Monin-Obukhov Length.
    """
    z=siteConst.variables_height
    d=siteConst.displacement_height
    zeta = (z-d)/L_m
    return zeta


def MonObuLen(theta_v_star, theta_v_mean, u_star, g=None):
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
    if g==None:
        g = constants.gravity
    kappa = constants.kappa
    Lm = - ( (u_star**2) * theta_v_mean) / (kappa *g* theta_v_star)
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
    ADD MIXED-LAYER CONVECTION SCALES FOR VELOCITY (w*) AND TEMPERATURE (t*) MAYBE
    """
    if notation_defs==None:
        defs=notation.get_notation()
    else:
        defs=notation_defs

    #---------
    # First we define the names of the columns according to notation
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
    #---------

    #-----------
    # Then we create the covariances matrix
    cov = data[[u,w,theta_v_fluc, qfluct] + solutesf ].cov()
    #-----------
    
    #---------
    # Now to calculate the characteristic lengths, scales and etc
    u_star  = np.sqrt(-cov.loc[u,w])
    u_std   = data[u].std()

    theta_v_star    = cov.loc[theta_v_fluc, w] / u_star
    theta_v_mean    = data[theta_v].mean()
    theta_v_std     = data[theta_v].std()

    q_star  = cov.loc[qfluct,w] / u_star
    q_mean  = data[q].mean()
    q_std   = data[qfluct].std()

    c_stars =[]
    c_stds  =[]
      #---------
      # The solutes have to be calculated separately
    for c in solutesf:
        c_stars.append( cov.loc[c, w] / u_star )
        c_stds.append( data[c].std() )
      #---------
    theta_mean=data[theta].mean()
    theta_std=data[theta].std()
    theta_star=(theta_v_star - 0.61*theta_mean*q_star)/(1.+0.61*q_mean)
    #---------

    #---------
    # Now we calculate the obukhov length and the similarity variable
    Lm=MonObuLen(theta_v_star, theta_v_mean, u_star, g=constants.gravity)
    zeta=MonObuSimVar(Lm, siteConst)
    #---------

    #---------
    # Finally we construct the output dataframe
    if output_as_df:
        namespace=locals()
        columns=['zeta', 'Lm', 'u_std', 'u_star', 'theta_v_std', 'theta_v_star', 'theta_std', 'theta_star', 'q_std', 'q_star']
        dic={ col : [namespace[col]] for col in columns }
        out=pd.DataFrame(dic, index=[data.index[0]])
        #-----------
        # We have to input the solutes separately
        for solute, c_star, c_std in zip(solutes, c_stars, c_stds):
                out[ '{}_std'.format(solute) ] = c_std
                out[ '{}_star'.format(solute)] = c_star
        #-----------
        return out
    else:
        return zeta, Lm, (u_std, u_star), (theta_v_std, theta_v_star), (theta_std, theta_star), (q_std, q_star), (c_std, c_star)
    #---------



def eddyCov(data, wpl=True,
        notation_defs=None, solutes=[]):
    """
    Get fluxes according to characteristic lengths
    
    WARNING! If using wpl=True, be sure that all masses are consistent!
        For example, if q = [g/g], rho_h2o = [g/m3] and rho_co2 = [g/m3] and so on.
        Avoid mixing kg/m3 with g/m3 (e.g. for co2 and h2o) and mg/kg with g/g (e.g. for
        co2 and h2o).

    TODO: add support for pint units

    Parameters:
    -----------
    data: pandas.DataFrame
        dataframe with the characteristic lengths calculated
    notation_defs: pymicra.notation
        object that holds the notation used in the dataframe
    wpl: boolean
        whether or not to apply WPL correction on the latent heat flux and solutes flux
    solutes: list
        list that holds every solute considered for flux
    """
    cp=constants.cp_dry

    #---------
    # Define useful notation to look for
    if notation_defs==None:
        defs=notation.get_notation()
    else:
        defs=notation_defs
    stap, stas = defs.star_preffix, defs.star_suffix
    #---------

    #---------
        # Define name of variables to look for based on the notation
    p           =   defs.pressure
    theta_mean  =   defs.mean_preffix + defs.thermodyn_temp + defs.mean_suffix
    theta_v     =   defs.virtual_temp
    theta_star  =   stap + defs.thermodyn_temp + stas
    theta_v_star=   stap + defs.virtual_temp + stas
    u_star      =   stap + defs.u + stas
    q_star      =   stap + defs.specific_humidity + stas
    #---------

    #---------
    # Define auxiliar variables
    rho_mean=data['rho_air_mean']
    rho_h2o_mean=data['rho_h2o_mean']
    rho_co2_mean=data['rho_co2_mean']
    rho_dry_mean=data['rho_dry_mean']
    c_stars = []
    for solute in solutes:
        c_stars.append( stap + solute + stas )
    #---------

    #---------
    # Calculate the fluxes
    out=pd.DataFrame(index=data.index)
    out['tau']=rho_mean* ( data[u_star]**2.)
    out['H'] = rho_mean* cp* data[u_star]* data[theta_star]
    out['Hv']= rho_mean* cp* data[u_star]* data[theta_v_star]
    out['E'] = rho_mean* data[u_star]* data[q_star]
    for solute, c_star in zip(solutes, c_stars):
        out[ 'F_{}'.format(solute) ] =  rho_mean* data[u_star]* data[c_star]
    #---------

    #------------------------
    # APPLY WPL CORRECTION. PAGES 34-35 OF MICRABORDA
    if wpl:
        mu=constants.R_spec['h2o']/constants.R_spec['dry']
        rv=rho_h2o_mean/rho_dry_mean
        rc=rho_co2_mean/rho_dry_mean
        out['E'] = (1. +mu*rv)*( out['E'] + rho_h2o_mean * (data[theta_star]*data[u_star]/data[theta_mean]) )
        for solute, c_star in zip(solutes, c_stars):
            out[ 'F_{}'.format(solute) ] = \
                out[ 'F_{}'.format(solute) ] + \
                rho_co2_mean*(1. + mu*rv)*(data[theta_star]*data[u_star])/data[theta_mean] + \
                mu*rc*out['E']
    #------------------------
    return out



def eddyCov2(data, wpl=True,
        notation2_defs=None, solutes=[]):
    """
    Get fluxes according to characteristic lengths
    
    WARNING! If using wpl=True, be sure that all masses are consistent!
        For example, if q = [g/g], rho_h2o = [g/m3] and rho_co2 = [g/m3] and so on.
        Avoid mixing kg/m3 with g/m3 (e.g. for co2 and h2o) and mg/kg with g/g (e.g. for
        co2 and h2o).

    TODO: add support for pint units

    Parameters:
    -----------
    data: pandas.DataFrame
        dataframe with the characteristic lengths calculated
    notation_defs: pymicra.notation
        object that holds the notation used in the dataframe
    wpl: boolean
        whether or not to apply WPL correction on the latent heat flux and solutes flux
    solutes: list
        list that holds every solute considered for flux
    """
    cp=constants.cp_dry

    #---------
    # Define useful notation to look for
    if notation2_defs==None:
        defs=notation.get_notation2()
    else:
        defs=notation_defs2
    fluct = defs.fluctuation
    mean = defs.mean
    #---------

    #---------
    # Define name of variables to look for based on the notation
    u               =   fluct % defs.u
    w               =   fluct % defs.w
    p               =   defs.pressure
    q_fluc          =   fluct % defs.specific_humidity
    theta_mean      =   mean % defs.thermodyn_temp
    theta_v_fluc    =   fluct % defs.virtual_temp
    solutesf        = [ fluct % solute for solute in solutes ]
    #---------

    #---------
    # Now we try to calculate or identify the fluctuations of theta
    try:
        theta_fluc  =   fluct % defs.thermodyn_temp
    except:
        #---------
        # We need the mean of the specific humidity
        try:
            data_q_mean  =   data[ mean % defs.specific_humidity ]
        except:
            data_q_mean  =   data[ defs.specific_humidity ].mean()
        #---------

        data_theta_fluc  =   ( data[theta_v_fluc] - 0.61*data[theta_mean]*data[q_fluc])/(1.+0.61*data_q_mean)
    #---------

    #---------
    # First we construct the covariance matrix
    cov = data[[u,w,theta_v_fluc, q_fluc] + solutesf ].cov()
    #---------

    #---------
    # Define auxiliar variables
    mu=constants.R_spec['h2o']/constants.R_spec['dry']
    rho_air_mean    =   defs.mean % defs.density % defs.moist_air
    rho_h2o_mean    =   defs.mean % defs.density % defs.h2o
    rho_co2_mean    =   defs.mean % defs.density % defs.co2
    rho_dry_mean    =   defs.mean % defs.density % defs.dry_air
    #---------

    #---------
    # Calculate the fluxes
    out=pd.DataFrame(index=data.index)
    out['tau']  = rho_air_mean* cov[u, w]
    out['H']    = rho_air_mean* cp* cov[theta_fluc, w]
    out['Hv']   = rho_air_mean* cp* cov[theta_v_fluc, w]
    out['E']    = rho_air_mean* cov[q_fluc, w]
    for solute, solutef in zip(solutes, solutesf):
        out[ 'F_{}'.format(solute) ] =  rho_air_mean* cov[solutef, w]
    #---------

    #------------------------
    # APPLY WPL CORRECTION. PAGES 34-35 OF MICRABORDA
    if wpl:
        rv=rho_h2o_mean/rho_dry_mean
        rc=rho_co2_mean/rho_dry_mean
        out['E'] = (1. +mu*rv)*( out['E'] + rho_h2o_mean * (data[theta_star]*data[u_star]/data[theta_mean]) )
        for solute, c_star in zip(solutes, c_stars):
            out[ 'F_{}'.format(solute) ] = \
                out[ 'F_{}'.format(solute) ] + \
                rho_co2_mean*(1. + mu*rv)*(data[theta_star]*data[u_star])/data[theta_mean] + \
                mu*rc*out['E']
    #------------------------
    return out





#------------------------------------
# CLASSES
#------------------------------------
class siteConstants(object):
    """
    Keeper of the characteristics and constants of an experiment.

    Attributes:
    -----------
        variables_height: float
            the main height of the instruments in meters. Generally the height of the sonic anemometer
        canopy height: float
            the mean height of the vegetation meters
        displacement_height: float
            also called zero-plane displacement. Will be estimated as 2/3*canopy_height if not given.
    """
    def __init__(self, variables_height, canopy_height,
             displacement_height=None, z0=None, description=None):

        self.description=description
        self.variables_height = variables_height    #meters
        self.canopy_height = canopy_height          #meters
        self.z0 = z0
        self.description = description
        if displacement_height==None:
            self.displacement_height = (2./3.)*self.canopy_height #meters
        else:
            self.displacement_height=displacement_height


