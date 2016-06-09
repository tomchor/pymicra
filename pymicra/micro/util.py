#!/usr/bin/python
from __future__ import print_function
"""
Maybe add

molecular_weight of moist air
water vapor partial pressure
water vapor partial pressure at saturation
relative humidity?
dew point temperature?
partial pressure of dry air
dry air molar volume?
dry air heat capacity at constant pressure?
water vapor heat capacity at constant pressure?
refining ambient temperature?
moist air heat capacity at constant pressure?
specific evaporation heat?
* water to dry air mixing ratio
"""

def preProcess(data, units, notation_defs=None,
        rho_air_from_theta_v=True, inplace=True, solutes=[]):
    '''
    Pre-processes data by calculating moist and dry air densities, specific humidity
    mass density and other important variables

    Parameters:
    -----------
    data: pandas.DataFrame
        dataframe with micrometeorological measurements
    units: dict
        units dictionary with the columns of data as keys
    notation_defs: pymicra.notation
        defining notation used in data
    rho_air_from_theta_v: bool
        whether to use theta_v to calculate air density or theta
    '''
    from .. import constants
    from .. import notation
    from .. import algs
    from .. import physics
    from .. import ureg

    data = data.copy()
    if not inplace:
        units = units.copy()

    Rh2o = constants.R_spec['h2o']
    Mh2o = constants.molar_mass['h2o']
    molar_mass_unit = constants.units['molar_mass']

    #---------
    # Define useful notation to look for
    if notation_defs==None:
        defs=notation()
    else:
        defs=notation_defs
    #---------

    #---------
    # First convert the temperature if it is still in Celsius
    if units[ defs.thermodyn_temp ] == ureg[ 'degC' ]:
        print('Converting theta to kelvin ... ', end='')
        data.loc[:, defs.thermodyn_temp ] = data[ defs.thermodyn_temp ].apply(physics.CtoK)
        units.update({ defs.thermodyn_temp : ureg['kelvin'] })
        print('Done!')

    if units[ defs.virtual_temp ] == ureg[ 'degC' ]:
        print('Converting theta_v to kelvin ... ', end='')
        data.loc[:, defs.virtual_temp ] = data[ defs.virtual_temp ].apply(physics.CtoK)
        units.update({ defs.virtual_temp : ureg['kelvin'] })
        print('Done!')
    #---------

    #---------
    # Check for h2o mass density
    if defs.h2o_density not in data.columns:
        print("Didn't locate mass density of h2o. Trying to calculate it ... ", end='')
        data.loc[:, defs.h2o_density ] = data.loc[:, defs.h2o_molar_density ]*Mh2o
        units.update({ defs.h2o_density : units[ defs.h2o_molar_density ]*molar_mass_unit })
        print("Done!")
    data.loc[:, defs.h2o_density] = algs.convert_to(data[ defs.h2o_density ], units, 'kg/m**3', inplace=True, key=defs.h2o_density)
    #---------

    #---------
    # Check for h2o molar density
    if defs.h2o_molar_density not in data.columns:
        print("Didn't locate molar density of h2o. Trying to calculate it ... ", end='')
        data.loc[:, defs.h2o_molar_density ] = data.loc[:, defs.h2o_density ]/Mh2o
        units.update({ defs.h2o_molar_density : units[ defs.h2o_density ]/molar_mass_unit })
        print("Done!")
    #---------

    #---------
    # Calculation of rho_air is done here
    if (defs.moist_air_density not in data.columns):
        if rho_air_from_theta_v:
            print('Calculating rho_air = p/(Rdry * theta_v) ... ', end='')
            data = physics.airDensity_from_theta_v(data, units, notation=defs, inplace=True)
            print('Done!')
        else:
            print('calculation of air density from theta has to be implemented')
            pass
            #algs.airDensity_from_theta()
    #---------

    #---------
    # Calculation of rho_dry is done here
    if (defs.dry_air_density not in data.columns):
        print('Calculating rho_dry_air = rho_air - rho_h2o ... ', end='')
        data.loc[:, defs.dry_air_density ] = data[ defs.moist_air_density ] - data[ defs.h2o_density ]
        units.update({ defs.dry_air_density : units[ defs.moist_air_density ] })
        print('Done!')
    #---------
 
    #---------
    # Calculation of specific humidity is done here
    if (defs.specific_humidity not in data.columns):
        print('Calculating specific humidty = rho_h2o / rho_air ... ', end='')
        data.loc[:, defs.specific_humidity] = data[ defs.h2o_density ] / data[ defs.moist_air_density ]
        units.update({ defs.specific_humidity : units[ defs.h2o_density ] / units[ defs.moist_air_density ] })
        print('Done!')
        

    #---------
    # Calculation of mixing ratio is done here
    if (defs.h2o_mixing_ratio not in data.columns):
        print('Calculating h2o mixing ratio = rho_h2o / rho_dry ... ', end='')
        data.loc[:, defs.h2o_mixing_ratio] = data[ defs.h2o_density ] / data[ defs.dry_air_density ]
        units.update({ defs.h2o_mixing_ratio : units[ defs.h2o_density ] / units[ defs.dry_air_density ] })
        print('Done!')
    #---------

    #---------
    # Calculation of h2o molar concentration is done here
    #if (defs.h2o_molar_concentration not in data.columns):
    #    print('Calculating h2o molar concentration = mrho_h2o / mrho_air ... ', end='')
    #    data.loc[:, defs.h2o_molar_concentration] = data[ defs.h2o_molar_density ] / data[ defs.moist_air_molar_density ]
    #    units.update({ defs.h2o_molar_concentration : units[ defs.h2o_molar_density ] / units[ defs.moist_air_molar_density ] })
    #    print('Done!')
    #---------

    #---------
    # Here we deal with the SOLUTES!
    for solute in solutes:
        sol_density = eval('defs.{}_density'.format(solute))
        sol_molar_density = eval('defs.{}_molar_density'.format(solute))
        sol_concentration = eval('defs.{}_concentration'.format(solute))
        sol_mixing_ratio = eval('defs.{}_mixing_ratio'.format(solute))
        M_sol = constants.molar_mass[solute]

        #---------
        # Check for solute mass density
        if sol_density not in data.columns:
            print("Didn't locate mass density of {}. Trying to calculate it ... ".format(solute), end='')
            data.loc[:, sol_density ] = data.loc[:, sol_molar_density ]*M_sol
            units.update({ sol_density : units[ sol_molar_density ]*molar_mass_unit })
            print("Done!")
        data.loc[:, sol_density] = algs.convert_to(data[ sol_density ], units, 'kg/m**3', inplace=True, key=sol_density)
        #---------
    
        #---------
        # Check for solute molar density
        if sol_molar_density not in data.columns:
            print("Didn't locate molar density of {}. Trying to calculate it".format(solute))
            data.loc[:, sol_molar_density ] = data.loc[:, sol_density ]/M_sol
            units.update({ sol_molar_density : units[ sol_density ]/molar_mass_unit })
            print("Done!")
        #---------

        #---------
        # Calculation of SOLUTE concentration (g/g) is done here
        if (sol_concentration not in data.columns):
            print('Calculating {0} concentration (g/g) = rho_{0} / rho_air ... '.format(solute), end='')
            data.loc[:, sol_concentration] = data[ sol_density ] / data[ defs.moist_air_density ]
            units.update({ sol_concentration : units[ sol_density ] / units[ defs.moist_air_density ] })
            print("Done!")
            
    
        #---------
        # Calculation of SOLUTE mixing ratio is done here
        if (sol_mixing_ratio not in data.columns):
            print('Calculating {0} mixing ratio = rho_{0} / rho_dry ... '.format(solute), end='')
            data.loc[:, sol_mixing_ratio] = data[ sol_density ] / data[ defs.dry_air_density ]
            units.update({ sol_mixing_ratio : units[ sol_density ] / units[ defs.dry_air_density ] })
            print("Done!")
        #---------
    #---------
        
    
    if inplace:
        return data
    else:
        return data, units


def eddyCov(data, wpl=True,
        notation_defs=None, solutes=[], from_fluctuations=True, from_scales=False):
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
    from .. import constants
    from ..core import notation
    import pandas as pd

    if from_fluctuations:
        return eddyCov3(data, wpl=wpl, notation_defs=notation_defs, solutes=solutes)

    cp=constants.cp_dry
    lamb = constants.latent_heat_water

    #---------
    # Define useful notation to look for
    if notation_defs==None:
        defs=notation()
    else:
        defs=notation_defs
    star = defs.star
    #---------

    #---------
    # Define name of variables to look for based on the notation
    p           =   defs.pressure
    theta_mean  =   defs.mean % defs.thermodyn_temp
    theta_v     =   defs.virtual_temp
    theta_star  =   star % defs.thermodyn_temp
    theta_v_star=   star % defs.virtual_temp
    u_star      =   star % defs.u
    q_star      =   star % defs.specific_humidity
    #---------

    #---------
    # Define auxiliar variables
    rho_mean=data['rho_air_mean']
    rho_h2o_mean=data['rho_h2o_mean']
    rho_co2_mean=data['rho_co2_mean']
    rho_dry_mean=data['rho_dry_mean']
    c_stars = []
    for solute in solutes:
        c_stars.append( star % solute )
    #---------

    #---------
    # Calculate the fluxes
    out=pd.DataFrame(index=data.index)
    out[ defs.momentum_flux ]               = rho_mean* ( data[u_star]**2.)
    out[ defs.sensible_heat_flux ]          = rho_mean* cp* data[u_star] * data[theta_star]
    out[ defs.virtual_sensible_heat_flux ]  = rho_mean* cp* data[u_star] * data[theta_v_star]
    out[ defs.water_vapor_flux]             = rho_mean* data[u_star] * data[q_star]
    out[ defs.latent_heat_flux ]            = lamb( data[theta_mean] ) * rho_mean * data[u_star] * data[q_star]
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
        notation_defs=None, solutes=[]):
    """
    Get fluxes from the turbulent fluctuations
    
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
    from .. import constants
    from ..core import notation
    import pandas as pd

    cp = constants.cp_dry
    lamb = constants.latent_heat_water

    #---------
    # Define useful notation to look for
    if notation_defs==None:
        defs=notation()
    else:
        defs=notation_defs
    fluct = defs.fluctuation
    mean = defs.mean
    #---------

    #---------
    # Define name of variables to look for based on the notation
    u               =   fluct % defs.u
    w               =   fluct % defs.w
    q_fluc          =   fluct % defs.specific_humidity
    theta_fluc      =   fluct % defs.thermodyn_temp
    theta_v_fluc    =   fluct % defs.virtual_temp
    solutesf        = [ fluct % solute for solute in solutes ]
    #---------

    #---------
    # Now we try to calculate or identify the fluctuations of theta
    data_theta_mean = data[ defs.thermodyn_temp ].mean()
    try:
        data_theta_fluc  =   data[ fluct % defs.thermodyn_temp ]
    except:
        #---------
        # We need the mean of the specific humidity and temperature
        data_q_mean  =   data[ defs.specific_humidity ].mean()
        #---------

        data_theta_fluc  =   ( data[theta_v_fluc] - 0.61*data_theta_mean*data[q_fluc])/(1.+0.61*data_q_mean)
        data[ theta_fluc ] = data_theta_fluc
    #---------

    #---------
    # First we construct the covariance matrix
    cov = data[[u,w,theta_v_fluc, theta_fluc, q_fluc] + solutesf ].cov()
    #---------

    #---------
    # Define auxiliar variables
    mu=constants.R_spec['h2o']/constants.R_spec['dry']
    rho_air_mean    =   data[ defs.density % defs.moist_air ].mean()
    rho_h2o_mean    =   data[ defs.density % defs.h2o ].mean()
    rho_co2_mean    =   data[ defs.density % defs.co2 ].mean()
    rho_dry_mean    =   data[ defs.density % defs.dry_air ].mean()
    #---------

    #---------
    # Calculate the fluxes
    out=pd.DataFrame(index=[data.index[0]])
    out[ defs.momentum_flux ]               = -rho_air_mean * cov[u][w]
    out[ defs.sensible_heat_flux ]          = rho_air_mean * cp * cov[theta_fluc][w]
    out[ defs.virtual_sensible_heat_flux ]  = rho_air_mean * cp * cov[theta_v_fluc][w]
    out[ defs.water_vapor_flux ]            = rho_air_mean * cov[q_fluc][w]
    out[ defs.latent_heat_flux ]            = lamb( data_theta_mean ) * rho_air_mean * cov[q_fluc][w]
        #---------
        # Calculate flux for each solute
    for solute, solutef in zip(solutes, solutesf):
        out[ defs.flux_of % solute ] =  rho_air_mean * cov[solutef][w]
        #---------
    #---------

    #------------------------
    # APPLY WPL CORRECTION. PAGES 34-35 OF MICRABORDA
    if wpl:
        rv=rho_h2o_mean/rho_dry_mean
        rc=rho_co2_mean/rho_dry_mean
        out[ defs.water_vapor_flux ] = (1. +mu*rv)*( out[ defs.water_vapor_flux ] + rho_h2o_mean * (cov[theta_fluc][w]/data_theta_mean) )
        for solute in solutes:
            out[ defs.flux_of % solute ] = \
                out[ defs.flux_of % solute ] + \
                rho_co2_mean*(1. + mu*rv)*(cov[theta_fluc][w])/data_theta_mean + \
                mu*rc*out[ defs.water_vapor_flux ]
    #------------------------
    return out



def eddyCov3(data, units, wpl=True,
        notation=None, inplace=True, solutes=[]):
    """
    Get fluxes from the turbulent fluctuations
    
    WARNING! If using wpl=True, be sure that all masses are consistent!
        For example, if q = [g/g], rho_h2o = [g/m3] and rho_co2 = [g/m3] and so on.
        Avoid mixing kg/m3 with g/m3 (e.g. for co2 and h2o) and mg/kg with g/g (e.g. for
        co2 and h2o).

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
    from .. import constants
    from .. import algs
    import pandas as pd

    cp = constants.cp_dry
    lamb = constants.latent_heat_water

    defs = algs.get_notation(notation)
    data = data.copy()
    units = units.copy()
    cunits = constants.units

    #---------
    # Define name of variables to look for based on the notation
    u               =   defs.u_fluctuations
    w               =   defs.w_fluctuations
    mrho_h2o_fluc   =   defs.h2o_molar_density_fluctuations
    theta_fluc      =   defs.thermodyn_temp_fluctuations
    theta_v_fluc    =   defs.virtual_temp_fluctuations
    q_fluc          =   defs.specific_humidity_fluctuations
    solutesf        = [ defs.molar_density % defs.fluctuations % solute for solute in solutes ]
    #---------

    #---------
    # Now we try to calculate or identify the fluctuations of theta
    theta_mean = data[ defs.thermodyn_temp ].mean()
    if theta_fluc not in data.columns:
        print('Fluctuations of theta not found. Will try to calculate it ... ', end='')
        #---------
        # We need the mean of the specific humidity and temperature
        data_q_mean =   data[ defs.specific_humidity ].mean()
        data[ defs.thermodyn_temp_fluctuations ] = (data[theta_v_fluc] - 0.61*theta_mean*data[q_fluc])/(1.+0.61*data_q_mean)
        print('done!')
        #---------
    #---------

    #---------
    # First we construct the covariance matrix (slower but more readable than doing it separately
    cov = data[[u, w, theta_v_fluc, mrho_h2o_fluc, theta_fluc] + solutesf ].cov()
    #---------

    #---------
    # Define auxiliar variables
    mu = constants.R_spec['h2o']/constants.R_spec['dry']
    rho_air_mean    =   data[ defs.moist_air_density ].mean()
    rho_dry_mean    =   data[ defs.dry_air_density ].mean()
    #---------

    #---------
    # Calculate the fluxes
    out = pd.DataFrame(index=[data.index[0]])
    out[ defs.momentum_flux ]               = -rho_air_mean * cov[u][w]
    out[ defs.sensible_heat_flux ]          = rho_air_mean * cp * cov[theta_fluc][w]
    out[ defs.virtual_sensible_heat_flux ]  = rho_air_mean * cp * cov[theta_v_fluc][w]
    out[ defs.water_vapor_flux ]            = cov[mrho_h2o_fluc][w]
    out[ defs.latent_heat_flux ]            = lamb( theta_mean ) * cov[mrho_h2o_fluc][w]
    #---------

    #-----------------
    # And then the flux units
    fluxunits = {}
    fluxunits[ defs.sensible_heat_flux ] = units[ defs.moist_air_density ] * cunits['cp_water'] * units['theta']*units[ w ]
    fluxunits[ defs.virtual_sensible_heat_flux ] = units[ defs.moist_air_density ] * cunits['cp_water'] * units[ defs.virtual_temp ]*units[ w ]
    fluxunits[ defs.water_vapor_flux ] = units[ defs.molar_density % defs.h2o ]*units[ w ]
    fluxunits[ defs.latent_heat_flux ] = cunits[ 'latent_heat_water' ]*units[ defs.molar_density % defs.h2o ]*units[ w ]
    #-----------------

    #---------
    # Calculate flux for each solute
    for solute, solutef in zip(solutes, solutesf):
        out[ defs.flux_of % solute ] =  cov[solutef][w]
        fluxunits[ defs.flux_of % solute ] = units[ defs.molar_density % solute ]*units['w']
    #---------

    #------------------------
    # APPLY WPL CORRECTION. PAGES 34-35 OF MICRABORDA
    if wpl:
        rho_h2o_mean = data[ defs.h2o_density ].mean()
        #---------
        # If water vapor mixing ratio is present, use it. Otherwise we try to calculate it
        if defs.h2o_mixing_ratio in data.columns:
            r_h2o = data[ defs.h2o_mixing_ratio ].mean()
        else:
            r_h2o = rho_h2o_mean/rho_dry_mean
        #---------
        
        out[ defs.water_vapor_flux ] = (1. +mu*r_h2o)*( out[ defs.water_vapor_flux ] + rho_h2o_mean * (cov[theta_fluc][w]/data_theta_mean) )
        for solute in solutes:
            rho_sol_mean    =   data[ defs.molar_density % defs.co2 ].mean()
            rc=rho_sol_mean/rho_dry_mean
            out[ defs.flux_of % solute ] = \
                out[ defs.flux_of % solute ] + \
                rho_co2_mean*(1. + mu*r_h2o)*(cov[theta_fluc][w])/data_theta_mean + \
                mu*rc*out[ defs.water_vapor_flux ]
    #------------------------

    #-----------------
    # Here we calculate the flux units

    #fluxunits['LE'] = fluxunits['LE'] * pm.constants.units['molar_mass']
    #coef = (fluxunits['LE'].to(pm.ureg('watts/meter**2'))).magnitude
    #fluxunits['LE'] = pm.ureg('watts/meter**2')
    out.loc[:, 'LE'] = algs.convert_to(out['LE'], fluxunits, 'watts/meter**2', inplace=True, key='LE')
    out.loc[:, 'Hv'] = algs.convert_to(out['Hv'], fluxunits, 'watts/meter**2', inplace=True, key='Hv')
    out.loc[:, 'H'] = algs.convert_to(out['H'], fluxunits, 'watts/meter**2', inplace=True, key='H')

    #-----------------

    units.update(fluxunits)
    return out, units



