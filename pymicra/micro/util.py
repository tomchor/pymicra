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

def preProcess(data, units, notation=None, use_means=False,
        rho_air_from_theta_v=True, inplace=True, theta=None, theta_unit=None, solutes=[]):
    '''
    Pre-processes data by calculating moist and dry air densities, specific humidity
    mass density and other important variables

    Parameters:
    -----------
    data: pandas.DataFrame
        dataframe with micrometeorological measurements
    units: dict
        units dictionary with the columns of data as keys
    notation: pymicra.notation
        defining notation used in data
    rho_air_from_theta_v: bool
        whether to use theta_v to calculate air density or theta
    inplace: bool
        treat units inplace or not
    theta: pandas.Series
        auxiliar theta measurement to be used if rho_air_from_theta_v==False
    theta_unit: pint.quantity
        auxiliar theta measurement's unit to be used if rho_air_from_theta_v==False
    solutes: list
        list of string where each string is a solute to be considered

    Returns:
    --------
    data: pandas.DataFrame
        dataframe with original columns and new calculated ones
    '''
    from .. import constants
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
    defs = algs.get_notation(notation)
    #---------

    print('Beginning of pre-processing ...')
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
    if defs.h2o_mass_density not in data.columns:
        print("Didn't locate mass density of h2o. Trying to calculate it ... ", end='')
        data.loc[:, defs.h2o_mass_density ] = data.loc[:, defs.h2o_molar_density ]*Mh2o
        units.update({ defs.h2o_mass_density : units[ defs.h2o_molar_density ]*molar_mass_unit })
        print("Done!")
    data.loc[:, defs.h2o_mass_density] = algs.convert_to(data[ defs.h2o_mass_density ], units, 'kg/m**3', inplace=True, key=defs.h2o_mass_density)
    #---------

    #---------
    # Check for h2o molar density
    if defs.h2o_molar_density not in data.columns:
        print("Didn't locate molar density of h2o. Trying to calculate it ... ", end='')
        data.loc[:, defs.h2o_molar_density ] = data.loc[:, defs.h2o_mass_density ]/Mh2o
        units.update({ defs.h2o_molar_density : units[ defs.h2o_mass_density ]/molar_mass_unit })
        print("Done!")
    #---------

    #---------
    # Calculation of rho_air is done here
    if (defs.moist_air_density not in data.columns):
        print('Moist air density not present in dataset')
        if rho_air_from_theta_v:
            print('Calculating rho_air = p/(Rdry * theta_v) ... ', end='')
            data = physics.airDensity_from_theta_v(data, units, notation=defs, inplace=True, use_means=use_means)
        else:
            if theta:
                print('Trying to calculate rho_air using auxiliar theta measurement ... ', end='')
                data = physics.airDensity_from_theta(data, units, notation=defs, inplace=True, use_means=use_means, theta=theta, theta_unit=theta_unit)
            else:
                print('Trying to calculate rho_air using theta from this dataset ... ', end='')
                data = physics.airDensity_from_theta(data, units, notation=defs, inplace=True, use_means=use_means, theta=None)
        print('Done!')
    #---------

    #---------
    # Calculation of dry air mass density is done here
    if (defs.dry_air_mass_density not in data.columns):
        print('Calculating dry_air mass_density = rho_air - rho_h2o ... ', end='')
        data.loc[:, defs.dry_air_mass_density ] = algs.add([ data[defs.moist_air_mass_density], -data[defs.h2o_mass_density] ], 
                        [ units[defs.moist_air_mass_density], units[defs.h2o_mass_density] ], inplace=True, unitdict=units, key=defs.dry_air_mass_density)
        print('Done!')
    #---------

    #---------
    # Calculation of dry air molar density is done here
    if (defs.dry_air_molar_density not in data.columns):
        print('Dry air molar density not in dataset')
        if defs.dry_air_mass_density not in data.columns:
            print("Can't calculate it. Dry air mass density not present")
        print('Calculating dry_air molar_density = rho_dry / dry_air_molar_mass ... ', end='')
        data.loc[:, defs.dry_air_molar_density ] = data[ defs.dry_air_mass_density ]/constants.molar_mass['dry']
        units.update({ defs.dry_air_molar_density : (units[ defs.dry_air_mass_density ]/molar_mass_unit).to_base_units() })
        print('Done!')
    #---------
 
    #---------
    # Calculation of specific humidity is done here
    if (defs.specific_humidity not in data.columns):
        print('Calculating specific humidity = rho_h2o / rho_air ... ', end='')
        data.loc[:, defs.specific_humidity] = data[ defs.h2o_density ] / data[ defs.moist_air_density ]
        units.update({ defs.specific_humidity : units[ defs.h2o_density ] / units[ defs.moist_air_density ] })
        print('Done!')
    #---------

    #---------
    # Calculation of h2o mass mixing ratio is done here
    if (defs.h2o_mass_mixing_ratio not in data.columns):
        print('Calculating h2o mass mixing ratio = rho_h2o / rho_dry ... ', end='')
        data.loc[:, defs.h2o_mass_mixing_ratio] = data[ defs.h2o_mass_density ] / data[ defs.dry_air_mass_density ]
        units.update({ defs.h2o_mass_mixing_ratio : units[ defs.h2o_mass_density ] / units[ defs.dry_air_mass_density ] })
        print('Done!')
    #---------

    #---------
    # Calculation of h2o mass mixing ratio is done here
    if (defs.h2o_molar_mixing_ratio not in data.columns):
        print('Calculating h2o mass mixing ratio = rho_h2o / rho_dry ... ', end='')
        data.loc[:, defs.h2o_molar_mixing_ratio] = data[ defs.h2o_molar_density ] / data[ defs.dry_air_molar_density ]
        units.update({ defs.h2o_molar_mixing_ratio : units[ defs.h2o_molar_density ] / units[ defs.dry_air_molar_density ] })
        print('Done!')
    #---------

    #---------
    # Converting dry_air_molar_density and dry_air_molar_mixing_ratio to mol/meter**3 and mole/mole, respectively
    convert_to={defs.h2o_molar_mixing_ratio : 'mole/mole',
                defs.dry_air_molar_density : 'mole/meter**3'}
    data = data.convert_cols(convert_to, units, inplace=True)
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
        sol_mass_density = eval('defs.{}_mass_density'.format(solute))
        sol_molar_density = eval('defs.{}_molar_density'.format(solute))
        sol_mass_concentration = eval('defs.{}_mass_concentration'.format(solute))
        sol_mass_mixing_ratio = eval('defs.{}_mass_mixing_ratio'.format(solute))
        sol_molar_mixing_ratio = eval('defs.{}_molar_mixing_ratio'.format(solute))
        M_sol = constants.molar_mass[solute]

        #---------
        # Check for solute mass density
        if sol_mass_density not in data.columns:
            print("Didn't locate mass density of {}. Trying to calculate it ... ".format(solute), end='')
            data.loc[:, sol_mass_density ] = data.loc[:, sol_molar_density ]*M_sol
            units.update({ sol_mass_density : units[ sol_molar_density ]*molar_mass_unit })
            print("Done!")
        data.loc[:, sol_mass_density] = algs.convert_to(data[ sol_mass_density ], units, 'kg/m**3', inplace=True, key=sol_mass_density)
        #---------
    
        #---------
        # Check for solute molar density
        if sol_molar_density not in data.columns:
            print("Didn't locate molar density of {}. Trying to calculate it".format(solute))
            data.loc[:, sol_molar_density ] = data.loc[:, sol_mass_density ]/M_sol
            units.update({ sol_molar_density : units[ sol_mass_density ]/molar_mass_unit })
            print("Done!")
        #---------

        #---------
        # Calculation of SOLUTE MASS concentration (g/g) is done here
        if (sol_mass_concentration not in data.columns):
            print('Calculating {0} mass concentration (g/g) = rho_{0} / rho_air ... '.format(solute), end='')
            data.loc[:, sol_mass_concentration] = data[ sol_mass_density ] / data[ defs.moist_air_mass_density ]
            units.update({ sol_mass_concentration : units[ sol_mass_density ] / units[ defs.moist_air_mass_density ] })
            print("Done!")
        #---------
            
    
        #---------
        # Calculation of SOLUTE MASS mixing ratio is done here
        if (sol_mass_mixing_ratio not in data.columns):
            print('Calculating {0} mass mixing ratio = rho_{0} / rho_dry ... '.format(solute), end='')
            data.loc[:, sol_mass_mixing_ratio] = data[ sol_mass_density ] / data[ defs.dry_air_mass_density ]
            units.update({ sol_mass_mixing_ratio : units[ sol_mass_density ] / units[ defs.dry_air_mass_density ] })
            print("Done!")
        #---------

        #---------
        # Calculation of SOLUTE MOLAR mixing ratio is done here
        if (sol_molar_mixing_ratio not in data.columns):
            print('Calculating {0} molar mixing ratio = mrho_{0} / mrho_dry ... '.format(solute), end='')
            data.loc[:, sol_molar_mixing_ratio] = data[ sol_molar_density ] / data[ defs.dry_air_molar_density ]
            units.update({ sol_molar_mixing_ratio : units[ sol_molar_density ] / units[ defs.dry_air_molar_density ] })
            print("Done!")
        #---------

        #---------
        # Converting molar, mass mixing ratio and mass concentration to unitless
        convert_to={sol_molar_mixing_ratio : 'mole/mole',
                    sol_mass_mixing_ratio : 'g/g',
                    sol_mass_concentration : 'g/g'}
        data = data.convert_cols(convert_to, units, inplace=True)
        #---------
 
    #---------
        
    print('Pre-processing complete.\n')
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
    #from ..core import notation
    import pandas as pd
    from .. import algs

    if from_fluctuations:
        return eddyCov3(data, wpl=wpl, notation_defs=notation_defs, solutes=solutes)

    cp=constants.cp_dry
    lamb = constants.latent_heat_water

    #---------
    # Define useful notation to look for
    defs=algs.get_notation(notation_defs)
#    if notation_defs==None:
#        defs=notation()
#    else:
#        defs=notation_defs
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
    #from ..core import notation
    import pandas as pd
    from .. import algs

    cp = constants.cp_dry
    lamb = constants.latent_heat_water

    #---------
    # Define useful notation to look for
    defs = algs.get_notation(notation_defs)
#    if notation_defs==None:
#        defs=notation()
#    else:
#        defs=notation_defs
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
    units: dict
        units dictionary
    wpl: boolean
        whether or not to apply WPL correction on the latent heat flux and solutes flux
    notation: pymicra.notation
        object that holds the notation used in the dataframe
    inplace: bool
        whether or not to treat the units inplace
    solutes: list
        list that holds every solute considered for flux
    """
    from .. import constants
    from .. import algs
    import pandas as pd
    from .. import ureg

    cp = constants.cp_dry
    lamb = constants.latent_heat_water

    defs = algs.get_notation(notation)
    data = data.copy()
    if (not inplace) and units:
        units = units.copy()
    cunits = constants.units

    print('Beginning Eddy Covariance method...')
    #---------
    # Define name of variables to look for based on the notation
    u_fluc          =   defs.u_fluctuations
    w_fluc          =   defs.w_fluctuations
    mrho_h2o_fluc   =   defs.h2o_molar_density_fluctuations
    rho_h2o_fluc    =   defs.h2o_density_fluctuations
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
        if not (units[ theta_v_fluc ]==ureg['kelvin'] and units[ defs.thermodyn_temp ]==ureg['kelvin']):
            raise TypeError('Units for both the virtual temp fluctuations and the thermodynamic temperature must be Kelvin')
        data_q_mean =   data[ defs.specific_humidity ].mean()
        data[ theta_fluc ] = (data[theta_v_fluc] - 0.61*theta_mean*data[q_fluc])/(1.+0.61*data_q_mean)
        theta_fluc_unit = units[ theta_v_fluc ]
        print('done!')
        #---------
    #---------

    #---------
    # First we construct the covariance matrix (slower but more readable than doing it separately)
    # maybe figure out later a way that is both faster and more readable
    cov = data[[u_fluc, w_fluc, theta_v_fluc, mrho_h2o_fluc, rho_h2o_fluc, theta_fluc] + solutesf ].cov()
    #---------

    #---------
    # Define auxiliar variables
    rho_air_mean    =   data[ defs.moist_air_density ].mean()
    rho_dry_mean    =   data[ defs.dry_air_density ].mean()
    #---------

    #---------
    # Calculate the fluxes
    print('Calculating fluxes ... ', end='')
    out = pd.DataFrame(index=[ data.index[0] ])
    out[ defs.momentum_flux ]               = -rho_air_mean * cov[ u_fluc ][ w_fluc ]
    out[ defs.sensible_heat_flux ]          = rho_air_mean * cp * cov[theta_fluc][w_fluc]
    out[ defs.virtual_sensible_heat_flux ]  = rho_air_mean * cp * cov[theta_v_fluc][w_fluc]
    out[ defs.water_vapor_flux ]            = cov[mrho_h2o_fluc][w_fluc]
    out[ defs.latent_heat_flux ]            = lamb( theta_mean ) * cov[rho_h2o_fluc][w_fluc]
    print('done!')
    #---------

    #-----------------
    # And then the flux units
    fluxunits = {}
    fluxunits[ defs.momentum_flux ]         = units[ defs.moist_air_density ] * units[ u_fluc ]*units[ w_fluc ]
    fluxunits[ defs.sensible_heat_flux ]    = units[ defs.moist_air_density ] * cunits['cp_water'] * theta_fluc_unit * units[ w_fluc ]
    fluxunits[ defs.virtual_sensible_heat_flux ] = units[ defs.moist_air_density ] * cunits['cp_water'] * units[ theta_v_fluc ]*units[ w_fluc ]
    fluxunits[ defs.water_vapor_flux ]      = units[ defs.h2o_molar_density ]*units[ w_fluc ]
    fluxunits[ defs.latent_heat_flux ]      = cunits[ 'latent_heat_water' ]*units[ rho_h2o_fluc ]*units[ w_fluc ]
    #-----------------

    #---------
    # Calculate flux for each solute
    for solute, solutef in zip(solutes, solutesf):
        out[ defs.flux_of % solute ] =  cov[ solutef ][ w_fluc ]
        fluxunits[ defs.flux_of % solute ] = units[ solutef ]*units[ w_fluc ]
    #---------

    #------------------------
    # APPLY WPL CORRECTION. PAGES 34-35 OF MICRABORDA
    if wpl:
        print('Applying WPL correction for water vapor flux ... ', end='')
        mu = constants.R_spec['h2o']/constants.R_spec['dry']
        mrho_h2o_mean = data[ defs.h2o_molar_density ].mean()

        #---------
        # If there are solutes to correct we save the original E to use in the calculation
        if solutes:
            E_orig = out[ defs.water_vapor_flux ].copy()
            E_orig_unit = fluxunits[ defs.water_vapor_flux ].copy()
        #---------

        #---------
        # If water vapor mixing ratio is present, use it. Otherwise we try to calculate it
        if defs.h2o_molar_mixing_ratio in data.columns:
            mr_h2o = data[ defs.h2o_molar_mixing_ratio ].mean()
        elif units[ defs.h2o_molar_density ] == units[ defs.dry_air_molar_density ]:
            mr_h2o = mrho_h2o_mean/mrho_dry_medan
        else:
            raise TypeError('Either water mixing ratio should be provided, or dry air and water density should be the same')
        #---------
        
        #---------
        # We calculate WPL little by little to make it easy to handle the units
        aux1 = out[ defs.water_vapor_flux ]
        unt1 = fluxunits[ defs.water_vapor_flux ]

        aux2 = mrho_h2o_mean * (cov[theta_fluc][w_fluc]/theta_mean)
        unt2 = units[ defs.h2o_molar_density ] * units[ w_fluc ]

        aux3, unt3 = algs.add([aux1, aux2], [unt1, unt2], inplace=False)

        out.loc[:, defs.water_vapor_flux ] = (1. +mu*mr_h2o)* aux3
        fluxunits[ defs.water_vapor_flux ] = unt3
        print('done!')
        #---------

        #---------
        # Now we re-calculate LE based on the corrected E
        print('Applying WPL correction for latent heat flux ... ', end='')
        out.loc[:, defs.latent_heat_flux ] = lamb(theta_mean) * out[ defs.water_vapor_flux ] * constants.molar_mass['h2o']
        fluxunits[ defs.latent_heat_flux ] = cunits[ 'latent_heat_water' ] * fluxunits[ defs.water_vapor_flux ] * cunits['molar_mass']
        print('done!')
        #---------

        #---------
        # NEEDS IMPROVEMENT
        #---------
        for solute in solutes:
            print('Applying WPL correction for {} ... '.format(solute), end='')
            sol_molar_density_mean    =   data[ defs.molar_density % solute ].mean()

            #---------
            # Obtaining solute mass mixing ratio
            if defs.molar_mixing_ratio % solute in data.columns:
                mr_sol = data[ defs.molar_mixing_ratio % solute ].mean()
                mr_sol_unit = units[ defs.molar_mixing_ratio % solute ]
            elif units[ defs.molar_density % solute ] == units[ defs.dry_air_molar_density ]:
                mr_sol = sol_molar_density_mean/rho_dry_mean
                mr_sol_unit = 'mole/mole'
            else:
                raise TypeError('Either water mixing ratio should be provided, or dry air and water density should be the same')
            #---------

            #---------
            # Again, we calculate WPL little by little to make it easy to handle the units
            aux1 = out[ defs.flux_of % solute ]
            unt1 = fluxunits[ defs.flux_of % solute ]

            aux2 = sol_molar_density_mean*(1. + mu*mr_h2o)*(cov[theta_fluc][w_fluc])/theta_mean
            unt2 = units[ defs.molar_density % solute ] * units[ w_fluc ]

            aux3 = mu * mr_sol* E_orig
            unt3 = mr_sol_unit * E_orig_unit
 
            out.loc[:, defs.flux_of % solute ] = algs.add([aux1, aux2, aux3], [unt1, unt2, unt3],
                                                    inplace=True, unitdict=fluxunits, key=defs.flux_of % solute)
            print('done!')
            #---------
    #------------

    #------------
    # Here we convert the units to watts/m2
    convert_to ={defs.latent_heat_flux : 'watts/meter**2',
                defs.virtual_sensible_heat_flux : 'watts/meter**2',
                defs.sensible_heat_flux : 'watts/meter**2'}
    out = out.convert_cols(convert_to, fluxunits, inplace=True)
    #------------

    print('Done with Eddy Covariance.\n')
    if inplace:
        units.update(fluxunits)
        return out
    else:
        return out, fluxunits



def fluxes_from_scales(data, units, wpl=True,
        notation=None, inplace=True, solutes=[]):
    """
    Get fluxes from the turbulent scales
    
    WARNING! If using wpl=True, be sure that all masses are consistent!
        For example, if q = [g/g], rho_h2o = [g/m3] and rho_co2 = [g/m3] and so on.
        Avoid mixing kg/m3 with g/m3 (e.g. for co2 and h2o) and mg/kg with g/g (e.g. for
        co2 and h2o).

    Parameters:
    -----------
    data: pandas.DataFrame
        dataframe with the characteristic lengths calculated
    units: dict
        dictionary with variable and their units
    notation: pymicra.notation
        object that holds the notation used in the dataframe
    wpl: boolean
        whether or not to apply WPL correction on the latent heat flux and solutes flux
    inplace: bool
        whether or not to treat units inplace
    solutes: list
        list that holds every solute considered for flux
    """
    from .. import constants
    from .. import algs
    import pandas as pd
    from .. import ureg

    cp = constants.cp_dry
    lamb = constants.latent_heat_water

    defs = algs.get_notation(notation)
    data = data.copy()
    if not inplace:
        units = units.copy()
    cunits = constants.units

    print('Beginning retrieval of fluxes ...')
    #---------
    # Define name of variables to look for based on the notation
    u_star          =   defs.u_star
    mrho_h2o_star   =   defs.h2o_molar_density_star
    rho_h2o_star    =   defs.h2o_density_star
    theta_star      =   defs.thermodyn_temp_star
    theta_v_star    =   defs.virtual_temp_star
    q_star          =   defs.specific_humidity_star
    solute_stars    = [ defs.molar_density % defs.star % solute for solute in solutes ]
    #---------

    #---------
    # Now we try to calculate or identify the turbulent scale of theta
    theta_mean = data[ defs.thermodyn_temp ].mean()
    if theta_star not in data.columns:
        print('Fluctuations of theta not found. Will try to calculate it ... ', end='')
        #---------
        # We need the mean of the specific humidity and temperature
        if not (units[ theta_v_star ]==ureg['kelvin'] and units[ defs.thermodyn_temp ]==ureg['kelvin']):
            raise TypeError('Units for both the virtual temp and the thermodynamic temperature turbulent scales must be the same')
        data_q_mean =   data[ defs.specific_humidity ].mean()
        data[ theta_star ] = (data[theta_v_star] - 0.61*theta_mean*data[q_star])/(1.+0.61*data_q_mean)
        theta_star_unit = units[ theta_v_star ]
        print('done!')
        #---------
    #---------

    #---------
    # Define auxiliar variables
    rho_air_mean    =   data[ defs.moist_air_density ].mean()
    rho_dry_mean    =   data[ defs.dry_air_density ].mean()
    #---------

    #---------
    # Calculate the fluxes
    print('Calculating fluxes ... ', end='')
    out = pd.DataFrame(index=[ data.index[0] ])
    out[ defs.momentum_flux ]               = -rho_air_mean * data[ u_star ] * data[ u_star ]
    out[ defs.sensible_heat_flux ]          = rho_air_mean * cp * data[ theta_star ] * data[ u_star ]
    out[ defs.virtual_sensible_heat_flux ]  = rho_air_mean * cp * data[ theta_v_star] * data[ u_star ]
    out[ defs.water_vapor_flux ]            = data[ mrho_h2o_star ] * data[ u_star ]
    out[ defs.latent_heat_flux ]            = lamb(theta_mean) * data[ rho_h2o_star ] * data[ u_star ]
    print('done!')
    #---------

    #-----------------
    # And then the flux units
    fluxunits = {}
    fluxunits[ defs.momentum_flux ]         = units[ defs.moist_air_density ] * units[ u_fluc ]*units[ w_fluc ]
    fluxunits[ defs.sensible_heat_flux ]    = units[ defs.moist_air_density ] * cunits['cp_water'] * theta_fluc_unit * units[ w_fluc ]
    fluxunits[ defs.virtual_sensible_heat_flux ] = units[ defs.moist_air_density ] * cunits['cp_water'] * units[ theta_v_fluc ]*units[ w_fluc ]
    fluxunits[ defs.water_vapor_flux ]      = units[ defs.h2o_molar_density ]*units[ w_fluc ]
    fluxunits[ defs.latent_heat_flux ]      = cunits[ 'latent_heat_water' ]*units[ rho_h2o_fluc ]*units[ w_fluc ]
    #-----------------

    #---------
    # Calculate flux for each solute
    for solute, solutef in zip(solutes, solutesf):
        out[ defs.flux_of % solute ] =  cov[ solutef ][ w_fluc ]
        fluxunits[ defs.flux_of % solute ] = units[ solutef ]*units[ w_fluc ]
    #---------

    #------------------------
    # APPLY WPL CORRECTION. PAGES 34-35 OF MICRABORDA
    if wpl:
        print('Applying WPL correction for water vapor ... ', end='')
        mu = constants.R_spec['h2o']/constants.R_spec['dry']
        rho_h2o_mean = data[ defs.h2o_density ].mean()

        #---------
        # If water vapor mixing ratio is present, use it. Otherwise we try to calculate it
        if defs.h2o_mixing_ratio in data.columns:
            r_h2o = data[ defs.h2o_mixing_ratio ].mean()
        elif units[ defs.h2o_density ] == units[ defs.dry_air_density ]:
            r_h2o = rho_h2o_mean/rho_dry_mean
        else:
            raise TypeError('Either water mixing ratio should be provided, or dry air and water density should be the same')
        #---------
        
        out.loc[:, defs.water_vapor_flux ] = (1. +mu*r_h2o)*( out[ defs.water_vapor_flux ] + rho_h2o_mean * (cov[theta_fluc][w_fluc]/theta_mean) )
        print('done!')

        for solute in solutes:
            print('Applying WPL correction for {} ... '.format(solute), end='')
            #---------
            # Check if mixing ratio units are correct
            if units[ defs.density % solute ] == units[ defs.dry_air_density ]:
                rho_sol_mean    =   data[ defs.molar_density % solute ].mean()
                rc = rho_sol_mean/rho_dry_mean
            else:
                raise TypeError('Either water mixing ratio should be provided, or dry air and water density should be the same')
            #---------

            out.loc[:, defs.flux_of % solute ] = \
                out[ defs.flux_of % solute ] + \
                rho_sol_mean*(1. + mu*r_h2o)*(cov[theta_fluc][w_fluc])/theta_mean + \
                mu*rc*out[ defs.water_vapor_flux ]
            print('done!')
    #------------------------

    #-----------------
    # Here we convert the units to watts/m2
    out.loc[:, defs.latent_heat_flux ] = algs.convert_to(out[ defs.latent_heat_flux ], fluxunits, 'watts/meter**2', inplace=True, key='LE')
    out.loc[:, defs.virtual_sensible_heat_flux ] = algs.convert_to(out['Hv'], fluxunits, 'watts/meter**2', inplace=True, key='Hv')
    out.loc[:, defs.sensible_heat_flux ] = algs.convert_to(out['H'], fluxunits, 'watts/meter**2', inplace=True, key='H')
    #-----------------

    print('Done with Eddy Covariance.\n')
    if inplace:
        units.update(fluxunits)
        return out
    else:
        return out, fluxunits



