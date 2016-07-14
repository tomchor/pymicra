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
    # First convert any temperature if it is still in Celsius
    temps = { col:'kelvin' for col in data.columns if col in [defs.thermodyn_temp, defs.virtual_temp, defs.sonic_temp, defs.potential_temp] }
    print('Converting {} to kelvin ... '.format(' and '.join(temps.keys())), end='')
    data = data.convert_cols(temps, units, inplace=True)
    print("Done!")
    #---------

    #---------
    # Check for h2o mass density
    if defs.h2o_mass_density not in data.columns:
        print("Didn't locate mass density of h2o. Trying to calculate it ... ", end='')
        data.loc[:, defs.h2o_mass_density ] = data.loc[:, defs.h2o_molar_density ]*Mh2o
        units.update({ defs.h2o_mass_density : units[ defs.h2o_molar_density ]*molar_mass_unit })
        print("Done!")
    #data.loc[:, defs.h2o_mass_density] = algs.convert_to(data[ defs.h2o_mass_density ], units, 'kg/m**3', inplace=True, key=defs.h2o_mass_density)
    data = data.convert_cols({defs.h2o_mass_density:'kg/m**3'}, units, inplace=True)
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
    if (defs.moist_air_mass_density not in data.columns):
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
        data.loc[:, defs.specific_humidity] = data[ defs.h2o_mass_density ] / data[ defs.moist_air_mass_density ]
        units.update({ defs.specific_humidity : units[ defs.h2o_mass_density ] / units[ defs.moist_air_mass_density ] })
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
        data = data.convert_cols({sol_mass_density:'kg/m**3'}, units, inplace=True)
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


def eddyCovariance(data, units, wpl=True, get_turbulent_scales=True, site_config=None, output_as_df=True,
        notation=None, theta_fluct_from_theta_v=True, inplace=True, solutes=[]):
    """
    Get fluxes from the turbulent fluctuations
    
    Parameters:
    -----------
    data: pandas.DataFrame
        dataframe with the characteristic lengths calculated
    units: dict
        units dictionary
    wpl: boolean
        whether or not to apply WPL correction on the latent heat flux and solutes flux
    get_turbulent_scales: bool
        whether or not to use getScales to return turbulent scales
    site_config: pymica.siteConfig
        siteConfig object to pass to getScales if get_turbulent_scales==True
    notation: pymicra.Notation
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
    defsdic = defs.__dict__

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
    rho_h2o_fluc    =   defs.h2o_mass_density_fluctuations
    theta_fluc      =   defs.thermodyn_temp_fluctuations
    theta_v_fluc    =   defs.virtual_temp_fluctuations
    q_fluc          =   defs.specific_humidity_fluctuations
    solutesf        = [ defsdic['%s_molar_density_fluctuations' % solute] for solute in solutes ]
    solutefluxes    = [ defsdic[ '%s_flux' % solute ] for solute in solutes ]
    #---------

    #---------
    # Now we try to calculate or identify the fluctuations of theta
    theta_mean = data[ defs.thermodyn_temp ].mean()
    if (theta_fluc not in data.columns) or theta_fluct_from_theta_v:
        print('Fluctuations of theta not found. Will try to calculate it ... ', end='')
        #---------
        # We check the units of theta_v and theta
        if not (units[ theta_v_fluc ]==ureg['kelvin'] and units[ defs.thermodyn_temp ]==ureg['kelvin']):
            raise TypeError('Units for both the virtual temp fluctuations and the thermodynamic temperature must be Kelvin')
        #---------

        #---------
        # Estimate theta fluctuations from theta_v
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
    rho_air_mean    =   data[ defs.moist_air_mass_density ].mean()
    rho_dry_mean    =   data[ defs.dry_air_mass_density ].mean()
    #---------

    #---------
    # Calculate the fluxes
    print('Calculating fluxes ... ', end='')
    idx0 = data.index[0]
    out = pd.Series(name=idx0)
    out[ defs.momentum_flux ]               = -rho_air_mean * cov[ u_fluc ][ w_fluc ]
    out[ defs.sensible_heat_flux ]          = rho_air_mean * cp * cov[theta_fluc][w_fluc]
    out[ defs.virtual_sensible_heat_flux ]  = rho_air_mean * cp * cov[theta_v_fluc][w_fluc]
    out[ defs.water_vapor_flux ]            = cov[mrho_h2o_fluc][w_fluc]
    out[ defs.latent_heat_flux ]            = lamb(theta_mean) * cov[rho_h2o_fluc][w_fluc]
    print('done!')
    #---------

    #-----------------
    # And then the flux units
    fluxunits = {}
    fluxunits[ defs.momentum_flux ]         = units[ defs.moist_air_mass_density ] * units[ u_fluc ]*units[ w_fluc ]
    fluxunits[ defs.sensible_heat_flux ]    = units[ defs.moist_air_mass_density ] * cunits['cp_water'] * theta_fluc_unit * units[ w_fluc ]
    fluxunits[ defs.virtual_sensible_heat_flux ] = units[ defs.moist_air_mass_density ] * cunits['cp_water'] * units[ theta_v_fluc ]*units[ w_fluc ]
    fluxunits[ defs.water_vapor_flux ]      = units[ defs.h2o_molar_density ]*units[ w_fluc ]
    fluxunits[ defs.latent_heat_flux ]      = cunits[ 'latent_heat_water' ]*units[ rho_h2o_fluc ]*units[ w_fluc ]
    #-----------------

    #---------
    # Calculate flux for each solute
    for sol_flux, solute, solutef in zip(solutefluxes, solutes, solutesf):
        out[ sol_flux ] =  cov[ solutef ][ w_fluc ]
        fluxunits[ sol_flux ] = units[ solutef ]*units[ w_fluc ]
    #---------

    #------------------------
    # APPLY WPL CORRECTION. PAGES 34-35 OF MICRABORDA
    if wpl:
        print('Applying WPL correction for water vapor flux ... ', end='\n')
        mu = constants.R_spec['h2o']/constants.R_spec['dry']
        mrho_h2o_mean = data[ defs.h2o_molar_density ].mean()

        #---------
        # If there are solutes to correct we save the original E to use in the calculation
        E_orig = out[ defs.water_vapor_flux ].copy()
        E_orig_unit = fluxunits[ defs.water_vapor_flux ]
        #---------

        #---------
        # If water vapor mixing ratio is present, use it. Otherwise we try to calculate it
        if defs.h2o_molar_mixing_ratio in data.columns:
            mr_h2o = data[ defs.h2o_molar_mixing_ratio ].mean()
            mr_h2o = (mr_h2o * units[ defs.h2o_molar_mixing_ratio ]).to('dimensionless').magnitude

        elif defs.dry_air_molar_density in data.columns:
            mr_h2o = (mrho_h2o_mean*units[ defs.h2o_molar_density ] / 
                    (data[ defs.dry_air_molar_density ].mean()*units[ defs.dry_air_molar_density ])).to('dimensionless').magnitude

        else:
            raise TypeError('Either water molar mixing ratio should be provided, or dry air and water density should be the same')
        #---------
        
        #---------
        # We calculate WPL little by little to make it easy to handle the units
        aux1 = E_orig
        unt1 = E_orig_unit

        aux2 = mrho_h2o_mean * (cov[theta_fluc][w_fluc]/theta_mean)
        unt2 = units[ defs.h2o_molar_density ] * units[ w_fluc ]

        aux3 = aux1*unt1 + aux2*unt2
        unt3 = aux3/aux3.magnitude
        aux3 = aux3.magnitude

        out.loc[ defs.water_vapor_flux ] = (1. +mu*mr_h2o)* aux3
        fluxunits[ defs.water_vapor_flux ] = unt3
        #---------

        #---------
        # Now we re-calculate LE based on the corrected E
        print('Applying WPL correction for latent heat flux ... ', end='\n')
        out[ defs.latent_heat_flux ] = lamb(theta_mean) * out[ defs.water_vapor_flux ] * constants.molar_mass['h2o']
        fluxunits[ defs.latent_heat_flux ] = cunits[ 'latent_heat_water' ] * fluxunits[ defs.water_vapor_flux ] * cunits['molar_mass']
        #---------

        #---------
        # Recalculating covs with WPL
        print("Re-calculating cov(%s, w') according to WPL correction ... " % mrho_h2o_fluc, end='')
        wplcov = cov.copy()
        w_h2o_units = units[ w_fluc ] * units[ mrho_h2o_fluc ]
        wplcov.loc[ mrho_h2o_fluc, w_fluc ] = (out[ defs.water_vapor_flux ] * fluxunits[ defs.water_vapor_flux ]).to(w_h2o_units).magnitude
        wplcov.loc[ w_fluc, mrho_h2o_fluc ] = wplcov.loc[ mrho_h2o_fluc, w_fluc ]
        print('done!')
        #---------

        #---------
        # NEEDS IMPROVEMENT
        for sol_flux, solutef, solute in zip(solutefluxes, solutesf, solutes):

            print('Applying WPL correction for {} ... '.format(sol_flux), end='\n')
            sol_molar_density_mean    =   data[ defs.molar_density % solute ].mean()

            #---------
            # Obtaining solute mass mixing ratio
            if defs.molar_mixing_ratio % solute in data.columns:
                mr_sol = data[ defs.molar_mixing_ratio % solute ].mean()
                mr_sol_unit = units[ defs.molar_mixing_ratio % solute ]

            elif defs.dry_air_molar_density in data.columns:
                mr_sol = sol_molar_density_mean/data[ defs.dry_air_molar_density ]
                mr_sol_unit = units[ defs.molar_density % solute ] / units[ defs.dry_air_molar_density ]

            else:
                raise TypeError('Either water mixing ratio should be provided, or dry air molar density')
            #---------

            #---------
            # Again, we calculate WPL little by little to make it easy to handle the units
            aux1 = out[ defs.flux_of % solute ]
            unt1 = fluxunits[ defs.flux_of % solute ]

            aux2 = sol_molar_density_mean * (1. + mu*mr_h2o) * (cov[theta_fluc][w_fluc])/theta_mean
            unt2 = units[ defs.molar_density % solute ] * units[ w_fluc ]

            aux3 = mu * mr_sol* E_orig
            unt3 = mr_sol_unit * E_orig_unit
 
            aux4 = aux1*unt1 + aux2*unt2 + aux3*unt3
            unt4 = aux4/aux4.magnitude
            aux4 = aux4.magnitude

            out[ sol_flux ] = aux4
            fluxunits[ sol_flux ] = unt4

            print("Re-calculating cov(%s, w') according to WPL correction ... " % solutef, end='')
            w_sol_units = units[ w_fluc ] * units[ solutef ]
            wplcov.loc[ solutef, w_fluc ] = (out.loc[ sol_flux ] * fluxunits[ sol_flux ]).to(w_sol_units).magnitude
            wplcov.loc[ w_fluc, solutef ] = wplcov.loc[ solutef, w_fluc ]
            print('done!')
            #---------
        #---------
    #------------

    #------------
    # Here we convert the units to watts/m**2
    convert_to ={defs.latent_heat_flux : 'watts/meter**2',
                defs.virtual_sensible_heat_flux : 'watts/meter**2',
                defs.sensible_heat_flux : 'watts/meter**2'}
    out = out.convert_indexes(convert_to, fluxunits, inplace=True)
    #------------

    #------------
    # We calculate the turbulent scales
    if get_turbulent_scales:
        assert site_config is not None, 'Must provide site_config keyword if get_turbulent_scales==True'
        from ..micro import turbulentScales
        theta_v_mean = data[ defs.virtual_temp ].mean()
        scales, scaleunits = turbulentScales(cov, site_config, units, notation=defs, inplace=False, 
                                solutes=solutes, theta_v_mean=theta_v_mean, output_as_df=False)
        out = pd.concat([out, scales])
        fluxunits.update(scaleunits)
    #------------

    #------------
    # Crate a one-row dataframe if output_as_df is True
    if output_as_df:
        out = out.to_frame().T
    #------------

    print('Done with Eddy Covariance.\n')
    if inplace:
        units.update(fluxunits)
        return out
    else:
        return out, fluxunits


def rotateCoor(data, notation=None, how='2d'):
    """
    """
    from .. import data as pmdata

    if how=='2d':
        return pmdata.rotate2D(data, notation=notation, how=how)
    else:
        return None


