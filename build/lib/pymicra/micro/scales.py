from __future__ import absolute_import, print_function, division

def MonObuVar(L_m, siteConf):
    """Redirects to stabilityParam()"""
    return stabilityParam(L_m, siteConf)

def stabilityParam(L_m, siteConf):
    """
    Calculates the Monin-Obukhov Similarity Variable
    defined as

    zeta = (z-d)/Lo
    where d is the displacement height or zero-plane displacement
    and L_m is the Monin-Obukhov Length.
    """

    z=siteConf.measurement_height
    d=siteConf.displacement_height
    zeta = (z-d)/L_m
    return zeta

def MonObuLen(theta_v_star, theta_v_mean, u_star, g=None):
    """Redirects to obukhovLen()"""
    return obukhovLen(theta_v_star, theta_v_mean, u_star)


def obukhovLen(data, units, theta_v_mean=None, theta_v_mean_unit=None, notation=None, inplace_units=True):
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
    from .. import constants
    from .. import algs

    defs = algs.get_notation(notation)
    data = data.copy()

    g = constants.gravity
    kappa = constants.kappa

    if not theta_v_mean:
        if not theta_v_mean_unit:
            raise ValueError('Must provide theta_v_mean_unit keyword if theta_v_mean is provided')
        theta_v_mean = data[ defs.mean_virtual_temperature ]

    num = (data[ defs.u_star ]**2.)*theta_v_mean
    num_unit = (units[ defs.u_star ]**2.)*theta_v_mean_unit

    denom = kappa*g*data[ defs.virtual_temperature_star ]
    denom_unit = constants.units['gravity']*units[ defs.virtual_temperature_star ]

    Lo = - num/denom
    Lo_unit = num_unit/denom_unit

    if inplace_units:
        units.update({defs.obukhov_length:Lo_unit})
        return Lo
    else:
        return Lo, Lo_unit



def turbulentScales(data, siteConf, units, notation=None, theta_v_mean=None, theta_v_mean_unit=None,
        theta_fluct_from_theta_v=True, solutes=[], output_as_df=True, inplace_units=True):
    """
    Calculates characteristic lengths for data

    The names of the variables are retrived out the dictionary. You can update the dictionary
    and change the names by using the notation_defs keyworkd, which is a notation object

    Parameters
    ----------
    data: pandas.DataFrame
        dataset to be used. It must either be the raw and turbulent data, or the covariances of such data
    siteConf: pymicra.siteConfig object
        has the site configurations to calculate the obukhovLen
    units: dict
        dict units for the input data
    output_as_df: boolean
        True if you want the output to be a one-line pandas.DataFrame. A pd.Series
        will be output if False.
    inplace_units: bool
        whether or not to update the units dict in place

    Returns
    -------
    pandas.Series or pandas.Dataframe
        depending on return_as_df
    """
    from .. import constants
    from .. import algs
    from ..  import ureg
    import pandas as pd
    import numpy as np

    defs = algs.get_notation(notation)
    defsdic = defs.__dict__
    data = data.copy()
    outunits = {}
    cunits = constants.units

    print('Beginning to extract turbulent scales...')

    #---------
    # First we define the names of the columns according to notation
    u_fluc          =   defs.u_fluctuations
    w_fluc          =   defs.w_fluctuations
    mrho_h2o_fluc   =   defs.h2o_molar_density_fluctuations
    rho_h2o_fluc    =   defs.h2o_mass_density_fluctuations
    theta_fluc      =   defs.thermodyn_temp_fluctuations
    theta_v_fluc    =   defs.virtual_temp_fluctuations
    q_fluc          =   defs.specific_humidity_fluctuations
    solutesf        = [ defsdic['%s_molar_density_fluctuations' % solute] for solute in solutes ]
    solutestars     = [ defsdic['%s_molar_density_star' % solute] for solute in solutes ]
    concsolutesf    = [ defsdic['%s_mass_concentration_fluctuations' % solute] for solute in solutes ]
    #---------

    #---------
    # If data is already covariances we go from there
    if (data.shape[0] == data.shape[1]) and all(data.index == data.columns):
        print('Data seems to be covariances. Will it use as covariances ...')
        cov = data.copy()
        outname = None
    #---------

    #---------
    # If data is raw data, calculate covariances.
    else:
        print('Data seems to be raw data. Will calculate covariances ...')

        #---------
        # Now we try to calculate or identify the fluctuations of theta
        outname = data.index[0]
        theta_mean = data[ defs.thermodyn_temp ].mean()
        if (theta_fluc not in data.columns) or theta_fluct_from_theta_v:
            print('Fluctuations of theta not found. Will try to calculate it ... ', end='')
            #---------
            # We need the mean of the specific humidity and temperature
            if not (units[ theta_v_fluc ]==ureg('kelvin') and units[ defs.thermodyn_temp ]==ureg('kelvin')):
                raise TypeError('\nUnits for both the virtual temperature fluctuations and the thermodynamic temperature fluctuations must be Kelvin')
            data_q_mean =   data[ defs.specific_humidity ].mean()
            data[ theta_fluc ] = (data[theta_v_fluc] - 0.61*theta_mean*data[q_fluc])/(1. + 0.61*data_q_mean)
            theta_fluc_unit = units[ theta_v_fluc ]
            print('done!')
            #---------
        #---------
    
        #-----------
        # First we construct the covariance matrix (slower but more readable than doing it separately)
        # maybe figure out later a way that is both faster and more readable
        print('Calculating the covariances ... ', end='')
        cov = data[[u_fluc, w_fluc, theta_v_fluc, theta_fluc, q_fluc, mrho_h2o_fluc] + solutesf ].cov()
        print('done!')
        #-----------
    #---------


    #---------
    # Now to calculate the characteristic lengths, scales and etc
    print('Calculating the turbulent scales of wind, temperature and humidity ... ', end='')
    out = pd.Series(name=outname)

    u_star  = np.sqrt(-cov.loc[u_fluc, w_fluc])
    out[ defs.u_star ]  = u_star

    theta_v_star = cov.loc[theta_v_fluc, w_fluc] / u_star
    out[ defs.virtual_temp_star ]   = theta_v_star

    out[ defs.thermodyn_temp_star ] = cov.loc[theta_fluc, w_fluc] / u_star
    out[ defs.h2o_molar_density_star ] = cov.loc[ mrho_h2o_fluc, w_fluc ] / u_star
    print('done!')
    #---------

    #---------
    # Now we set the units of the legths
    outunits = {}
    outunits[ defs.u_star ]  = units[ u_fluc ]
    outunits[ defs.virtual_temp_star ]   = units[ theta_v_fluc ]
    outunits[ defs.thermodyn_temp_star ] = units[ theta_v_fluc ]
    outunits[ defs.h2o_molar_density_star ] = units[ mrho_h2o_fluc ]
    #---------

    #---------
    # The solutes have to be calculated separately
    for sol_star, sol_fluc, sol in zip(solutestars, solutesf, solutes):
        print('Calculating the turbulent scale of %s ... ' % sol, end='')
        out[ sol_star ] = cov.loc[sol_fluc, w_fluc] / u_star
        outunits[ sol_star  ] = units[ sol_fluc ]
        print('done!')
    #---------

    #---------
    # We check for the mean virtual temperature
    if not theta_v_mean:
        if defs.mean_virtual_temperature in data.columns:
            theta_v_mean = data[ defs.mean_virtual_temperature ].mean()
            theta_v_mean_unit = units[defs.mean_virtual_temperature]
        else:
            theta_v_mean = data[ defs.virtual_temperature ].mean()
            theta_v_mean_unit = units[defs.virtual_temperature]
    #---------

    #---------
    # Now we calculate the obukhov length and the similarity variable
    print('Calculating Obukhov length and stability parameter ... ', end='')
    Lo = obukhovLen(out, outunits, theta_v_mean=theta_v_mean, theta_v_mean_unit=theta_v_mean_unit, inplace_units=True)
    out[ defs.obukhov_length ]      = Lo
    out[ defs.stability_parameter ] = stabilityParam(Lo, siteConf)

    outunits[ defs.obukhov_length ] = (outunits[ defs.u_star ]**2.)/cunits[ 'gravity' ]
    outunits[ defs.stability_parameter ] = ureg.meter/outunits[ defs.obukhov_length ]
    print('done!')
    #---------

    #---------
    # Create a one-row dataframe if output_as_df is True
    if output_as_df:
        out = out.to_frame().T
    #---------

    #---------
    # Finally we construct the output dataframe
    if inplace_units:
        units.update(outunits)
        return out
    else:
        return out, outunits
    #---------



