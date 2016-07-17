from __future__ import print_function

def MonObuVar(L_m, siteConf):
    """
    Redirects to stabilityParam
    """
    return stabilityParam(L_m, siteConf)

def stabilityParam(L_m, siteConf):
    """
    Calculates the Monin-Obukhov Similarity Variable
    defined as

    zeta = (z-d)/Lm
    where d is the displacement height or zero-plane displacement
    and L_m is the Monin-Obukhov Length.
    """

    z=siteConf.measurement_height
    d=siteConf.displacement_height
    zeta = (z-d)/L_m
    return zeta

def MonObuLen(theta_v_star, theta_v_mean, u_star, g=None):
    """
    """
    return obukhovLen(theta_v_star, theta_v_mean, u_star, g=g)

def obukhovLen(theta_v_star, theta_v_mean, u_star, units=None):
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

    g = constants.gravity
    kappa = constants.kappa
    Lm = - ( (u_star**2) * theta_v_mean) / (kappa *g* theta_v_star)
    return Lm



def turbulentScales(data, siteConf, units, notation=None, theta_v_mean=None,
        theta_fluct_from_theta_v=True, solutes=[], output_as_df=True, inplace=True):
    """
    Calculates characteristic lengths for data

    The names of the variables are retrived out the dictionary. You can update the dictionary
    and change the names by using the notation_defs keyworkd, which is a notation object

    Parameters:
    -----------
    data: pandas.DataFrame
        dataset to be used. It must either be the raw and turbulent data, or the covariances of such data
    siteConf: pymicra.siteConfig object
        has the site configurations to calculate the obukhovLen
    units: dict
        dict units for the input data
    output_as_df: boolean
        True if you want the output to be a one-line pandas.DataFrame. A pd.Series
        will be output if False.

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
    solutesf        = [ defs.molar_density % defs.fluctuations % solute for solute in solutes ]
    solutestars     = [ defsdic[ '%s_molar_density_star' % solute ] for solute in solutes ]
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
            if not (units[ theta_v_fluc ]==ureg['kelvin'] and units[ defs.thermodyn_temp ]==ureg['kelvin']):
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

    out[ defs.specific_humidity_star ] = cov.loc[ q_fluc, w_fluc ] / u_star
    print('done!')
    #---------

    #---------
    # Now we set the units of the legths
    outunits = {}
    outunits[ defs.u_star ]  = units[ u_fluc ]
    outunits[ defs.virtual_temp_star ]   = units[ theta_v_fluc ]
    outunits[ defs.thermodyn_temp_star ] = units[ theta_v_fluc ]
    outunits[ defs.h2o_molar_density_star ] = units[ mrho_h2o_fluc ]
    outunits[ defs.specific_humidity_star ] = units[ q_fluc ]
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
        else:
            theta_v_mean = data[ defs.virtual_temperature ].mean()
    #---------

    #---------
    # Now we calculate the obukhov length and the similarity variable
    print('Calculating Obukhov length and stability parameter ... ', end='')
    Lm = obukhovLen(theta_v_star, theta_v_mean, u_star)
    out[ defs.obukhov_length ]      = Lm
    out[ defs.stability_parameter ] = stabilityParam(Lm, siteConf)

    outunits[ defs.obukhov_length ] = (outunits[ defs.u_star ]**2.)/cunits[ 'gravity' ]
    outunits[ defs.stability_parameter ] = ureg['meter']/outunits[ defs.obukhov_length ]
    print('done!')
    #---------

    #---------
    # Create a one-row dataframe if output_as_df is True
    if output_as_df:
        out = out.to_frame().T
    #---------

    #---------
    # Finally we construct the output dataframe
    print()
    if inplace:
        units.update(outunits)
        return out
    else:
        return out, outunits
    #---------



def get_scales(dataframe, siteConf, notation_defs=None,
  output_as_df=True, theta_from_theta_v=True, solutes=[]):
    """
    Calculates characteristic lengths for data

    The names of the variables are retrived out the dictionary. You can update the dictionary
    and change the names by using the notation_defs keyworkd, which is a notation object

    Parameters:
    -----------
    data: pandas.DataFrame
        dataset to be used. Must have a minimum of columns in order for it to work
    siteConf: pymicra.siteConfig object
        currently not used
    updt: dictionary
        use this dictionary to change a small amount of the names that the variables have in the dataframe
    output_as_df: boolean
        True if you want the output to be a one-line pandas.DataFrame

    CHECKLIST:
    ADD MIXED-LAYER CONVECTION SCALES FOR VELOCITY (w*) AND TEMPERATURE (t*) MAYBE
    """
    from .. import algs
    from .. import constants
    from .. import physics
    import pandas as pd
    import numpy as np

    #if notation_defs==None:
    #    defs=notation()
    #else:
    #    defs=notation_defs
    defs = algs.get_notation(notation)
    data = dataframe.copy()
    
    #---------
    # First we define the names of the columns according to notation
    fluc = defs.fluctuation
    u       =   fluc % defs.u
    v       =   fluc % defs.v
    w       =   fluc % defs.w
    p       =   defs.pressure
    theta   =   defs.thermodyn_temp
    theta_v =   defs.virtual_temp
    theta_v_fluc= fluc % defs.virtual_temp
    q           = defs.specific_humidity
    qfluct      = fluc % defs.specific_humidity
    solutesf    = [ fluc % el for el in solutes ]
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
    theta_v_std     = data[theta_v_fluc].std()

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
    if theta_from_theta_v:
        theta_star  = (theta_v_star - 0.61*theta_mean*q_star)/(1.+0.61*q_mean)
        theta_std   = physics.theta_std_from_theta_v(data[theta_v_fluc], data[qfluct], theta_v_mean, q_mean, theta_mean)
    else:
        theta_star  = cov.loc[theta_fluc, w] / u_star
        theta_std   = data[theta].std()
    #---------

    #---------
    # Now we calculate the obukhov length and the similarity variable
    Lm = obukhovLen(theta_v_star, theta_v_mean, u_star, g=constants.gravity)
    zeta = stabilityParam(Lm, siteConf)
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


