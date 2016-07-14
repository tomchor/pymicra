#!/usr/bin/python
"""
Author: Tomas Chor

Module that contains physical functions for general use

TO DO LIST:

* ADD GENERAL SOLAR ZENITH CALCULATION
* ADD FOOTPRINT CALCULATION?
"""


def theta_from_theta_s(data, units, notation=None, return_df=True):
    """
    From Schotanus, Nieuwstadt, de Bruin; DOI 10.1007/BF00164332

    theta_s = theta (1 + 0.51 q) (1 - (vn/c)**2)**0.5
    theta_s ~ theta (1 + 0.51 q)

    Parameter:
    ----------
    data: pandas.dataframe
        dataset
    units: dict
        units dictionary
    notation: pymicra.Notation
    """
    from . import algs

    defs = algs.get_notation(notation)
    data = data.copy()
    theta = data[ defs.sonic_temp ]/(1. + 0.51*data[ defs.specific_humidity ])
    if return_df:
        data.loc[:, defs.thermodyn_temp ] = theta
        return data
    else:
        return theta


def theta_from_theta_v(data, units, notation=None, return_df=True):
    """

    theta_v ~ theta (1 + 0.61 q)

    Parameter:
    ----------
    data: pandas.dataframe
        dataset
    units: dict
        units dictionary
    notation: pymicra.Notation
    """
    from . import algs

    defs = algs.get_notation(notation)
    data = data.copy()
    theta = data[ defs.sonic_temp ]/(1. + 0.61*data[ defs.specific_humidity ])
    if return_df:
        data.loc[:, defs.thermodyn_temp ] = theta
        return data
    else:
        return theta



def theta_std_from_theta_v_fluc(data, units, notation=None, return_df=True):
    """
    Derived from theta_v = theta(1 + 0.61 q)

    Parameters:
    -----------
    theta_v: array
    q: array
    theta_v_mean: float
    q_mean: float
    theta_mean: float
    """
    from . import algs

    defs = algs.get_notation(notation)
    data = data.copy()

    denom = (1. + 0.61*data[ defs.mean_specific_humidity ])**2.
    num1 = data[ defs.virtual_temp ]**2.

    if defs.mean_thermodyn_temp in data.columns:
        theta_mean = data[ defs.mean_thermodyn_temp ]
    else:
        theta_mean = data[ defs.thermodyn_temp ].mean()

    num2 = 2.*0.61*theta_mean*np.nanmean(theta_v * q)
    num3 = ((0.61*theta_mean)**2.) * np.nanmean(q*q)
    var = (num1 - num2 + num3)/denom
    return np.sqrt(var)



def theta_std_from_theta_v(theta_v, q, theta_v_mean, q_mean, theta_mean):
    """
    Derived from theta_v = theta(1 + 0.61 q)

    Parameters:
    -----------
    theta_v: array
    q: array
    theta_v_mean: float
    q_mean: float
    theta_mean: float
    """
    import numpy as np
    denom = (1. + 0.61*q_mean)**2.
    num1 = np.nanmean(theta_v*theta_v) 
    num2 = 2.*0.61*theta_mean*np.nanmean(theta_v * q)
    num3 = ((0.61*theta_mean)**2.) * np.nanmean(q*q)
    var = (num1 - num2 + num3)/denom
    return np.sqrt(var)


def ppxv2density(ser, T, p, units, solute=None):
    """
    ser should be a series!
    concentration in ser should be ppmv, which will be transformed to g/m3
    p should be in kPa!
    """
    from . import constants
    from . import algs
    from . import ureg

    units = units.copy()
    ser = ser.copy()

    out = p/(T*constants.R_spec[solute])
    out = ser.multiply(out, axis=0)
    unit = units[solute]*units['p']/(units['T']*constants.units['R_spec'])
    out, unit = algs.convert_to(out, unit, 'kg/m**3')
    units.update({solute:unit})
    return out, units


def airDensity_from_theta_v(data, units, notation=None, inplace=True, use_means=False):
    """
    Calculates moist air density using p = rho R_dry T_virtual
    """
    from . import algs
    from . import constants

    defs = algs.get_notation(notation)
    data = data.copy()
    if not inplace:
        units = units.copy()
    R_dry = constants.R_spec['dry']

    #---------
    # Calculates air density and its unit with theta_v and R_dry
    if use_means:
        data.loc[:, defs.moist_air_mass_density ] = data[ defs.pressure ].mean()/(R_dry * data[ defs.virtual_temp ].mean())
    else:
        data.loc[:, defs.moist_air_mass_density ] = data[ defs.pressure ]/(R_dry * data[ defs.virtual_temp ])
    units.update({ defs.moist_air_mass_density : units[ defs.pressure ]/(constants.units['R_spec']*units[ defs.virtual_temp ]) })
    #---------

    #---------
    # We pass it to the standard mass density unit: kg/m**3
    data = data.convert_cols({defs.moist_air_mass_density:'kg/m**3'}, units, inplace=True)
    #---------

    if inplace:
        return data
    else:
        return data, units

def dryAirDensity_from_p(data, units, notation=None, inplace=True):
    """
    Calculates dry air density
    NEEDS IMPROVEMENT REGARDING HANDLING ON UNITS
    """
    from . import algs
    from . import constants

    defs = algs.get_notation(notation)
    data = data.copy()
    if not inplace:
        units = units.copy()

    Rdry = constants.R_spec['dry']
    Rh2o = constants.R_spec['h2o']
    Runit = constants.units['R_spec']

    p_h2o = data[ defs.h2o_mass_density ]*Rh2o*data[ defs.thermodyn_temp ]
    p_h2o_unit = units[ defs.h2o_mass_density ]*Runit*units[ defs.thermodyn_temp ]

    #-----------
    # Before adding or subtracting we need to make sure both units are the same
    p_h2o, p_h2o_unit = algs.convert_to(p_h2o, p_h2o_unit, 'kPa')
    p_air, p_air_unit = algs.convert_to(data[ defs.pressure ], units[ defs.pressure ], 'kPa')

    p_dry = p_air - p_h2o
    #-----------

    #-----------
    data.loc[:, defs.dry_air_mass_density ] = p_dry/(Rdry * data[ defs.thermodyn_temp ])
    units.update({ defs.dry_air_mass_density : p_air_unit/(Runit * units[ defs.thermodyn_temp ])})
    #-----------

    if inplace:
        return data
    else:
        return data, units


def airDensity_from_theta(data, units, notation=None, inplace=True, use_means=False, theta=None, theta_unit=None):
    """
    Calculates moist air density using theta measurements

    Parameters:
    -----------
    data: pandas.DataFrame
        dataset to add rho_air
    units: dict
        units dictionary
    notation: pymicra.notation
        notation to be used
    inplace: bool
        whether or not to treat units inplace
    use_means: bool
        use the mean of theta or not when calculating
    theta: pandas.Series
        auxiliar theta measurement
    theta_unit: pint.quantity
        auxiliar theta's unit
    """
    from . import algs
    from . import constants

    defs = algs.get_notation(notation)
    data = data.copy()
    outunits = {}

    R_dry = constants.R_spec['dry']
    R_h2o = constants.R_spec['h2o']
    R_spec_unit = constants.units['R_spec']

    if not theta:
        theta = data[ defs.thermodyn_temp ]
        theta_unit = units[ defs.thermodyn_temp ]

    if use_means:
        theta = theta.mean()

    #-----------
    # We calculate the water vapor partial pressure
    p_h2o = data[ defs.h2o_mass_density ] * R_h2o * theta
    p_h2o_unit = units[ defs.h2o_mass_density ]*R_spec_unit*theta_unit
    #pv, pvu = algs.multiply([ data[defs.h2o_mass_density], R_h2o, theta ], [ units[defs.h2o_mass_density], R_spec_unit, theta_unit ])
    #print(p_h2o, p_h2o_unit)
    #-----------

    #-----------
    # We subtract p_h2o from p_dry using pymicra.algs.add because addition with different units is tricky
    p_dry, p_dry_unit = algs.add([ data[defs.pressure], -p_h2o ], [ units[defs.pressure], p_h2o_unit ], inplace=False)
    #-----------

    #-----------
    # We calculate dry air mass density
    data[ defs.dry_air_mass_density ] = p_dry/(R_dry*theta)
    outunits[ defs.dry_air_mass_density ] = p_dry_unit/(R_spec_unit*theta_unit)
    data = data.convert_cols({defs.dry_air_mass_density : 'kg/m**3'}, outunits, inplace=True)
    #-----------

    #-----------
    # We add rho_dry to rho_h2o using pymicra.algs.add because addition with different units is tricky
    data.loc[:, defs.moist_air_mass_density ] = algs.add([data[defs.dry_air_mass_density], data[defs.h2o_mass_density]],
                        [outunits[defs.dry_air_mass_density], units[defs.h2o_mass_density]], inplace=True, unitdict=outunits, key=defs.moist_air_mass_density)
    #-----------

    #-----------
    # Adjust the units
    conversions = { defs.moist_air_mass_density : 'kg/m**3',
                    defs.dry_air_mass_density : 'kg/m**3'}
    data = data.convert_cols(conversions, outunits, inplace=True)
    #-----------

    if inplace:
        units.update(outunits)
        return data
    else:
        return data, outunits


def solarZenith(date, lat=-3.1300, lon=-60.016667, lon0 = -63., negative=False, dr=None):
    """
    Calculates the solar zenith angle at any given day

    needs validation and needs to work without lon0

    Parameters:
    -----------
    date: datetime object
        the date and time for which the zenith angle has to be calculated
    lat: float
        latitude in degrees
    lon: float
        longitude in degrees
    lon0: float
        DEFINE IT BETTER

    dr is the julian day of the solstice. Default is to get from dictionary

    Returns:
    --------
    zen_ang: float
        the zenith angle in degrees
    """
    from math import pi,sin,cos,acos,radians,degrees
    from calendar import isleap
    from datetime import timedelta

    #----------
    # Finds out which is the day of the solstice for that year and puts it into dr variable
    if dr==None:
        from constants import sumsolstice
        dr=sumsolstice[date.year].timetuple().tm_yday
    #----------

    #------------
    # Checks if it's a leapyear and defines the da (days in year)
    if isleap(date.year):
        da=366
    else:
        da=365
    #------------

    tt=date.timetuple()
    j_day=tt.tm_yday
    gamma = 2.*pi/365.*(j_day - 1.)
    et = 1./60.* (229.18*(0.000075 + 0.001868*cos(gamma) - 0.032077*sin(gamma)-0.014615*cos(2.*gamma) - 0.04089*sin(2.*gamma)))
    y = (lon-lon0)/15.*60. # minutes
    soma = y + 60.*et # minutes
    solar_date = date + timedelta(soma/(24.*60.)) # add minutes to the official time
    delta = 0.409*cos((2.*pi*(j_day-dr)/da))
    Hsv = solar_date.hour + solar_date.minute/60. + solar_date.second/3600.
    h = 2.*pi/24.*(Hsv - 12.)
    zen_ang = degrees(acos(sin(radians(lat))*sin(delta) + cos(radians(lat))*cos(delta)*cos(h)))

    #----------
    # If negative==True then the angles before noon will be negative
    if negative:
        if solar_date.hour<=12:
            return zen_ang
        else:
            return -zen_ang
    else:
        return zen_ang
    #----------


def latent_heat_water(T):
    """
    Calculates the latent heat of evaporation for water

    Receives T in Kelvin and returns the latent heat in J/g
    """
    return 2500.827 -2.360*(T-273.15)


def CtoK(T):
    """
    Return temp in Kelvin given temp T in Celsius
    """
    return T + 273.15


def satWaterPressure(T, unit='kelvin'):
    """
    Returns the saturated water vapor pressure according eq (3.97) of Wallace and Hobbes, page 99.

    e0, b, T1 and T2 are constants specific for water vapor

    Parameters:
    -----------

    T: float
        thermodynamic temperature

    Returns:
    --------
        saturated vapor pressure of water (in kPa)
    """
    from math import exp
    e0=0.61094
    b=17.2694
    if unit=='kelvin':
        T1=273.16
        T2=35.86
    elif units=='celsius':
        T1=0.
        T2=243.04
    else:
        raise TypeError('Check your units')
    brackets=b*(T-T1)/(T-T2)
    return e0*exp(brackets)


def perfGas(p=None, rho=None, R=None, T=None, gas=None):
    """
    Returns the only value that is not provided in the ideal gas law

    P.S.: I'm using type to identify None objects because this way it works
    againt pandas objects
    """
    from constants import R_spec
    
    if R==None:
        if gas != None:
            R=R_spec[gas]
        elif gas==None:
            R=R_spec['dry']

    if type(p) == type(None):
        return rho*R*T
    elif type(rho) == type(None):
        return p / (R*T)
    elif type(T) == type(None):
        return p / (R*rho)
    elif type(R) == type(None):
        return p / (rho*T)
    return



def R_moistAir(q):
    """
    Calculates the gas constant for umid air from the specific humidity q

    Parameters:
    -----------
    q: float
        the specific humidity in g(water)/g(air)

    Returns:
    --------
    R_air: float
        the specific gas constant for humid air in J/(g*K)
    """
    from constants import R_spec

    return q* R_spec['h2o'] + (1.0 - q)*R_spec['dry_air']


def dewPointTemp(theta, e):
    """
    Calculates the dew point temperature.
    theta has to be in Kelvin and e in kPa
    """
    import numpy as np
    ln = np.log(e / 0.611)
    coef = ln / (17.502 - ln)
    return 240.97*coef + 273.16

