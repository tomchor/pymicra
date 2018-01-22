"""
Module that contains physical functions. They are all general use, but
most are specially frequent in micrometeorology.

TO DO LIST:
 * ADD GENERAL SOLAR ZENITH CALCULATION
 * ADD FOOTPRINT CALCULATION?
"""
from __future__ import absolute_import, print_function, division

def specific_humidity_from_ppxv(data, units, notation=None, return_full_df=True, inplace_units=True):
    """Calculates the specific humidity q from values of molar concentration of
    water (ppmv, ppthv and etc).

    The equation is
                                      mv x
                             q = ----------------
                                 (mv - md) x + md
    where x is the molar concentration in ppxv.

    Parameters
    ----------
    data: pandas.dataframe
        dataset
    units: dict
        units dictionary
    notation: pymicra.Notation
        notation to be used
    return_full_df: bool
        whether to return only the calculated series or the full df
    inplace_units: bool
        whether to return only a dict with the units of the new variables or include them in "units"

    Returns
    -------
    outdata: pandas.Series, pandas.DataFrame
        specific humidity
    """

    q = (mv*x)/(md+(mv-md)*x)
    return

def theta_s_from_c(data, units, notation=None, return_full_df=True, inplace_units=True):
    r"""Calculates sonic temperature using speed of sound

    From Schotanus, Nieuwstadt, de Bruin; DOI 10.1007/BF00164332

    theta_s = 1/403 * c**2

    :math:`theta_s \approx 1/403 c^2`

    Parameters
    ----------
    data: pandas.dataframe
        dataset
    units: dict
        units dictionary
    notation: pymicra.Notation

    Returns
    -------
    pandas.Series
        sonic temperature
    """
    from . import algs

    defs = algs.get_notation(notation)
    data = data.copy()
    
    theta_s = 1/403 * data[ defs.sound_speed ]**2
    theta_s_unit = units['K']

    if return_full_df:
        data.loc[:, defs.sonic_temp ] = theta_s
        out = data
    else:
        out = theta

    if inplace_units:
        units.update({ defs.sonic_temp : theta_s_unit })
    else:
        out = (out, theta_s_unit)

    return out

def theta_from_theta_s(data, units, notation=None, return_full_df=True, inplace_units=True):
    r"""Calculates thermodynamic temperature using sonic temperature measurements

    From Schotanus, Nieuwstadt, de Bruin; DOI 10.1007/BF00164332

    theta_s = theta (1 + 0.51 q) (1 - (vn/c)**2)^0.5

    :math:`theta_s \approx theta (1 + 0.51 q)`

    Parameters
    ----------
    data: pandas.dataframe
        dataset
    units: dict
        units dictionary
    notation: pymicra.Notation

    Returns
    -------
    pandas.Series
        thermodynamic temperature
    """
    from . import algs

    defs = algs.get_notation(notation)
    data = data.copy()
    theta = data[ defs.sonic_temp ]/(1. + 0.51*data[ defs.specific_humidity ])
    theta_unit = units[ defs.sonic_temp ]

    if return_full_df:
        data.loc[:, defs.thermodyn_temp ] = theta
        out = data
    else:
        out = theta

    if inplace_units:
        units.update({ defs.thermodyn_temp : theta_unit })
    else:
        out = (out, theta_unit)

    return out


def theta_from_theta_v(data, units, notation=None, return_full_df=True, inplace_units=True):
    r"""Calculates thermodynamic temperature from virtual temperature measurements

    :math:`theta_v \approx theta (1 + 0.61 q)`

    Parameters
    ----------
    data: pandas.dataframe
        dataset
    units: dict
        units dictionary
    notation: pymicra.Notation

    Returns
    -------
    pandas.DataFrame or Series
        virtual temperature
    """
    from . import algs

    defs = algs.get_notation(notation)
    data = data.copy()
    theta = data[ defs.virtual_temp ]/(1. + 0.61*data[ defs.specific_humidity ])
    theta_unit = units[ defs.virtual_temp ]

    if return_full_df:
        data.loc[:, defs.thermodyn_temp ] = theta
        out = data
    else:
        out = theta

    if inplace_units:
        units.update({ defs.thermodyn_temp : theta_unit })
    else:
        out = (out, theta_unit)

    return out


def theta_fluc_from_theta_v_fluc(data, units, notation=None, return_full_df=True, inplace_units=True):
    """
    Derived from theta_v = theta(1 + 0.61 q)

    Parameters
    ----------
    data: pandas.dataframe
        dataframe with q, q', theta, theta_v'
    units: dict
        units dictionary
    notation: pymicra.Notation
        Notation object or None

    Returns
    -------
    float
        standard deviation of the thermodynamic temperature
    """
    import numpy as np
    from . import algs

    defs = algs.get_notation(notation)
    data = data.copy()

    #-----------
    # Consider the existence of mean q
    if defs.mean_specific_humidity in data.columns:
        q_mean = data[ defs.mean_specific_humidity ]
    else:
        q_mean = data[ defs.specific_humidity ].mean()
    #-----------

    #-----------
    # Consider existence of mean theta
    if defs.mean_thermodyn_temp in data.columns:
        theta_mean = data[ defs.mean_thermodyn_temp ]
    else:
        theta_mean = data[ defs.thermodyn_temp ].mean()
    #-----------

    theta_fluc = (data[ defs.virtual_temperature_fluctuations ] - 0.61*theta_mean*data[ defs.specific_humidity_fluctuations ])/(1.+0.61*q_mean)
    theta_fluc_unit = units[ defs.virtual_temperature_fluctuations ]

    if return_full_df:
        data.loc[:, defs.thermodyn_temp_fluctuations ] = theta_fluc
        out = data
    else:
        theta_fluc.name = defs.thermodyn_temp_fluctuations
        out = theta_fluc

    if inplace_units:
        units.update({ defs.thermodyn_temp_fluctuations : theta_fluc_unit })
    else:
        out = (out, theta_fluc_unit)

    return out


def theta_std_from_theta_v_fluc(data, units, notation=None):
    """
    Derived from theta_v = theta(1 + 0.61 q)

    Parameters
    ----------
    data: pandas.dataframe
        dataframe with q, q', theta, theta_v'
    units: dict
        units dictionary
    notation: pymicra.Notation
        Notation object or None

    Returns
    -------
    float
        standard deviation of the thermodynamic temperature
    """
    import numpy as np
    from . import algs

    defs = algs.get_notation(notation)
    data = data.copy()
    units = units.copy()

    #-----------
    # Consider the existence of mean q
    if defs.mean_specific_humidity in data.columns:
        q_mean = data[ defs.mean_specific_humidity ]
    else:
        q_mean = data[ defs.specific_humidity ].mean()
    #-----------

    #-----------
    # Consider existence of mean theta
    if defs.mean_thermodyn_temp in data.columns:
        theta_mean = data[ defs.mean_thermodyn_temp ]
    else:
        theta_mean = data[ defs.thermodyn_temp ].mean()
    #-----------

    denom = (1. + 0.61*q_mean)**2.

    num1 = (data[ defs.virtual_temp_fluctuations ]**2.).mean()
    num2 = 2. * 0.61 * theta_mean * (data[ defs.virtual_temp_fluctuations] * data[ defs.specific_humidity_fluctuations ]).mean()
    num3 = ((0.61*theta_mean)**2.) * ( data[ defs.specific_humidity_fluctuations ]**2.).mean()

    var = (num1 - num2 + num3)/denom
    theta_fluc = np.sqrt(var)

    return theta_fluc



def ppxv2density(data, units, notation=None, inplace_units=True, solutes=[]):
    r"""
    Calculates density of solutes based on their molar concentration (ppmv, ppbv and etc), not
    to be confused with mass concentration (ppm, ppb and etc).

    Uses the relation
    :math:`\rho_x = \frac{C p}{\theta R_x}`

    Parameters
    ----------
    data: pandas.DataFrame
        dataset of micromet variables
    units: dict
        dict of pint units
    notation: pymicra.Notation
        notation to be used here
    inplace_units: bool
        whether or not to treat the dict units in place
    solutes: list or tuple
        solutes to consider when doing this conversion

    Returns
    -------
    pandas.DataFrame
        input data plus calculated density columns
    """
    from . import constants
    from . import algs

    if not solutes:
        raise ValueError('Should specify for which solutes')

    defs = algs.get_notation(notation)
    defsdic = vars(defs)
    data = data.copy()
    cunits = constants.units
    sol_units = {}
    convert_to={}

    p = data[ defs.pressure ]
    theta = data[ defs.thermodyn_temp ]

    for solute in solutes:
        sol_mconc = defsdic[ '%s_molar_concentration' % solute ]
        rho_sol = defsdic[ '%s_mass_density' % solute ]
        data[ rho_sol ] = data[ sol_mconc ]*p/(theta*constants.R_spec[solute])
        unit = units[ sol_mconc ]*units[ defs.pressure ]/(units[ defs.thermodyn_temp ]*cunits['R_spec'])

        sol_units.update({ rho_sol : unit })
        convert_to.update({ rho_sol : 'kg/m**3' })
    data = data.convert_cols(convert_to, sol_units, inplace_units=True)

    if inplace_units:
        units.update(sol_units)
        return data
    else:
        return data, unit


def airDensity_from_theta_v(data, units, notation=None, inplace_units=True, use_means=False, return_full_df=True):
    """
    Calculates moist air density using p = rho R_dry T_virtual

    Parameters
    ----------
    data: pandas.DataFrame
        data to use to calculate air density
    units: dict
        dictionary of units
    notation: pymicra.Notation
        notation to be used
    inplace_units: bool
        whether or not to update the units inplace. If False, units are returns too
    use_means: bool
        whether or not to use averages of pressure and virtual temperature, instead of the means plus fluctuations

    Returns
    -------
    
    """
    from . import algs
    from . import constants

    defs = algs.get_notation(notation)
    data = data.copy()
    if not inplace_units:
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
    data = data.convert_cols({defs.moist_air_mass_density:'kg/m**3'}, units, inplace_units=True)
    #---------

    if inplace_units:
        return data
    else:
        return data, units


def dryAirDensity_from_p(data, units, notation=None, inplace_units=True):
    """
    Calculates dry air density
    NEEDS IMPROVEMENT REGARDING HANDLING OF UNITS
    """
    from . import algs
    from . import constants

    defs = algs.get_notation(notation)
    data = data.copy()
    if not inplace_units:
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

    if inplace_units:
        return data
    else:
        return data, units


def airDensity_from_theta(data, units, notation=None, inplace_units=True, use_means=False, theta=None, theta_unit=None):
    """
    Calculates moist air density using theta measurements

    Parameters
    -----------
    data: pandas.DataFrame
        dataset to add rho_air
    units: dict
        units dictionary
    notation: pymicra.notation
        notation to be used
    inplace_units: bool
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
    #-----------

    #-----------
    # We subtract p_h2o from p_dry using pymicra.algs.add because addition with different units is tricky
    p_dry, p_dry_unit = algs.add([ data[defs.pressure], -p_h2o ], [ units[defs.pressure], p_h2o_unit ], inplace_units=False)
    #-----------

    #-----------
    # We calculate dry air mass density
    data[ defs.dry_air_mass_density ] = p_dry/(R_dry*theta)
    outunits[ defs.dry_air_mass_density ] = p_dry_unit/(R_spec_unit*theta_unit)
    data = data.convert_cols({defs.dry_air_mass_density : 'kg/m**3'}, outunits, inplace_units=True)
    #-----------

    #-----------
    # We add rho_dry to rho_h2o using pymicra.algs.add because addition with different units is tricky
    data.loc[:, defs.moist_air_mass_density ] = algs.add([data[defs.dry_air_mass_density], data[defs.h2o_mass_density]],
                        [outunits[defs.dry_air_mass_density], units[defs.h2o_mass_density]], inplace_units=True, unitdict=outunits, key=defs.moist_air_mass_density)
    #-----------

    #-----------
    # Adjust the units
    conversions = { defs.moist_air_mass_density : 'kg/m**3',
                    defs.dry_air_mass_density : 'kg/m**3'}
    data = data.convert_cols(conversions, outunits, inplace_units=True)
    #-----------

    if inplace_units:
        units.update(outunits)
        return data
    else:
        return data, outunits


def latent_heat_water(T):
    """
    Calculates the latent heat of evaporation for water

    Receives T in Kelvin and returns the latent heat in J/g
    """
    return 2500.827 -2.360*(T-273.15)


def satWaterPressure(T, unit='kelvin'):
    """
    Returns the saturated water vapor pressure according eq (3.97) of Wallace and Hobbes, page 99.

    e0, b, T1 and T2 are constants specific for water vapor

    Parameters
    ----------
    T: float
        thermodynamic temperature

    Returns
    -------
        saturated vapor pressure of water (in kPa)
    """
    from np import exp
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
    from .constants import R_spec
    
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

    Parameters
    ----------
    q: float
        the specific humidity in g(water)/g(air)

    Returns
    -------
    R_air: float
        the specific gas constant for humid air in J/(g*K)
    """
    from .constants import R_spec

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


def _solarZenith(date, lat=-3.1300, lon=-60.016667, lon0 = -63., negative=False, dr=None):
    """
    Calculates the solar zenith angle at any given day

    needs validation and needs to work without lon0

    Parameters
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

    Returns
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
        from .constants import sumsolstice
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


