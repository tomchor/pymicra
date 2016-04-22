#!/usr/bin/python
"""
Author: Tomas Chor

Module that contains physical functions for general use

TO DO LIST:

* ADD GENERAL SOLAR ZENITH CALCULATION
* MAYBE ADD PINT FUNCTIONALITY
* ADD FOOTPRINT CALCULATION

"""
from constants import *
from datetime import timedelta


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
    if gas==None and R==None:
        R=R_spec[gas]

    if type(p) == type(None):
        return rho*R*T
    elif type(rho) == type(None):
        return p / (R*T)
    elif type(T) == type(None):
        return p / (R*rho)
    elif type(R) == type(None):
        return p / (rho*T)
    return


def wetAirDens(p=None, T=None, q=None):
    """
    From R. S. Davis, Equation for the Determination of the Density of Moist Air (1981/91).
    Available at http://www.nist.gov/calibrations/upload/metv29i1p67-2.pdf
    """
    aux=constants()
    R_spec=aux.R_spec
    R_dry=R_spec['dry_air']
    R_h2o=R_spec['h2o']
    rho_wet=(p/(R_dry*T)) * (1. - q*(1. - R_dry/R_h2o))
    return rho_wet


def virtualTemp(T, tpres, ppres):
    """
    Gets virtual temperature from thermodynamic temperature, total pressure and water vapor pressure
    mu=R_dry/R_water_vapor is approx 0.622
    """
    virt_temp=temp/(1. -(ppres/tpres)*(1.-mu))
    return virt_temp


def R_umidAir(q):
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
    return q* R_spec['h2o'] + (1.0 - q)*R_spec['dry_air']


def dewPointTemp(theta, e):
    """
    Calculates the dew point temperature.
    theta has to be in Kelvin and e in kPa
    """
    ln = np.log(e / 0.611)
    coef = ln / (17.502 - ln)
    return 240.97*coef + 273.16

