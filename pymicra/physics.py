#!/usr/bin/python
"""
Author: Tomas Chor

Module that contains physical functions for general use

TO DO LIST:

* ADD GENERAL SOLAR ZENITH CALCULATION
* ADD PINT FUNCTIONALITY
* ADD FOOTPRINT CALCULATION


"""
from constants import *
from datetime import timedelta
from math import pi,sin,cos,acos,radians,degrees

def gradients(data, levels, order='Crescent'):
    """
    Calculates the gradients for data considering the levels provided.
    UNDER DEVELOPMENT
    
    Parameters
    ----------

    data: pandas DataFrame/ timeSeries
    the data for which the gradients should be calculated

    levels: list
    the columns considered to calculate the gradients

    """
    from algs import combine
    flux=pd.DataFrame(index=data.index)
    for pair in combine(levels, order=order):
        a,b=pair
        grad=data[a]-data[b]
        flux[str(a)+'-'+str(b)]= rho_mean* C_f* u_mean* grad
    return flux


def solarZenith(date, lat=-3.1300, lon=-60.016667, lon0 = -63., negative=False, dr=173):
    """
    Calculates the solar zenith angle at any given day

    needs validation and needs to work without dr or lon0

    Parameters
    ---------
    dr is the julian day of the solstice. Default is June 21
    """
    from calendar import isleap
    if isleap(date.year):   # checks if its a leapyear
        da=366
    else:
        da=365
    #da=lenYear(date.year)   # Number of days in the year
    tt=date.timetuple()
    j_day=tt.tm_yday
    gamma = 2.*pi/365.*(j_day - 1.)
    et = 1./60.* (229.18*(0.000075 + 0.001868*cos(gamma) - 0.032077*sin(gamma)-0.014615*cos(2.*gamma) - 0.04089*sin(2.*gamma)))
    y = (lon-lon0)/15.*60. # Calculo em minutos
    soma = y + 60.*et # Calculo em minutos
    solar_date = date + timedelta(soma/(24.*60.)) # adiciono os minutos ao dia oficial!
    delta = 0.409*cos((2.*pi*(j_day-dr)/da))
    Hsv = solar_date.hour + solar_date.minute/60. + solar_date.second/3600.
    h = 2.*pi/24.*(Hsv - 12.) # Descobrir em que unidade entrar essa hora solar aqui!!
    zen_ang = degrees(acos(sin(radians(lat))*sin(delta) + cos(radians(lat))*cos(delta)*cos(h)))
    if negative:
        if solar_date.hour<=12:
            return zen_ang
        else:
            return -zen_ang
    else:
        return zen_ang




def get_pressure_with_elevation(h, 
  Ps=standard_pressure, Ts=standard_temperature, 
  Tl=temperature_lapse_rate, Hb=0.0, R=R_spec['dry_air'], 
  g=gravity, M=earth_atmosphere_molar_mass):
    """
    This function returns an estimate of the pressure in pascals as a function of
    elevation above sea level
    NOTES:" 
      * This equation is only accurate up to 11000 meters
      * results might be odd for elevations below 0 (sea level), like Dead Sea.
    h=elevation relative to sea level (m)
    Ps= static pressure (pascals)
    Ts= temperature (kelvin)
    Tl= temperature lapse rate (kelvin/meter)
    Hb= height at the bottom of the layer
    R= universal gas constant for air
    g= gravitational acceleration
    M= Molar mass of atmosphere
    P = Ps * (Ts / ((Ts + Tl) * (h - Hb))) ^ ((g * M)/(R * Tl))
    returns pressure in pascals
    """
    if h > 11000.0 :
          print "Warning: Elevation used exceeds the recommended maximum elevation for this function (11,000m)"
    return  Ps * (Ts / (Ts + Tl * (h - Hb))) ** ((g * M) / (R * Tl))



def get_temperature_with_elevation(h, Ts=standard_temperature, Tl=temperature_lapse_rate):
    """This function returns an estimate of temperature as a function above sea level.
    NOTES:
      * This equation is only accurate up to 11,000 meters
      * results might be odd for elevations below 0 (sea level), like Dead Sea.
    Ts= temperature (kelvin)
    Tl= temperature lapse rate (kelvin/meter)
    returns temp in kelvin
    """
    return Ts + h *Tl


def perfGas(p=None, rho=None, R=None, T=None, gas='dry_air'):
    """
    Returns the only value that is not provided in the ideal gas law
    """
    aux=constants()
    R_spec=aux.R_spec[gas]
    if p == None:
        return rho*Ri_spec*T
    elif rho == None:
        return p / (R_spec*T)
    elif T == None:
        return p / (R_spec*rho)
    elif R == None:
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


