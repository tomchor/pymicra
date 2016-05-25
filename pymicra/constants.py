#!/usr/bin/python
"""
Defines some useful constants
"""
import algs
units={}    # creates units dictionary so it can be updated everytime a constant is added

#---------------------------------------
# Gas constants
#---------------------------------------
molar_mass={'dry' : 28.9645,
            'o3'  : 47.99820,
            'h2o' : 18.0153,
            'co2' : 44.0095,
            'co'  : 28.0101,
            'ch4' : 16.0425,
            'n2o' : 44.01280,
            'o' : 15.99940,
            'n' : 14.00670}
units.update({'molar_mass' : 'g/mol'})
   
R=8.3144621    # universal gas constant J/(mol.K)
units.update({'R' : 'J/(mol * K)'})

R_spec={}
for key, val in molar_mass.iteritems():
    R_spec.update( {key : R/val} )
del key, val
units.update({'R_spec' : 'J/(g * K)'})

units.update({'mu':'1'})

cp_dry=1.0035  # specific heat of dry air at constant pressure
units.update({'cp_dry' : 'J/(g * K)'})

cp_h2o=4.1813  # specific heat of water vapor at constant pressure
units.update({'cp_water' : 'J/(g * K)'})

from physics import latent_heat_water
units.update({'latent_heat_water' : 'J/g'})

from physics import satWaterPressure
units['satWaterPressure'] = 'kPa'


#---------------------------------------
# physical constants
#---------------------------------------
gravity=9.80665      # gravity
units.update({'gravity':'m/(s**2)'})

omega = 7.29212E-5  # angular velocity of the earth
units.update({'omega' : '1/s'})

earth_radius=6378140.0 # meters
units.update({'earth_radius' :'m'})

standard_pressure = 101325.00 # pascals
units.update({'standard_pressure':'Pa'})

standard_temperature = 288.15 # kelvin
units.update({'standard_temperature':'K'})

temperature_lapse_rate = -0.0065 # change in temperature with height, kelvin/metre
units.update({'temperature_lapse_rate' : 'K/m'})

earth_atmosphere_molar_mass =28.9644 # g/mol
units.update({'earth_atmosphere_molar_mass' : 'g/mol' })


#---------------------------------------
# OTHER CONSTANTS
#---------------------------------------
mole = 6.022140857e23

kappa=0.4

greek_alphabet = {
    u'\u0391': 'Alpha',
    u'\u0392': 'Beta',
    u'\u0393': 'Gamma',
    u'\u0394': 'Delta',
    u'\u0395': 'Epsilon',
    u'\u0396': 'Zeta',
    u'\u0397': 'Eta',
    u'\u0398': 'Theta',
    u'\u0399': 'Iota',
    u'\u039A': 'Kappa',
    u'\u039B': 'Lamda',
    u'\u039C': 'Mu',
    u'\u039D': 'Nu',
    u'\u039E': 'Xi',
    u'\u039F': 'Omicron',
    u'\u03A0': 'Pi',
    u'\u03A1': 'Rho',
    u'\u03A3': 'Sigma',
    u'\u03A4': 'Tau',
    u'\u03A5': 'Upsilon',
    u'\u03A6': 'Phi',
    u'\u03A7': 'Chi',
    u'\u03A8': 'Psi',
    u'\u03A9': 'Omega',
    u'\u03B1': 'alpha',
    u'\u03B2': 'beta',
    u'\u03B3': 'gamma',
    u'\u03B4': 'delta',
    u'\u03B5': 'epsilon',
    u'\u03B6': 'zeta',
    u'\u03B7': 'eta',
    u'\u03B8': 'theta',
    u'\u03B9': 'iota',
    u'\u03BA': 'kappa',
    u'\u03BB': 'lamda',
    u'\u03BC': 'mu',
    u'\u03BD': 'nu',
    u'\u03BE': 'xi',
    u'\u03BF': 'omicron',
    u'\u03C0': 'pi',
    u'\u03C1': 'rho',
    u'\u03C3': 'sigma',
    u'\u03C4': 'tau',
    u'\u03C5': 'upsilon',
    u'\u03C6': 'phi',
    u'\u03C7': 'chi',
    u'\u03C8': 'psi',
    u'\u03C9': 'omega',
}

from datetime import datetime
sumsolstice={
2010:'21 11:28',
2011:'21 17:16',
2012:'20 23:09',
2013:'21 05:04',
2014:'21 10:51',
2015:'21 16:38',
2016:'20 22:34',
2017:'21 04:24'
}
sumsolstice={ key : datetime.strptime('{0}-06-{1}'.format(key, val), '%Y-%m-%d %H:%M') for key, val in sumsolstice.iteritems() }

 
#--------------------------------------
# CLEAN DUMMY VARIABLES
#--------------------------------------
try:
    del val, key
except:
    pass

#---------------
# Parse units to pint!
units = algs.parseUnits(units)
#---------------
