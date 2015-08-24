#!/usr/bin/python
"""
Author: Tomas Chor

Defines some useful constants for physical calculations
"""
# Gas constants
molar_mass={'dry_air' : 28.9645,
            'o3'  : 47.99820,
            'h2o' : 18.0153,
            'co2' : 44.0095,
            'co'  : 28.0101,
            'ch4' : 16.0425,
            'n2o' : 44.01280,
            'o' : 15.99940,
            'n' : 14.00670}
    
R=8.3144621    # universal gas constant J/(mol.K)
R_spec={}
for key, val in molar_mass.iteritems():
    R_spec.update( {key : R/val} )
mu=R_spec['dry_air']/R_spec['h2o']
cp_dry=1003.5  # specific heat of dry air at constant pressure
cp_h2o=1850.0  # specific heat of water vapor at constant pressure

# physical constants
gravity=9.80665      # gravity
omega = 7.29212E-5  # angular valocity of the earth
earth_radius=6378140.0 # meters
standard_pressure = 101325.00 # pascals
standard_temperature = 288.15 # kelvin
celsius_offset = 273.15 # subtract from kelvin to get deg C, add to deg C to get kelvin
temperature_lapse_rate = -0.0065 # change in temperature with height, kelvin/metre
earth_atmosphere_molar_mass = 0.0289644 # kg/mol


# units
units={'molar_mass' : 'g/mol',
            'R' : 'J/(mol * K)',
            'R_spec' : 'J/(g * K)'}


