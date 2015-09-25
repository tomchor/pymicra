#!/usr/bin/python
"""
Author: Tomas Chor

Defines some useful constants for physical calculations
"""
units={}    # creates units dictionary so it can be updated everytime a constant is added

#---------------------------------------
# Gas constants
#---------------------------------------
molar_mass={'dry_air' : 28.9645,
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
units.update({'R_spec' : 'J/(g * K)'})

mu=R_spec['dry_air']/R_spec['h2o']
units.update({'mu':'1'})

cp_dry=1.0035  # specific heat of dry air at constant pressure
units.update({'cp_dry' : 'J/(g * K)'})

cp_h2o=4.1813  # specific heat of water vapor at constant pressure
units.update({'cp_water' : 'J/(g * K)'})

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
kappa=0.4


#--------------------------------------
# CLEAN DUMMY VARIABLES
#--------------------------------------
del val, key
