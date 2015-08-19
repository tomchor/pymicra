#!/usr/bin/python
"""
Author: Tomas Chor

Module that contains physical functions for general use
"""

molar_mass={'dry_air' : 28.9645,
	'o3'  : 47.99820,
	'h2o' : 18.0153,
	'co2' : 44.0095,
	'co'  : 28.0101,
	'ch4' : 16.0425,
	'n2o' : 44.01280,
	'o' : 15.99940,
	'n' : 14.00670}

R=8.3144621     #universal gas constant J/(mol.K)
R_spec={}
for key, val in molar_mass.iteritems():
	R_spec.update( {key : R/val} )


virt_temp=temp/(1. -(e/p)*(1.-eps))

