#!/usr/bin/python
"""
Author: Tomas Chor

Module that contains physical functions for general use

ADD PINT FUNCTIONALITY LATER
"""

class constants(object):
	def __init__(self):
		self.molar_mass={'dry_air' : 28.9645,
			'o3'  : 47.99820,
			'h2o' : 18.0153,
			'co2' : 44.0095,
			'co'  : 28.0101,
			'ch4' : 16.0425,
			'n2o' : 44.01280,
			'o' : 15.99940,
			'n' : 14.00670}
	
		self.R=8.3144621     #universal gas constant J/(mol.K)
		R_spec={}
		for key, val in self.molar_mass.iteritems():
			R_spec.update( {key : self.R/val} )
		self.R_spec=R_spec
		self.units={'molar_mass' : 'g/mol',
			'R' : 'J/(mol * K)',
			'R_spec' : 'J/(g * K)'}
	
	
temp,e,p,eps=[1]*4
virt_temp=temp/(1. -(e/p)*(1.-eps))

