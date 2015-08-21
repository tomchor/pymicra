#!/usr/bin/python
"""
Author: Tomas Chor
Date: 2015-08-07
-------------------------

This module works with micrometeorological data using pandas, numpy, datetime and several other packages

-------------------------

Modifications:


"""
import pandas as pd
import numpy as np

def rotCoor(data, varNames=['u','v','w']):
	"""
	Rotates the coordinates of wind data
	----------------------

	Parameters
	----------
	data: pandas DataFrame 
	the dataFrame to be rotated

	varNames: list
	a list with the names of columns which define the wind speed, typically ['u','v','w']

	"""
	from numpy import cos,sin,zeros,dot
	from math import atan2, sqrt
	DC= np.zeros((3,3))
	wind_vec = data[varNames].mean().values
	m_u, m_v, m_w = wind_vec
	alfa = atan2(m_v,m_u)
	beta =-atan2(m_w, sqrt((m_u**2.)+(m_v**2.)))
	# definition of rotation matrix
	DC[0,0],DC[0,1],DC[0,2] = np.cos(alfa)*np.cos(beta), np.cos(beta)*np.sin(alfa),-np.sin(beta)
	DC[1,0],DC[1,1],DC[1,2] =-np.sin(alfa)          , np.cos(alfa)          , 0.
	DC[2,0],DC[2,1],DC[2,2] = np.cos(alfa)*np.sin(beta), np.sin(alfa)*np.sin(beta), np.cos(beta)
	# application of rotation as a matrix product
	data[[4,5,6]] = np.dot(DC, data[[4,5,6]].values.T).T
	return data


def average(data, mode='moving', rule='10min', **kwargs):
	"""
	Wrapper to return the moving/block average of a given data
	-------------

	Parameters
	----------

	data: pandas DataFrame 
	the dataFrame to be rotated

	mode: string
	mode of average to apply. Currently {'moving', 'block'}.

	rule: string
	pandas offset string to define the block in the block average. Default is "10min".

	"""
	if mode.lower()=='moving':
		return pd.rolling_mean(data, **kwargs)
	elif mode.lower()=='block':
		return data.resample(rule, how='mean', **kwargs)
	else:
		raise KeyError('Mode defined is not correct. Options are "moving" and "block".')

def spectrum(data, col, frequency=10, absolute=True):
	"""
	Author: Tomas Chor

	Calculates the spectrum for a set of data
	"""
	sig=data[col].values
	spec= np.fft.rfft(sig)
	freq= np.fft.rfftfreq(len(sig), d=1./frequency)
	if absolute==True:
		spec=map(abs,spec)
	aux=pd.DataFrame( data={'spectrum':spec}, index=freq )
	aux.index.name='frequencies'
	return aux

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
		flux[str(a)+'-'+str(b)]=data[a]-data[b]
        return flux


#------------------------------------
#
#------------------------------------
class siteConstants(object):
	"""
	Author: Tomas Chor

	Keeper of the characteristics and constants of an experiment.

	Attributes:
	-----------
		variables height: should be a dict with the keys being the timeSeries' variables
		canopy height: should be a float
		displacement_height: should be a float (will be estimated as 2/3*canopy_height if not given)
		variables_units: a dict with the units for every variable to which units are aplicable
		gravity: the acceleration of gravity in m/(s*s) (assumed to be 9.81 if not provided)
		R: Universal gas constant
		Rs: Specific gas constant for dry air
		Rv: Specific gas constant for water vapor
		mu: Rv/Rs
	"""
	from physics import constants

	def __init__(self, variables_height, canopy_height,
		     displacement_height=None,
		     variables_units=None,
		     gravity=9.81,
		     Rs=287.04, 
		     Rv=461.50, 
		     cp=1003.5):
		# Checks if the variables_height is actually a Dict
		if not isinstance(variables_height, dict):
			raise TypeError('variables_height should be a dictionary. Ex.: {"u" : 10, "v" : 10, "theta" : 12 }')
		self.variables_height = variables_height #meters
		self.canopy_height = canopy_height         #meters
		if displacement_height==None:
			self.displacement_height = (2./3.)*self.canopy_height #meters
		else:
			self.displacament_height=displacement_height
		self.gravity = gravity        #meters/(s**2)
		self.Rs = Rs	#specific gas constant for dry air J/(kg.K)
		self.Rv = Rv	#J/(kg.K)
		self.mu = self.Rv/self.Rs
		self.cp = cp     #specific heat for constant pressure J/(kg.K)


