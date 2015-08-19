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
from sys import exit

def readDataFile(fname, varNames=None, **kwargs):
	"""
	Author: Tomas Chor

	Parameters
	----------
	kwargs: dict
		dictionary with kwargs of pandas' read_csv function
		see http://pandas.pydata.org/pandas-docs/stable/generated/pandas.read_csv.html for more detail

	varNames: list
		list containing the names of each variable in the file.
		
	Returns
	---------
	dataFrame: pandas.DataFrame object
	"""
	data=pd.read_csv(fname, **kwargs)
	if varNames:
		data.columns=varNames + list(data.columns[len(varNames):])
	return data


def readDataFiles(flist, **kwargs):
	"""
	Author: Tomas Chor
	** needs to be tested! **

	-------------------------------------

	* kwargs are readDataFile kwargs

	* returns one pandas.DataFrame
	"""
	data=pd.DataFrame()
	for f in flist:
		subdata=readDataFile(f, **kwargs)
		data=pd.concat( [data, subdata], ignore_index=True)
	return data

def parseDates(data, date_cols, separator='-', first_time_skip=1):
	"""
	Author: Tomas Chor
	date: 2015-08-10
	This routine parses the date from a pandas DataFrame when it is divided into several columns
	----------------------------------

	data: pandas DataFrame

	date_cols: list of strings
	A list of the names of the columns in which the date is divided
	the naming of the date columns must be in accordance with the datetime directives,
	so if the first column is only the year, its name must be `%Y` and so forth.
	see https://docs.python.org/2/library/datetime.html#strftime-and-strptime-behavior
	"""
	from datetime import timedelta,datetime
	#------------------------------------
	# joins the names of the columns, which must match the datetime directive (see __doc__)
	#------------------------------------
	date_format=separator.join(date_cols)
	auxformat='%Y-%m-%d %H:%M:%S.%f'
	#-------------------------------------
	# joins the appropriate pandas columns because pandas can read only one column into datetime
	#-------------------------------------
	aux=data[date_cols[0]].astype(str)
	for col in date_cols[1:]:
		aux+=separator + data[col].astype(str)
	dates=pd.to_datetime(aux, format=date_format)
	#-------------------------------------
	# The next steps are there to check if there are fractions that are not expressed in the datetime convention
	# and it assumes that the lowest time period expressed is the minute
	#-------------------------------------
	first_date=dates.unique()[1]
	n_fracs=len(dates[dates.values==first_date])
	print 'Identified that each minute contains', n_fracs, 'fractions'
	dates=[ date.strftime(auxformat) for date in dates ]
	aux=dates[0]
	cont=first_time_skip
	for i,date in enumerate(dates):
		if date==aux:
			pass
		else:
			cont=0
			aux=date
		dates[i]=datetime.strptime(date, auxformat) + timedelta(minutes=cont/float(n_fracs))
		cont+=1
	#-------------------------------------
	# setting new dates list as the index
	#-------------------------------------
	data=data.set_index([dates])
	#-------------------------------------
	# removing the columns used to generate the date
	#-------------------------------------
	data=data[ [col for col in data.columns if col not in date_cols] ]
	return data


#------------------------------------
#
#------------------------------------
class dataloggerConf(object):
	"""
	This class defines a specific configuration of a datalogger output file
	--------------------------

	Parameters:
	----------

	varNames: list of strings
	should be a list of strings with the names of the variables. If the variable
	is part if the date, then it should be provided as a datetime directive,
	so if the columns is only the year, its name must be `%Y` and so forth. While
	if it is the date in YYYY/MM/DD format, it should be `%Y/%m/%d`. For more info
	see https://docs.python.org/2/library/datetime.html#strftime-and-strptime-behavior

	date_cols: list of strings
	should be the subset of varNames that corresponds to the variables that compose
	the timestamp. If it is not provided the program will try to guess by getting
	all variable names that have a percentage sign.

	date_connector: string
	generally not really necessary. It is used to join and then parse the date_cols.

	columns_separator: string
	"""
	def __init__(self, varNames,
		     date_cols=None,
		     frequency=None,
		     date_connector='-', 
		     columns_separator=',',
		     header_lines=None,
		     description='generic datalogger configuration file'):

		self.varNames=varNames
		if date_cols:
			self.date_cols=date_cols
		else:	#tries to guess the date columns by assuming that no other columns has a % sign on their name
			date_cols=[ el for el in varNames if '%' in el ]
		self.frequency=frequency
		self.date_connector=date_connector
		self.columns_separator=columns_separator
		self.header_lines=header_lines
		self.file_names=file_names
		self.description=description


class timeSeries(object):
	"""
	Creates a micrometeorological time series from a file or list of files.
	It needs a dataloggerConf object.
	"""
	def __init__(self, dataloggerConf, index_by_date=True):
		self.description='Data read with dataloggerConf: ',dataloggerConf.description
