#!/usr/bin/python
"""

"""
import numpy as np
def combine(levels, order='Crescent'):
	""" Combines the elements of the levels in non-repeting pairs
	for a n length array, gives n*(n-1)/2 combinations

	Parameters
	----------
	levels: list
	list whose elements should be combined
	"""
	print order
	if order.lower()=='decrescent':
		levels=levels[::-1]
	combinations=[]
	for i in range(0, len(levels)):
		for j in range(i+1, len(levels)):
			if j>i:
				combinations.append([levels[j],levels[i]])
        return combinations


