from .. import decorators

def multiply(elems, units, inplace=False, unitdict=None, key=None):
    """
    Multiply elements considering their units
    """
    return operate(elems, units, inplace=inplace, unitdict=unitdict, key=key, operation='*')

def add(elems, units, inplace=False, unitdict=None, key=None):
    """
    Add elements considering their units
    """
    return operate(elems, units, inplace=inplace, unitdict=unitdict, key=key, operation='+')

def divide(elems, units, inplace=False, unitdict=None, key=None):
    """
    Divide elements considering their units
    """
    return operate(elems, units, inplace=inplace, unitdict=unitdict, key=key, operation='/')



def operate(elems, units, inplace=False, unitdict=None, key=None, operation='+'):
    """
    Operate on elements considering their units

    Parameters:
    -----------
    elems: list, tuple
        list of pandas.Series
    units: list, tuple
        list of pint.units ordered as the elems list
    inplace: bool
        sets dictionary inplace
    unitdict: dict
        dict to be set inplace
    key: str
        name of variables to be set inplace as dict key
    """
    import pandas as pd
    import numpy as np
    from .. import ureg

    idx = elems[0].index

    if operation=='+':
        result = elems[0].values*units[0]
        for elem, unit in zip(elems[1:], units[1:]):
            if type(elem) == pd.Series:
                elem = elem.reindex(idx)
                result += elem.values*unit
            else:
                result += elem*unit

    if operation=='*':
        result = elems[0].values*units[0]
        for elem, unit in zip(elems[1:], units[1:]):
            if type(elem) == pd.Series:
                elem = elem.reindex(idx)
                result *= elem.values*unit
            else:
                result *= elem*unit

    if operation=='/':
        result = elems[0].values*units[0]
        for elem, unit in zip(elems[1:], units[1:]):
            if type(elem) == pd.Series:
                elem = elem.reindex(idx)
                result /= elem.values*unit
            else:
                result /= elem*unit

    out = pd.Series(result.magnitude)
    funit = ureg.Quantity(1, result.units)

    if inplace==True:
        unitdict.update({key : funit})
        return out
    else:
        return out, funit



def parseUnits(unitstr):
    """
    Gets unit from string, list of strings, or dict's values, using the UnitRegistry
    defined in __init__.py
    """
    from .. import ureg

    if isinstance(unitstr, str):
        return ureg[unitstr].u
    elif isinstance(unitstr, list):
        return [ ureg[el].u for el in unitstr ]
    elif isinstance(unitstr, dict):
        return { key: ureg[el].u for key, el in unitstr.items() }


def convert_to(data, inunit, outunit, inplace=False, key=None): 
    ''' 
    Converts data from one unit to the other 
 
    Parameters: 
    ----------- 
    data: pandas.series 
        to be chanhed from one unit to the other 
    inunit: pint.quantity or dict 
        unit(s) that the data is in 
    outunit: str 
        convert to this unit 
    inplace: bool 
        if inunit is a dict, the dict is update in place. "key" keyword must be provided 
    key: str 
        if inunit is a dict, it is the name of the variable to be changed 
    ''' 
    from .. import Q_ 
 
    if key: 
        Q = inunit[key].to(outunit) 
    else: 
        Q = inunit.to(outunit) 
 
    coef = Q.magnitude 
    outunit = Q_(1, Q.units) 
    if inplace: 
        inunit.update({key : outunit}) 
        return data*coef 
    else: 
        return data*coef, outunit 


#@decorators.pdgeneral(convert_out=True)
def convert_cols(data, guide, units, inplace=False):
    ''' 
    Converts data from one unit to the other 
 
    Parameters: 
    ----------- 
    data: pandas.DataFrame
        to be chanhed from one unit to the other 
    guide: dict
        {names of columns : units to converted to}
    units: dict
        units dictionary
    inplace: bool 
        if inunit is a dict, the dict is update in place. "key" keyword must be provided 
    ''' 
    from .. import algs
    from .. import Q_

    data = data.copy()

    #-------
    # An attempt to make it work with Series
    if len(data.columns)==1 and (type(guide) != dict):
        guide = { data.columns[0] : guide }
    guide = algs.parseUnits(guide)
    #-------

    #-------
    # We first turn it into a numpy array to make the conversion using pint natively
    for col, outunit in guide.iteritems():
        aux = Q_(data[ col ].values, units[ col ])
        aux = aux.to(outunit)
        data.loc[:, col] = aux
    #-------
        
    if inplace: 
        units.update(guide) 
        return data 
    else: 
        return data, guide 

def convert_indexes(data, guide, units, inplace=False):
    ''' 
    Converts data from one unit to the other 
 
    Parameters: 
    ----------- 
    data: pandas.Series
        to be chanhed from one unit to the other 
    guide: dict
        {names of columns : units to converted to}
    units: dict
        units dictionary
    inplace: bool 
        if inunit is a dict, the dict is update in place. "key" keyword must be provided 
    ''' 
    from .. import algs

    data = data.copy()
    guide = algs.parseUnits(guide)

    #-------
    # We first turn it into a numpy array to make the conversion using pint natively
    for idx, outunit in guide.iteritems():
        aux = data[ idx ] * units[ idx ]
        aux = aux.to(outunit)
        data.loc[ idx ] = aux.magnitude
    #-------
        
    if inplace: 
        units.update(guide) 
        return data 
    else: 
        return data, guide 

