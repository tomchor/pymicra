
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
            elem = elem.reindex(idx)
            result += elem.values*unit

    if operation=='*':
        result = elems[0].values*units[0]
        for elem, unit in zip(elems[1:], units[1:]):
            elem = elem.reindex(idx)
            result *= elem.values*unit

    if operation=='/':
        result = elems[0].values*units[0]
        for elem, unit in zip(elems[1:], units[1:]):
            elem = elem.reindex(idx)
            result /= elem.values*unit

    out = pd.Series(result.magnitude)
    funit = ureg.Quantity(1, result.units)

    if inplace==True:
        unitdict.update({key : funit})
        return out
    else:
        return out, funit



def parseUnits(unitstr):
    '''
    Gets unit from string, list of strings, or dict's values, using the UnitRegistry
    defined in __init__.py
    '''
    try:
        from .. import ureg, Q_
    except ImportError:
        print('You should have pint installed to use units!')
        return unitstr

    if isinstance(unitstr, str):
        return ureg[unitstr]
    elif isinstance(unitstr, list):
        return [ ureg[el] for el in unitstr ]
    elif isinstance(unitstr, dict):
        return { key: ureg[el] for key, el in unitstr.items() }


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
        if inunit is a dict, it is the name of the variable to be chamged 
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

