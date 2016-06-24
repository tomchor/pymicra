

def add(elems, units):
    """
    """
    import pandas as pd

    idx = elems[0].index
    result = pd.Series(index=idx)
    for elem, unit in zip(elems, units):
        elem = elem.reindex(idx)
        result += elem.values*



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

