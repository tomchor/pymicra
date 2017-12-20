"""
Defines useful decorators for Pymicra
"""
from __future__ import absolute_import, print_function, division
from functools import wraps


def pdgeneral_in(func):
    """Defines a decorator that transforms Series into DataFrame"""
    import pandas as pd
    @wraps(func)
    def wrapper(*args, **kwargs):
        if isinstance(args[0], pd.Series):
            ser = args[0]
            args = tuple([pd.DataFrame(ser, dtype=ser.dtype)]) + args[1:]
        elif not isinstance(args[0], pd.DataFrame):
            raise TypeError('Input must be a pandas.Series or pandas.DataFrame')
        result = func(*args, **kwargs)
        return result
    return wrapper


def pdgeneral_io(func):
    """
    If the input is a series transform it to a dtaframe then transform the output from dataframe
    back into a series. If the input is a series and the output is a one-element series, transform it to a float.

    Currently the output functionality works only when the output is one variable, not an array
    of elements.
    """
    import pandas as pd

    @wraps(func)
    def wrapper(*args, **kwargs):
        if isinstance(args[0], pd.Series):
            ser = args[0]
            args = tuple([pd.DataFrame(ser, dtype=ser.dtype)]) + args[1:]

            result = func(*args, **kwargs)

            if isinstance (result, pd.DataFrame):
                result = pd.Series(result.iloc[:, 0].values, index=result.index)

            elif isinstance (result, pd.Series):
                if result.size==1:
                    result = float(result)
                else:
                    pass

            elif type(result) in [ list, tuple ]:
                outser = pd.Series(results[0].iloc[:, 0].values, index=result[0].index)
                if isinstance(result, list):
                    result[0] = outser
                else:
                    result = tuple([outser]) + result[1:]

            return result

        elif not isinstance(args[0], pd.DataFrame):
            raise TypeError('Input must be a pandas.Series or pandas.DataFrame')
        else:
            result = func(*args, **kwargs)
            return result

    return wrapper

def pdgeneral(convert_out=True):
    """Defines a decorator to make functions work on both pandas.Series and DataFrames

    Parameters
    ----------
    convert_out: bool
        if True, also converts output back to Series if input is Series
    """
    if convert_out:
        return pdgeneral_io
    else:
        return pdgeneral_in




def autoassign(*names, **kwargs):
    """Decorator that automatically assigns keywords as atributes
 
    allow a method to assign (some of) its arguments as attributes of
    'self' automatically.  E.g.
    
    To restrict autoassignment to 'bar' and 'baz', write:
    
    @autoassign('bar', 'baz')
    def method(self, foo, bar, baz): ...

    To prevent 'foo' and 'baz' from being autoassigned, use:

    @autoassign(exclude=('foo', 'baz'))
    def method(self, foo, bar, baz): ...
    """
    from functools import wraps
    from inspect import getargspec, isfunction
    #------
    # This is only needed if using Python 2
    #from itertools import izip, ifilter as zip, filter
    #------
    from itertools import starmap

    if kwargs:
        exclude, f = set(kwargs['exclude']), None
        sieve = lambda l: filter(lambda nv: nv[0] not in exclude, l)
    elif len(names) == 1 and isfunction(names[0]):
        f = names[0]
        sieve = lambda l:l
    else:
        names, f = set(names), None
        sieve = lambda l: filter(lambda nv: nv[0] in names, l)
    def decorator(f):
        fargnames, _, _, fdefaults = getargspec(f)
        # Remove self from fargnames and make sure fdefault is a tuple
        fargnames, fdefaults = fargnames[1:], fdefaults or ()
        defaults = list(sieve(zip(reversed(fargnames), reversed(fdefaults))))
        @wraps(f)
        def decorated(self, *args, **kwargs):
            assigned = dict(sieve(zip(fargnames, args)))
            assigned.update(sieve(kwargs.items()))
            for _ in starmap(assigned.setdefault, defaults): pass
            self.__dict__.update(assigned)
            return f(self, *args, **kwargs)
        return decorated
    return f and decorator(f) or decorator
