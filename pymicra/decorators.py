from functools import wraps


def pdgeneral_in(func):
    """
    """
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
    """
    import pandas as pd
    @wraps(func)
    def wrapper(*args, **kwargs):
        if isinstance(args[0], pd.Series):
            ser = args[0]
            args = tuple([pd.DataFrame(ser, dtype=ser.dtype)]) + args[1:]

            result = func(*args, **kwargs)

            if isinstance (result, pd.DataFrame):
                result = pd.Series(result.iloc[:, 0].values)

            elif isinstance (result, pd.Series):
                pass

            elif isinstance(result, list):
                result[0] = pd.Series(result[0].iloc[:, 0].values)

            elif isinstance(result, tuple):
                result = tuple([pd.Series(result[0].iloc[:, 0].values)]) + result[1:]

            return result

        elif not isinstance(args[0], pd.DataFrame):
            raise TypeError('Input must be a pandas.Series or pandas.DataFrame')
        else:
            result = func(*args, **kwargs)
            return result
    return wrapper

def pdgeneral(convert_out=True):
    if convert_out:
        return pdgeneral_io
    else:
        return pdgeneral_in



#@pdgeneral(convert_out=True)
#def square(df, n, func=pow):
#    df = df.copy()
#    for c in df.columns:
#        df[c] = df[c]**n
#    return df, n


def autoassign(*names, **kwargs):
    """
    autoassign(function) -> method
    autoassign(*argnames) -> decorator
    autoassign(exclude=argnames) -> decorator
    
    allow a method to assign (some of) its arguments as attributes of
    'self' automatically.  E.g.
    
    >>> class Foo(object):
    ...     @autoassign
    ...     def __init__(self, foo, bar): pass
    ... 
    >>> breakfast = Foo('spam', 'eggs')
    >>> breakfast.foo, breakfast.bar
    ('spam', 'eggs')
    
    To restrict autoassignment to 'bar' and 'baz', write:
    
        @autoassign('bar', 'baz')
        def method(self, foo, bar, baz): ...

    To prevent 'foo' and 'baz' from being autoassigned, use:

        @autoassign(exclude=('foo', 'baz'))
        def method(self, foo, bar, baz): ...
    """
    from functools import wraps
    from inspect import getargspec, isfunction
    from itertools import izip, ifilter, starmap

    if kwargs:
        exclude, f = set(kwargs['exclude']), None
        sieve = lambda l:ifilter(lambda nv: nv[0] not in exclude, l)
    elif len(names) == 1 and isfunction(names[0]):
        f = names[0]
        sieve = lambda l:l
    else:
        names, f = set(names), None
        sieve = lambda l: ifilter(lambda nv: nv[0] in names, l)
    def decorator(f):
        fargnames, _, _, fdefaults = getargspec(f)
        # Remove self from fargnames and make sure fdefault is a tuple
        fargnames, fdefaults = fargnames[1:], fdefaults or ()
        defaults = list(sieve(izip(reversed(fargnames), reversed(fdefaults))))
        @wraps(f)
        def decorated(self, *args, **kwargs):
            assigned = dict(sieve(izip(fargnames, args)))
            assigned.update(sieve(kwargs.iteritems()))
            for _ in starmap(assigned.setdefault, defaults): pass
            self.__dict__.update(assigned)
            return f(self, *args, **kwargs)
        return decorated
    return f and decorator(f) or decorator