class ConvertArgumentTypes(object):
    """
    Converts function arguments to specified types.
    
    Example:
    @ConvertArgumentTypes(int, float, c=int, d=str)
    def z(a,b,c=0,d=""):
            return a + b, (c,d)
    """
    def __init__(self,*args, **kw):
        self.args = args
        self.kw = kw
    def __call__(self, f):
        def func(*args, **kw):
            nargs = [x[0](x[1]) for x in zip(self.args, args)]
            invalidkw = [x for x in kw if x not in self.kw]
            if len(invalidkw) > 0:
                raise TypeError, f.func_name + "() got an unexpected keyword argument '%s'" % invalidkw[0]
            kw = dict([(x,self.kw[x](kw[x])) for x in kw])
            v = f(*nargs, **kw)
            return v
        return func

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
