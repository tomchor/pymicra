import pandas as pd

def _integrate_series(self, how='trapz', dateindex=False, **kwargs):
    '''
    Numerically a given series.

    Parameters:
    -----------
    how: string
        the method to use (trapz by default)
        Available methods:
         * trapz - trapezoidal
         * cumtrapz - cumulative trapezoidal
         * simps - Simpson's rule
         * romb - Romberger's rule
    dateindex: bool
        whether or not to assume index is sequence of dates

    See http://docs.scipy.org/doc/scipy/reference/integrate.html for the method details.
    or the source code
    https://github.com/scipy/scipy/blob/master/scipy/integrate/quadrature.py
    '''
    import numpy as np
    from scipy import integrate

    available_rules = set(['trapz', 'cumtrapz', 'simps', 'romb'])
    if how in available_rules:
        rule = integrate.__getattribute__(how)
    else:
        print('Unsupported integration rule: %s' % (how))
        print('Expecting one of these sample-based integration rules: %s' % (str(list(available_rules))))
        raise AttributeError
    if dateindex:
        result = rule(self.values, self.index.astype(np.int64), **kwargs)
    else:
        result = rule(self.values, np.array(self.index), **kwargs)
        #result = rule(self.values, self.index.astype(np.float64), **kwargs)
    return result


#----------
# Definition of the .integrate method for dataframes
from .. import decorators
@decorators.pdgeneral(convert_out=True)
def _integrate_df(self, how='trapz', dateindex=False, **kwargs):
    '''
    Numerically integrate a given dataframe

    Parameters:
    -----------
    how: string
        the method to use (trapz by default)
        Available methods:
         * trapz - trapezoidal
         * cumtrapz - cumulative trapezoidal
         * simps - Simpson's rule
         * romb - Romberger's rule
    dateindex: bool
        whether or not to assume index is sequence of dates

    If you get a ValueError when using cumtrapz scheme, use the initial keyword.

    See http://docs.scipy.org/doc/scipy/reference/integrate.html for the method details.
    or the source code
    https://github.com/scipy/scipy/blob/master/scipy/integrate/quadrature.py
    '''
    import numpy as np
    from scipy import integrate

    df=self.copy()
    available_rules = set(['trapz', 'cumtrapz', 'simps', 'romb'])
    if how in available_rules:
        rule = integrate.__getattribute__(how)
    else:
        print('Unsupported integration rule: %s' % (how))
        print('Expecting one of these sample-based integration rules: %s' % (str(list(available_rules))))
        raise AttributeError
    if dateindex:
        xaxis=df.index.astype(np.int64)
    else:
        xaxis=np.array(df.index)
        #xaxis=df.index.astype(np.float64)
    for c in df.columns:
        df[c] = rule(df[c].values, xaxis, **kwargs)
    if how in ['trapz', 'simps']:
        return df.iloc[0]
    else:
        return df
#pd.Series.integrate = _integrate_series
pd.Series.integrate = _integrate_df
pd.DataFrame.integrate = _integrate_df
#----------

#----------
# Definition of the .differentiate method for dataframes
def _diff_df(self, how='central', axis=0):
    """
    Differentiates self
    """
    import numpy as np
    import auxiliar as aux
    import pandas as pd

    #---------
    # Choose x axis by column or index, considering axis keyword
    if axis==0:
        x = np.array(self.index)
    elif axis==1:
        x = np.array(self.transpose().index)
    else:
        raise ValueError('Axis must be 0 or 1')
    dx = np.diff(x)
    #---------

    #---------
    if how=='central':
        ix = x[1:-1]
        diff=aux.diff_central
    elif how=='fwd':
        ix=x[:-1]
        diff=lambda x,y: y/x
    else:
        raise NameError
    #---------

    #---------
    # Choose what to iterate based on axis argument
    if axis==0:
        out = pd.DataFrame(index=ix, columns=self.columns)
        for c in self.columns:
            dydx = diff(x, np.array(self[c]))
            out[c] = dydx
    elif axis==1:
        out = pd.DataFrame(columns=self.index, index=ix)
        for c in self.transpose().columns:
            dydx = diff(x, np.array(self.transpose()[c]))
            out[c] = dydx
        out = out.transpose()
    else:
        raise ValueError('Axis must be 0 or 1')
    #---------
    return out
pd.DataFrame.differentiate = _diff_df
#----------

del pd
