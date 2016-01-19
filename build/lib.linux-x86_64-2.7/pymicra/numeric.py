from scipy import integrate
import pandas as pd
import numpy as np
import notation
import algs

def integrate_series(self, how='trapz', dateindex=False, **kwargs):
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


def integrate_df(self, how='trapz', dateindex=False, **kwargs):
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

    See http://docs.scipy.org/doc/scipy/reference/integrate.html for the method details.
    or the source code
    https://github.com/scipy/scipy/blob/master/scipy/integrate/quadrature.py
    '''
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
pd.Series.integrate = integrate_series
pd.DataFrame.integrate = integrate_df

def diff_df(self, how='central'):
    x= np.array(self.index)
    dx = np.diff(x)
    if how=='central':
        ix = x[1:-1]
        diff=algs.diff_central
    elif how=='fwd':
        ix=x[:-1]
        diff=lambda x,y: y/x
    else:
        raise NameError
    out = pd.DataFrame(index=ix, columns=self.columns)
    for c in self.columns:
        dydx = diff(x, np.array(self[c]))
        out[c] = dydx
    return out

pd.DataFrame.derivate = diff_df

def pcolor():
    print _notation.pressure
