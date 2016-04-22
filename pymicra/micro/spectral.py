
def hfc_Dias_ea_16(cross_spec, T):
    """
    Applies correction to high frequencies using the quadrature

    Parameters:
    -----------
    cross_spec: series of dataframe
        the cross spectrum whose coespectrum you'd like to correct
    T: float
        response time to use
    """
    import numpy as np

    Co = cross_spec.apply(np.real)
    Qu = cross_spec.apply(np.imag)
    n = np.array(cross_spec.index)
    return Co + 2.*np.pi*T*Qu.multiply(n, axis=0)

def hfc_Massman_Ibrom_08(df):
    pass
    return


def hfc_zeroQuad(slow_spec, freqs, T):
    """
    Applies a correction factor to the spectrum of a slow-measured variable
    based on the response-time T. Quadrature must be analytically zero!

    Parameters:
    -----------
    slow_spec: numpy.array or series
        the spectrum to be corrected (not cross-spectrum!)
    freqs: numpy.array
        frequencies
    T: float
        the response-time

    Returns:
    --------
    spec: numpy.array
        the recovered array
    """
    import numpy as np 

    assert len(slow_spec)==len(freqs)
    return slow_spec*(1.0 + 4.*(np.pi**2.)*(freqs**2.) * (T**2.))


def Ogive(df, no_nan=True):
    """
    Integrates the Ogive from Coespectra

    Parameters:
    -----------
    df: dataframe
        cospectrum to be integrated
    """
    import scipy as sp
    import numpy as np
    import pandas as pd

    if no_nan:
        x=np.array(df.index)
        out=pd.DataFrame(index=x)
        for c in df.columns:
            y=np.array(df[c])
            og=sp.integrate.cumtrapz(y[::-1], x=x[::-1], initial=0)
            out.loc[:,c]=-og[::-1]
    else:
        dfs=[]
        for c in df.columns:
            ser=df[c].dropna()
            x=np.array(ser.index)
            y=np.array(ser)
            og=sp.integrate.cumtrapz(y[::-1], x=x[::-1], initial=0)
            dfs.append( pd.DataFrame(data=-og[::-1], columns=[c], index=x) )
        for d in dfs[1:]:
            dfs[0]=dfs[0].join(d, how='outer')
        out=dfs[0]
    return out


def recspe(slow_spec, freqs, T):
    """
    Applied a correction factor to the spectrum of a slow-measured variable
    based on the response-time T

    Parameters:
    -----------
    freqs: numpy.array
        frequencies
    slow_spec: numpy.array
        the spectrum
    T: float
        the response-time
    """
    import numpy as np

    assert len(slow_spec)==len(freqs)
    return slow_spec*(1.0 + 4.*(np.pi**2.)*(freqs**2.) * (T**2.))


def recspeAux(df, T):
    """
    Wrapper to make recspe work in a pandas.DataFrame
    """
    import numpy as np

    freqs=np.array(df.index)
    for c in df.columns:
        spec =np.array(df[c])
        df[c]=recspe(spec, freqs, T)
    return df


