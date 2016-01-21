import numpy as np
import algs

def hfc_Dias_ea_16(df):
    pass
    return

def hfc_Massman_Ibrom_08(df):
    pass
    return

def hfc_zeroQuad(slow_spec, freqs, T):
    """
    Applies a correction factor to the spectrum of a slow-measured variable
    based on the response-time T. Quadrature must be analytically zero!

    Parameters:
    -----------
    freqs: numpy.array
        frequencies
    slow_spec: numpy.array
        the spectrum
    T: float
        the response-time

    Returns:
    --------
    spec: numpy.array
        the recovered array
    """
    assert len(slow_spec)==len(freqs)
    return slow_spec*(1.0 + 4.*(np.pi**2.)*(freqs**2.) * (T**2.))
