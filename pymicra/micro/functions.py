#!/usr/bin/python
"""
Author: Tomas Chor
Date: 2015-08-07
-------------------------

This module works with micrometeorological data using pandas, numpy, datetime and several other packages

-------------------------

TO-DO LIST

Maybe add all these variables that EddyPro calculates:
http://www.licor.com/env/help/EddyPro3/Content/Topics/Calculating_Micromet_Variables.htm

"""
import numpy as np
from .. import notation


def phi_c(zeta, x=None):
    """
    The nondimensional standard deviation function, defined as:

    phi_c(zeta) = sigma_c / c_star

    From zahn.ea
    """
    if zeta<0:
        if x=='u' or x=='w':
            return 1.25*(1. - 3.*zeta)**(1./3.)
        else:
            zeta = abs(zeta)
            return 2.*(1. + 9.5*zeta)**(-1./3.)
    else:
        if x=='u' or x=='w':
            return 1.25
        else:
            return 2. 




def phi(zeta, x=None, C_unstable=None, C_stable=None):
    """
    The nondimensional gradients, defined as:

    phi_F(zeta) = kappa*(z-d)*dCdz/c_star
    phi_H(zeta) = kappa*(z-d)*dTdz/T_star
    phi_E(zeta) = kappa*(z-d)*dqdz/q_star

    Currently using Businger-Dyer eqs.

    TODO: generalize coefficients
    """
    if zeta<0:
        if x=='tau':
            return (1. - 16.*zeta)**(-1./4.)
        else:
            return (1. - 16.*zeta)**(-1./2.)
    if zeta>=0:
        if zeta<1.:
            return 1. + 5.*zeta
        if zeta>=1.:
            return 6.



def Psi(zeta, x='tau', zeta0=0.):
    """
    Integral Monin-Obukhov scale or deviation function, which is the deviation of a
    variable (x) in relation nto their logarithmic profiles due the stability zeta != 0

    Taken from Simpson.ea1998--the.validity.of.similarity.theory;.in.the.roughness.sublayer

    Parameters:
    -----------
    zeta: float
        the stability variable
    x: string
        the variable. Options are 'tau', 'H', 'E', 'F'
    zeta0: float
        value of zeta_{0 x}. Only used for the stable case
    """
    #---------
    # For unstable conditions
    if zeta < 0:
        b = (1. - 16.*zeta)**(1./4.)
        if x=='tau':
            psi = np.log((b*b + 1.)/2.) + 2.*np.log((b+1.)/2.) - 2.*np.arctan(b) + np.pi/2.
        if any([ x == var for var in ['H', 'E', 'F'] ]):
            psi = 2.*np.log((b*b + 1.)/2.)
        return psi
    #---------

    #---------
    # For stable conditions
    else:
        return -5.*(zeta - zeta0)
    #---------


def ste(data, w_fluctuations="w'"):
    """
    Returns the Symmetric Transfer Efficiency in the time domain, ste
    according to Cancelli, Dias, Chamecki. Dimensionless criteria for the production-dissipation equilibrium
    of scalar fluctuations and their implications for scalar similarity, Water Resources Research, 2012

    NEEDS TO BE VALIDATED!
    """
    from data import bulkCorr
    wcol=w_fluctuations
    df=data.copy()
    w=df[wcol]
    df=df.drop(wcol, axis=1)
    a,b=df.columns
    a,b=df[a],df[b]
    rwa=bulkCorr([w,a])
    rwb=bulkCorr([w,b])
    rwa=np.abs(rwa)
    rwb=np.abs(rwb)
    return 1. - np.abs( rwa - rwb )/( rwa + rwb )
 



#---------------------
# Not tested out yet
#---------------------
def _Cx(x, za, zb, d, z0, Lm, Psif=None):
    """
    Get Cx for gradient-flux method.
    In this code, zb means tau and za means anything else and the bottom level
    is the ground (height zero)

    Parameters:
    -----------

    x: string
        {'tau', 'H', 'E', 'c', 'F'}
    za: float
        height at which substance (water, c2o, temperature, etc) is measured
    zb: float
        height at which wind velocity is measured
    d: float
        displacement height
    z0: float
        roughness length
    Lm: float
        Monin-Obukhov Length
    Psif: function (optional)
        deviation function Psi=log(abs(zeta)) - Phi(zeta)
    """
    #--------
    # Checks for Psi and if it's a function
    if Psif==None:
        print 'Should get a repository for functions'
    elif hasattr(Psif, '__call__')==False:
        raise TypeError('Psi argument has to be a function')
    else:
        Psif = Psi
    #--------

    if x=='tau':
        cx=1./(auxLogMinPsi(zb, d, z0, Psif, Lm))**2.
    elif x=='H':
        pass
    elif x=='E':
        pass
    elif x=='F':
        cx = 1./( auxLogMinPsi(zb, d, z0, Psif, Lm) * auxLogMinPsi(za, d, z0, Psif, Lm) )
    else:
        raise NameError('x is not on the list. See help(Cx) for more information')
    cx*=kappa**2
    return cx



def _auxLogMinPsi(z_i, d, z0, Psi, Lm):
    """
    Auxiliar function that is used in the denominator of the C_tau-like coefficients

    Taken from Dias, Micrometeorological Approaches to to estimate fluxes of greenhouse gases between the surface and the
    atmosphere, (Chapter 7, eqs. 94-97, p. 22.)

    CONSIDER TAKING THIS SOMEWHERE ELSE!
    """
    if hasattr(Psi, '__call__')==False:
        raise TypeError('Psi argument has to be a function')
    a=np.log((z_i-d)/z0)
    b=Psi((z_i-d)/Lm)
    return a-b




