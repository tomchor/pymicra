"""
This module sets the default options in pymicra

STILL NEEDS TO BE TESTED AND OPTIMIZED
"""
#print locals()

mean_shit='original'

def set_defaults(color=None):
#global _notation
#color=None
#class _notation(object):
    global _notation
    class get_notation:
        pass
    _notation=get_notation()

    _notation.mean_preffix=''
    mean_suffix='_mean'
    
    fluctuation_preffix=""
    fluctuation_suffix="'"
    
    star_preffix=''
    star_suffix='_star'
    
    concentration_preffix='conc_'
    concentration_suffix=''
    
    density_preffix='rho_'
    density_suffix=''
    
    molar_density_preffix='mrho_'
    molar_density_suffix=''
    
    mixing_ratio_preffix='r_'
    mixing_ratio_suffix=''
    
    u='u'
    v='v'
    w='w'
    thermodyn_temp='theta'
    virtual_temp='theta_v'
    pressure='p'
    relative_humidity='rh'

    if color==None:
        _color='blue'
