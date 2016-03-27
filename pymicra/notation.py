"""
This module sets the default notation options in pymicra

"""


class get_notation2:
    mean='%s_mean'
    fluctuation="%s'"
    star='%s_star'
    concentration='conc_%s'
    density='rho_%s'
    molar_density='mrho_%s'
    mixing_ratio='r_%s'

    u='u'
    v='v'
    w='w'
    thermodyn_temp='theta'
    virtual_temp='theta_v'
    pressure='p'
    relative_humidity='rh'
    specific_humidity='q'

    h2o='h2o'
    co2='co2'
    ch4='ch4'
    o3 = 'o3'
    moist_air = 'air'
    dry_air   = 'dry'

    similarityVar = 'zeta'
    MonObuLen = 'Lm'
    momentum_flux = 'tau'
    sensible_heat_flux = 'H'
    virtual_sensible_heat_flux = 'Hv'
    latent_heat_flux = 'E'
    flux_of = 'F_%s'

class get_notation:
    mean_preffix=''
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
    specific_humidity='q'

    h2o='h2o'
    co2='co2'
    ch4='ch4'
    o3 = 'o3'
    moist_air = 'air'
    dry_air   = 'dry'

    similarityVar = 'zeta'
    MonObuLen = 'Lm'
    momentum_flux = 'tau'
    sensible_heat_flux = 'H'
    virtual_sensible_heat_flux = 'Hv'
    latent_heat_flux = 'E'
    solute_flux = 'F_{}'
