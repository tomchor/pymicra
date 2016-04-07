"""
This module sets the default notation options in pymicra

"""


class get_notation2:
    """
    This creates an object that holds the default notation for pymicra.
    Example of usage:

    notation = get_notation()

    fluctuation_of_co2_conc = notation.concentration % notation.fluctuation % notation.co2

    You should be careful with the order. The last argument should not have any '%' symbols
    or you'll get a "TypeError: not all arguments converted during string formatting" message.
    """
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
    water_vapor_flux = 'E'
    latent_heat_flux = 'LE'
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
    water_vapor_flux = 'E'
    latent_heat_flux = 'LE'
    solute_flux = 'F_{}'
