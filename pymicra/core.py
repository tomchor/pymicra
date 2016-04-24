#!/usr/bin/python
"""
"""

class dataloggerConf(object):
    """
    This class defines a specific configuration of a datalogger output file

    Parameters:
    ----------
    varNames: list of strings or dict
        If a list: should be a list of strings with the names of the variables. If the variable
        is part if the date, then it should be provided as a datetime directive,
        so if the columns is only the year, its name must be `%Y` and so forth. While
        if it is the date in YYYY/MM/DD format, it should be `%Y/%m/%d`. For more info
        see https://docs.python.org/2/library/datetime.html#strftime-and-strptime-behavior
        If a dict: the keys should be the numbers of the columns and the items should follow
        the rules for a list.

    date_cols: list of strings
        should be the subset of varNames that corresponds to the variables that compose
        the timestamp. If it is not provided the program will try to guess by getting
        all variable names that have a percentage sign (%).

    date_connector: string
        generally not really necessary. It is used to join and then parse the date_cols.

    columns_separator: string
        used to assemble the date. Should only be used if the default char creates conflict. If
        the file is tabular-separated then this should be "whitespace"

    header_lines: int
        up to which line of the file is a header. See pandas.read_csv header option.

    first_time_skip: int
        how many units of frequency the first line of the file is offset (generally zero)

    filename_format: string
        tells the format of the file with the standard notation for date and time and with variable
        parts as "?". E.g. if the files are 56_20150101.csv, 57_20150102.csv etc filename_format should be:
            ??_%Y%m%d.csv
        this is useful primarily for the quality control feature

    units: dictionary
        very important: a dictionary whose keys are the columns of the file and whose items are
        the units in which they appear.

    description: string
        brief description of the datalogger configuration file
    """
    def __init__(self, varNames,
            date_cols=None,
            frequency=None,
            date_connector='-', 
            columns_separator=',',
            header_lines=None,
            first_time_skip=False,
            units=None,
            skiprows=None,
            filename_format='CSV_?????.fluxo_??m_%Y_%j_%H%M.dat',
            description='generic datalogger configuration file'):
        #-------------------------
        # Makes sure that units is a dictionary type
        if units is not None:
            if not isinstance(units, dict):
                raise TypeError('units should be a dictionary. Ex.: {"u" : "m/s", "v" : "m/s", "theta" : "K" }')
        #-------------------------

        self.varNames=varNames
        #-------------------------
        # If date_cols is not provided, this will try to guess the date columns
        # by assuming that no other columns has a % sign on their name
        if date_cols:
            self.date_cols=date_cols
        else:
            if type(varNames) == list:
                self.date_cols = [ el for el in varNames if '%' in el ]
            if type(varNames) == dict:
                self.date_cols = { k : it for (k, it) in varNames.iteritems() if '%' in it }
        #-------------------------

        self.frequency=frequency
        self.date_connector=date_connector
        self.columns_separator=columns_separator
        self.header_lines=header_lines
        self.skiprows=skiprows
        self.first_time_skip=first_time_skip
        self.units=units
        self.description=description
        self.filename_format=filename_format


class siteConstants(object):
    """
    Keeper of the characteristics and constants of an experiment.

    Attributes:
    -----------
        variables_height: float
            the main height of the instruments in meters. Generally the height of the sonic anemometer
        canopy height: float
            the mean height of the vegetation meters
        displacement_height: float
            also called zero-plane displacement. Will be estimated as 2/3*canopy_height if not given.
    """
    def __init__(self, variables_height, canopy_height,
             displacement_height=None, z0=None, description=None):

        self.description=description
        self.variables_height = variables_height    #meters
        self.canopy_height = canopy_height          #meters
        self.z0 = z0
        self.description = description
        if displacement_height==None:
            self.displacement_height = (2./3.)*self.canopy_height #meters
        else:
            self.displacement_height=displacement_height


class get_notation(object):
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


