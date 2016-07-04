from __future__ import print_function

class dataloggerConfig(object):
    """
    This class defines a specific configuration of a datalogger output file

    Parameters:
    ----------
    from_file: str
        path of .dlc file (datalogger configuration file) to read from. This will ignore all other
        keywords.

    varNames: list of strings or dict
        If a list: should be a list of strings with the names of the variables. If the variable
        is part if the date, then it should be provided as a datetime directive,
        so if the columns is only the year, its name must be `%Y` and so forth. While
        if it is the date in YYYY/MM/DD format, it should be `%Y/%m/%d`. For more info
        see https://docs.python.org/2/library/datetime.html#strftime-and-strptime-behavior
        If a dict: the keys should be the numbers of the columns and the items should follow
        the rules for a list.

    date_cols: list of ints
        should be indexes of the subset of varNames that corresponds to the variables that compose
        the timestamp. If it is not provided the program will try to guess by getting
        all variable names that have a percentage sign (%).

    date_connector: string
        generally not really necessary. It is used to join and then parse the date_cols.

    columns_separator: string
        used to assemble the date. Should only be used if the default char creates conflict. If
        the file is tabular-separated then this should be "whitespace".

    header_lines: int
        up to which line of the file is a header. See pandas.read_csv header option.

    first_time_skip: int
        how many units of frequency the first line of the file is offset (generally zero).

    filename_format: string
        tells the format of the file with the standard notation for date and time and with variable
        parts as "?". E.g. if the files are 56_20150101.csv, 57_20150102.csv etc filename_format should be:
            ??_%Y%m%d.csv
        this is useful primarily for the quality control feature.

    units: dictionary
        very important: a dictionary whose keys are the columns of the file and whose items are
        the units in which they appear.

    description: string
        brief description of the datalogger configuration file.
    """

    def __init__(self,
            from_file=None,
            varNames=None,
            date_cols=None,
            frequency=None,
            date_connector=None,
            columns_separator=None,
            header_lines=None,
            first_time_skip=False,
            units=None,
            skiprows=None,
            filename_format=None,
            description='Generic datalogger configuration file. Type help(dataloggerConfig) to read intructions'):
        """
        Initiates the class
        """
        import algs as aux

        #-------------------------
        # Makes sure that units is a dictionary type
        if units is not None:
            if not isinstance(units, dict):
                raise TypeError('units should be a dictionary. Ex.: {"u" : "m/s", "v" : "m/s", "theta" : "K" }')
            else:
                units = aux.parseUnits(units)
        #-------------------------

        #-------------------------
        if from_file:
            from io import read_dlc
            dlconf = read_dlc(from_file)
            self.__dict__.update(dlconf.__dict__)
            return
        #-------------------------

        #-------------------------
        if (type(varNames)==list) or (type(varNames)==dict):
            self.varNames=varNames

            #-------------------------
            # If date_cols is not provided, this will try to guess the date columns
            # by assuming that no other columns has a % sign on their name
            if date_cols:
                import numpy as np
                self.date_cols=date_cols
                self.date_col_names = [ varNames[ idx ] for idx in date_cols ]
            else:
                if type(varNames) == list:
                    self.date_col_names = [ el for el in varNames if '%' in el ]
                    self.date_cols = aux.get_index(varNames, self.date_cols)
                if type(varNames) == dict:
                   date_cols = { k : it for (k, it) in varNames.iteritems() if '%' in it }
                   self.date_col_names = date_cols.values()
                   self.date_cols = date_cols.keys()
           #-------------------------

        else:
            self.varNames = varNames
            self.date_cols = date_cols
            self.date_col_names = None
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


    def __str__(self):
        return '<pymicra.dataloggerConfig>\n{}'.format(self.description)
    __repr__ = __str__


class siteConfig(object):
    """
    Keeper of the configurations and constants of an experiment. (such as height of instruments,
    location, canopy height and etc)
    """
    from . import decorators
    @decorators.autoassign
    def __init__(self, from_file=None,
             instruments_height=None, measurement_height=None, canopy_height=None,
             displacement_height=None, roughness_length=None,
             latitude=None, longitude=None, altitude=None, description=None):
        """
        Parameters:
        -----------
            from_file: str
                path to .site file with the configurations of the experiment. All atributes are taken from there.
            measurement_height: float
                the mainn height of the instruments in meters considered for calculations.
                Generally it's the height of the sonic anemometer.
            canopy_height: float
                the mean height of the vegetation meters
            displacement_height: float
                also called zero-plane displacement. Will be estimated as 2/3*canopy_height if not given.
            roughness_length:
                the average length of roughness elements.
            description:
                short description of the site
        """

        #---------
        # If from file keyword exists then we get it from there
        if from_file:
            from io import read_site
            siteconf = read_site(from_file)
            self.__dict__.update(siteconf.__dict__)
            return
        #---------

        #self.description = description
        #self.instruments_height = instruments_height    #meters
        #self.measurement_height = measurement_height    #meters
        #self.canopy_height = canopy_height              #meters
        #self.roughness_length = roughness_length

        #---------
        # If there's no displacement height, we try to calculate it
        if displacement_height != None:
            self.displacement_height=displacement_height
        else:
            if canopy_height:
                self.displacement_height = (2./3.)*self.canopy_height #meters
            else:
                self.displacement_height = None
        #---------
            
        #self.latitude = latitude
        #self.longitude = longitude
        #self.altitude = altitude


    def __str__(self):
        """
        Creates nice representation for printing siteConfig objects based on pandas.Series
        """
        import pandas as pd
        string = '<pymicra.siteConfig> object\n{}\n----\n'.format(self.description)
        string+= pd.Series(self.__dict__).__str__()
        return string
    __repr__ = __str__


class Notation(object):
    """
    This creates an object that holds the default notation for pymicra.
    Example of usage:

    notation = notation()

    fluctuations_of_co2_conc = notation.concentration % notation.fluctuations % notation.co2

    You should be careful with the order. The last argument should not have any '%' symbols
    or you'll get a "TypeError: not all arguments converted during string formatting" message.
    """
    mean='%s_mean'
    fluctuations="%s'"
    star='%s_star'
    std='%s_std'
    mass_concentration='conc_%s'
    molar_concentration='mconc_%s'
    mass_density='rho_%s'
    molar_density='mrho_%s'
    mass_mixing_ratio='r_%s'
    molar_mixing_ratio='mr_%s'

    density=mass_density
    concentration=mass_concentration
    mixing_ratio=mass_mixing_ratio

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
    cross_spectrum = 'Cr_%s_%s'
    spectrum = 'Sp_%s'
    co_spectrum = 'Co_%s_%s'
    quadrature_spectrum = 'Qu_%s_%s'

    def __init__(self):
        """
        When initialized, calls the build method to build the full notation
        """
        self.build()

    def build(self):
        """
        This useful method builds the full notation based on the base notation
        """
        for subst in ['h2o', 'co2', 'ch4', 'o3', 'moist_air', 'dry_air']:
            for mode, comp in zip(['', '_mean' ], ['', 'self.mean %']):
                exec('self.{0}{1}_mass_density = {2} self.mass_density % self.{0}'.format(subst, mode, comp))
                exec('self.{0}{1}_molar_density = {2} self.molar_density % self.{0}'.format(subst, mode, comp))
                exec('self.{0}{1}_mass_mixing_ratio = {2} self.mass_mixing_ratio % self.{0}'.format(subst, mode, comp))
                exec('self.{0}{1}_molar_mixing_ratio = {2} self.molar_mixing_ratio % self.{0}'.format(subst, mode, comp))
                exec('self.{0}{1}_mass_concentration = {2} self.mass_concentration % self.{0}'.format(subst, mode, comp))
                exec('self.{0}{1}_molar_concentration = {2} self.molar_concentration % self.{0}'.format(subst, mode, comp))

                exec('self.{0}{1}_density = {2} self.density % self.{0}'.format(subst, mode, comp))
                exec('self.{0}{1}_mixing_ratio = {2} self.mixing_ratio % self.{0}'.format(subst, mode, comp))
                exec('self.{0}{1}_concentration = {2} self.concentration % self.{0}'.format(subst, mode, comp))

            for mode in ['fluctuations', 'star', 'std']:
                exec('self.{0}_mass_density_{1} =  self.{1} % self.{0}_mass_density'.format(subst,mode))
                exec('self.{0}_molar_density_{1} = self.{1} % self.{0}_molar_density'.format(subst,mode))
                exec('self.{0}_mass_mixing_ratio_{1} = self.{1} % self.{0}_mass_mixing_ratio'.format(subst,mode))
                exec('self.{0}_molar_mixing_ratio_{1} = self.{1} % self.{0}_molar_mixing_ratio'.format(subst,mode))
                exec('self.{0}_mass_concentration_{1} = self.{1} % self.{0}_mass_concentration'.format(subst,mode))
                exec('self.{0}_molar_concentration_{1} = self.{1} % self.{0}_molar_concentration'.format(subst,mode))
    
                exec('self.{0}_density_{1} =  self.{1} % self.{0}_density'.format(subst,mode))
                exec('self.{0}_mixing_ratio_{1} = self.{1} % self.{0}_mixing_ratio'.format(subst,mode))
                exec('self.{0}_concentration_{1} = self.{1} % self.{0}_concentration'.format(subst,mode))
    


        for subst in ['co2', 'ch4', 'o3']:
            exec('self.{0}_flux = self.flux_of % self.{0}'.format(subst))

        for subst in ['u', 'v', 'w', 'thermodyn_temp', 'virtual_temp', 'specific_humidity', 'relative_humidity', 'pressure']:
            exec('self.{0}_fluctuations = self.fluctuations % self.{0}'.format(subst))
            exec('self.{0}_mean = self.mean % self.{0}'.format(subst))
            exec('self.{0}_star = self.star % self.{0}'.format(subst))
            exec('self.{0}_std = self.std % self.{0}'.format(subst))

    def __str__(self):
        import pandas as pd
        pd.options.display.max_rows=9999
        string = pd.Series(self.__dict__).__str__()
        pd.reset_option('max_rows')
        return string
    __repr__ = __str__

