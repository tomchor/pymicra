"""
Defines classes that are the basis of Pymicra
"""
from __future__ import print_function
from . import decorators as _decors


class fileConfig(object):
    """
    This class defines a specific configuration of a data file

    Parameters
    ----------
    from_file: str
        path of .cfg file (configuration file) to read from. This will ignore all other
        keywords.
    variables: list of strings or dict
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
        used to assemble the date. If the file is tabular-separated then this should be "whitespace".
    header_lines: int or list
        up to which line of the file is a header. See pandas.read_csv header option.
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
    varNames: DEPRECATED
        use variables now.
    """

    @_decors.autoassign
    def __init__(self,
            from_file=None,
            variables=None,
            date_cols=None,
            frequency=None,
            date_connector='-',
            columns_separator=None,
            header_lines=None,
            units=None,
            filename_format=None,
            skiprows=None,
            varNames=None,
            description='Type help(fileConfig) to read intructions', **read_csv_kwargs):
        """
        Reads the arguments, transforms them into attributes of the object and calcultes some
        other attributes.
        """
        from . import algs

        #-------------------------
        if from_file:
            from . import io
            dlconf = io.read_fileConfig(from_file)
            self.__dict__.update(dlconf.__dict__)
            return
        #-------------------------

        #-------------------------
        # Makes sure that units is a dictionary type
        if units is not None:
            if not isinstance(units, dict):
                raise TypeError('units should be a dictionary. Ex.: {"u" : "m/s", "v" : "m/s", "theta" : "K" }')
            else:
                self.units = algs.parseUnits(units)
        #-------------------------

        #------------------
        # Recuperates deprecated option varNames
        if (not variables) and varNames:
            self.variables = varNames
        #------------------

        self.get_date_cols()
        self._check_consistency()
        self._unite_header()


    def get_date_cols(self):
        """
        Guesses what are the columns that contain the dates by searching
        for percentage signs in them
        """
        if self.date_cols:
            self.date_col_names = [ self.variables[ idx ] for idx in self.date_cols ]
        else:
            date_cols = { k : it for (k, it) in self.variables.iteritems() if '%' in it }
            self.date_col_names = date_cols.values()
            self.date_cols = date_cols.keys()


    def _check_consistency(self):
        """
        Checks consistency of fileConfig
        Currently only checks every key of self.units dictionary agaisnt values of variables dict.
        """
        for key in self.units.keys():
            if key not in self.variables.values():
                print('fileConfig WARNING!:\n    {} is defined in "units" but not defined in "variables"!'.format(key))


    def _unite_header(self):
        """
        Checks if both skiprows and header_lines were given. Since header lines
        have to be skipped, we unite header_lines and skiprows into skiprows
        """
        if self.skiprows==None:
            self.skiprows = []
        elif isinstance(self.skiprows, int):
            self.skiprows = range(self.skiprows)

        if self.header_lines==None:
            self.header_lines = []
        elif isinstance(self.header_lines, int):
            self.header_lines = range(self.header_lines)

        self.skiprows = list(set(self.skiprows) | set(self.header_lines))


    def __str__(self):
        return '<pymicra.fileConfig>\n{}'.format(self.description)
    __repr__ = __str__


class siteConfig(object):
    """
    Keeps the configurations and constants of an experiment. (such as height of instruments,
    location, canopy height and etc)

    Check help(pm.siteConfig.__init__) for other parameters

    Parameters
    ----------
    from_file: str
        path to .site file which contais other keywords
    """
    @_decors.autoassign
    def __init__(self, from_file=None,
             instruments_height=None, measurement_height=None, canopy_height=None,
             displacement_height=None, roughness_length=None,
             latitude=None, longitude=None, altitude=None, description=None, **kwargs):
        """
        Parameters
        ----------
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


    def __str__(self):
        """
        Creates nice representation for printing siteConfig objects based on pandas.Series
        """
        import pandas as pd
        string = '<pymicra.siteConfig> object\n{}\n----\n'.format(self.description)
        string+= pd.Series(self.__dict__).drop('description').__str__()
        return string
    __repr__ = __str__



class Notation(object):
    """
    Holds the notation used in every function of pymicra except when told otherwise.
    """
    @_decors.autoassign
    def __init__(self, solutes=['h2o', 'co2', 'ch4', 'o3', 'n2o'], **kwargs):
        """
        When initialized, calls the build method to build the full notation

        Parameters
        ----------
        solutes: list
            list of solutes for which to build the notation (still to be implemented)

        Returns
        -------
        pymicra.Notation
            default Notation object
        """

        self.mean='%s_mean'
        self.fluctuations="%s'"
        self.star='%s_star'
        self.std='%s_std'
        self.mass_concentration='conc_%s'
        self.molar_concentration='mconc_%s'
        self.mass_density='rho_%s'
        self.molar_density='mrho_%s'
        self.mass_mixing_ratio='r_%s'
        self.molar_mixing_ratio='mr_%s'

        self.flux_of = 'F_%s'

        self.cross_spectrum = 'X_%s_%s'
        self.spectrum = 'sp_%s'
        self.cospectrum = self.cross_spectrum.replace('X', 'co')
        self.quadrature = self.cross_spectrum.replace('X', 'qu')

        self.u='u'
        self.v='v'
        self.w='w'
        self.rotated_u=self.u
        self.rotated_v=self.u
        self.rotated_w=self.u
    
        self.thermodynamic_temperature='theta'
        self.virtual_temperature='theta_v'
        self.sonic_temperature='theta_s'
        self.potential_temperature='theta_p'
        self.pressure='p'
        self.relative_humidity='rh'
        self.specific_humidity='q'
    
        self.h2o='h2o'
        self.co2='co2'
        self.ch4='ch4'
        self.o3 = 'o3'
        self.n2o='n2o'

        self.moist_air = 'air'
        self.dry_air   = 'dry'
    
        self.stability_parameter = 'zeta'
        self.obukhov_length = 'Lo'
        self.momentum_flux = 'tau'
        self.sensible_heat_flux = 'H'
        self.virtual_sensible_heat_flux = 'Hv'
        self.water_vapor_flux = 'E'
        self.latent_heat_flux = 'LE'
    
        self.build()

    def build(self, from_level=0):
        """
        This useful method builds the full notation based on the base notation.

        Given notation for means, fluctuations, and etc, along with names of variables, this 
        method builds the notation for mean h2o concentration, virtual temperature fluctuations
        and so on.

        Parameters
        ----------
        self: pymicra.Notation
            notation to be built
        from_level: int
            level from which to build. If 0, build everything from scratch and higher notations
            will be overwritten. If 1, skip one step in building process. Still to be implemented!

        Returns
        -------
        pymicra.Notation
            Notation object with built notation
        """

        #-------------
        # Iterates over every "air" variable
        for subst in ['h2o', 'co2', 'ch4', 'o3', 'n2o', 'moist_air', 'dry_air']:
            for mode, comp in zip(['', 'mean_' ], ['', 'self.mean %']):
                exec('self.{1}{0}_mass_density = {2} self.mass_density % self.{0}'.format(subst, mode, comp))
                exec('self.{1}{0}_molar_density = {2} self.molar_density % self.{0}'.format(subst, mode, comp))
                exec('self.{1}{0}_mass_mixing_ratio = {2} self.mass_mixing_ratio % self.{0}'.format(subst, mode, comp))
                exec('self.{1}{0}_molar_mixing_ratio = {2} self.molar_mixing_ratio % self.{0}'.format(subst, mode, comp))
                exec('self.{1}{0}_mass_concentration = {2} self.mass_concentration % self.{0}'.format(subst, mode, comp))
                exec('self.{1}{0}_molar_concentration = {2} self.molar_concentration % self.{0}'.format(subst, mode, comp))


            for mode in ['fluctuations', 'star', 'std']:
                exec('self.{0}_mass_density_{1} = self.{1} % self.{0}_mass_density'.format(subst, mode))
                exec('self.{0}_molar_density_{1} = self.{1} % self.{0}_molar_density'.format(subst, mode))
                exec('self.{0}_mass_mixing_ratio_{1} = self.{1} % self.{0}_mass_mixing_ratio'.format(subst, mode))
                exec('self.{0}_molar_mixing_ratio_{1} = self.{1} % self.{0}_molar_mixing_ratio'.format(subst, mode))
                exec('self.{0}_mass_concentration_{1} = self.{1} % self.{0}_mass_concentration'.format(subst, mode))
                exec('self.{0}_molar_concentration_{1} = self.{1} % self.{0}_molar_concentration'.format(subst, mode))
        #-------------
    


        dic = self.__dict__
        for subst in ['co2', 'ch4', 'o3', 'n2o']:
            dic['%s_flux' % subst] = self.flux_of % subst

        substances = ['u', 'v', 'w', 'thermodynamic_temperature', 'virtual_temperature', 
                    'sonic_temperature', 'potential_temperature', 'specific_humidity', 'relative_humidity', 'pressure']
        for subst in substances:
            exec('self.{0}_fluctuations = self.fluctuations % self.{0}'.format(subst))
            exec('self.mean_{0} = self.mean % self.{0}'.format(subst))
            exec('self.{0}_star = self.star % self.{0}'.format(subst))
            exec('self.{0}_std = self.std % self.{0}'.format(subst))

        #-------------
        # Here we apply the aliases on the attributes to make it easier to write some long names
        self._apply_aliases()
        #-------------


    def _apply_aliases(self):
        """
        Applies short names using aliases
        """
        aliases = { 'thermodynamic':'thermodyn',
                    'temperature':'temp',
                    'parameter':'variable' }

        dic = self.__dict__
        for key, val in dic.items():
            short=key
            for lon in aliases.keys():
                short = short.replace(lon, aliases[lon])
            if key!=short:
                dic.update({short:val})


    def __str__(self):
        """
        Currently displays Notations as a very long pandas.Series
        """
        import pandas as pd
        pd.options.display.max_rows=9999
        string = pd.Series(self.__dict__).__str__()
        pd.reset_option('max_rows')
        return string
    __repr__ = __str__



import pandas as _pd
class _myData(object):
    """
    Attempt to create a myData object
    """
    def __init__(self, df, dic):
        self.df = df.copy()
        self.dic = dic
#        self.df = df

    #def as_df(self):
    #    import pandas as pd
    #    return self.copy()

    def __repr__(self):
        return self.df.__repr__()

    #def __str__(self):
    #    return _pd.DataFrame.__str__(self)

    #def __div__(self):
    #    import pandas as pd
       # return pd.DataFrame(self)
