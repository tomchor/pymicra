from __future__ import print_function
from sys import path
path.insert(0, './..')
import pymicra as pm
import pandas as pd
from glob import glob
from matplotlib import pyplot as plt

#-----------
# We read the site and file configurations
siteconf = pm.siteConfig('lake.site')
fileconf = pm.fileConfig('lake.config')
#-----------

#-----------
# Select the files we want to analyse
filelist = sorted(glob('ex_data/*'))
#-----------

allresults=[]
for fname in filelist:
    data, units = pm.timeSeries(fname, fileconf, parse_dates=True)
    data = pm.micro.preProcess(data, units, expand_temperature=True, rotation='2d',
        use_means=False, rho_air_from_theta_v=True, skip_h2o=True,
        solutes=['co2'], inplace_units=True)
    ddata = data.detrend(how='linear', units=units, ignore=['theta', 'p'])
    data = data.join(ddata)

    print(data.with_units(units))
#    results = pm.micro.eddyCovariance(data, units, site_config=siteconf, get_turbulent_scales=True, wpl=True)
#    allresults.append(results)

#fluxes = pd.concat(allresults)
#print(fluxes.with_units(units))
