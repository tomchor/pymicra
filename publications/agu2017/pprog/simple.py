#from sys import path
#path.insert(0, '/home/tomas/repos/pymicra')
import pymicra as pm
import pandas as pd
from glob import glob
from matplotlib import pyplot as plt

fconfig = pm.fileConfig('tij_pr_qc.config')
fnames = sorted(glob('mydata/*.out'))



# First, quality control
pm.util.qc_replace(fnames, fconfig,
    file_lines=36000,
    lower_limits=dict(theta_v=10, mrho_h2o=0, mrho_co2=0),
    upper_limits=dict(theta_v=45),
    spikes_test=True,
    max_replacement_count=360, # replacement count test
    chunk_size=1200,
    outdir='out1',
    replaced_report='rrep.txt')

fnames2 = sorted(glob('out1/*.out'))
pm.util.qc_discard(fnames2, fconfig,
    std_limits = dict(u=0.03, v=0.03, w=0.01, theta_v=0.02),
    dif_limits = dict(u=4.0, v=4.0, w=1.0, theta_v=2.0),
    chunk_size=1200,
    outdir='out2',
    summary_file='discard_summary.csv',
    full_report='frep.txt')


# Second, get fluxes and length estimates
fconfig = pm.fileConfig('tij_pr.config')
sconfig = pm.siteConfig('tij_pr.site')
fnames = sorted(glob('out2/*out'))
allresults=[]
for fname in fnames:
    data, units = pm.timeSeries(fname, fconfig, parse_dates=False)
    data = pm.micro.preProcess(data, units, solutes=['co2'])
    fulldata = data.detrend(units=units, ignore=['p'], join_data=True)

    results = pm.micro.eddyCovariance(fulldata, units, site_config=sconfig, solutes=['co2'])
    allresults.append(results)

allresults = pd.concat(allresults, ignore_index=True)
print(allresults.with_units(units))


