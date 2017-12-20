from sys import path
path.insert(0, '/home/tomas/gits/pymicra')
import pymicra as pm
from glob import glob
from matplotlib import pyplot as plt

fconfig = pm.fileConfig('lib/tij-pr.config')
fnames = sorted(glob('ex_data/*.out'))


for fname in fnames:
    break
    data, units = pm.timeSeries(fname,fconfig,parse_dates=False)
    print(data.with_units(units))
    data[['p']].plot()
    plt.show()
# --------------------------------------------------------------------


# First, quality control
pm.util.qc_replace(fnames, fconfig,
        file_lines=36000,
        nans_test=True,
        lower_limits=dict(mrho_h2o=0, mrho_co2=0, theta_v=10),
        upper_limits=dict(theta_v=45),
        spikes_test=True,
        spikes_detrend=dict(how='linear'),
        spikes_func = myspike, 
        max_consec_spikes=10,
        max_replacement_count=360, # replacement count test
        chunk_size=1200,
        outdir='out_1',
        replaced_report='rrep.txt',
        )
fnames2 = sorted(glob('out_1/*.out'))
pm.util.qc_discard(fnames2, fconfig,
        std_limits={'u':0.03,'v':0.03,'w':0.01,'theta_v':0.02},
        std_detrend={'how':'movingmean', 'window':900},
        # MAX DIF test
        dif_limits={'u':4.0,'v':4.0,'w':1.0,'theta_v':2.0},
        maxdif_detrend={},
        maxdif_trend={'how':'movingmedian', 'window':1200},
        chunk_size=1200,
        outdir='out_2',
        summary_file='filter_summary.csv',
        full_report='rrep2.txt')


# Second, get fluxes and length estimates

