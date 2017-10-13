import pymicra as pm
from glob import glob

config=pm.fileConfig('examples/lake.config')
files=sorted(glob('examples/ex_data/*csv'))

print(files)
pm.qcontrol(files,config, file_lines=72000, spikes_test=True, spikes_detrend=True, 
            replaced_report=None,
            spikes_func=lambda x: abs( x - x.mean() ) > 5*x.std(),
            std_stat_limits={'u':0, 'v':0},
            lower_limits={'theta_v':0.0},
            upper_limits={'theta_v':50.0},
            dif_limits={'u':2.0,'v':4.0,'w':1.0,'theta_v':2.0},
            maxdif_trend=True,
            maxdif_trend_kw={'how':'movingmedian','window':9600},
            maxdif_detrend=True,
            maxdif_detrend_kw={'how':'linear'},
            RAT=False,
            nans_test=True,
            )


