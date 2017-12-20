
# coding: utf-8
#\begin{block}{Preamble}
# In[1]:


import pymicra
import pandas as pd
from glob import glob
import matplotlib.pyplot as plt

fconfig = pymicra.fileConfig('tij_pr_qc.config')

#\end{block}

#\begin{block}{Basic quality control example}
# In[2]:


fnames = sorted(glob('mydata/*.out'))
# Prints reports on screen and writes further info to file
pymicra.util.qc_replace(fnames, fconfig,
    file_lines=36000,
    lower_limits=dict(theta_v=10, mrho_h2o=0, mrho_co2=0),
    upper_limits=dict(theta_v=45),
    spikes_test=True,
    max_replacement_count=360,
    chunk_size=1200,
    outdir='out1',
    replaced_report='rrep.txt')

fnames2 = sorted(glob('out1/*.out'))
# Prints reports on screen and writes further info to file
pymicra.util.qc_discard(fnames2, fconfig,
    std_limits = dict(u=0.03, v=0.03, w=0.01, theta_v=0.02),
    dif_limits = dict(u=4.0, v=4.0, w=1.0, theta_v=2.0),
    chunk_size=1200,
    outdir='out2',
    summary_file='discard_summary.csv',
    full_report='frep.txt')

#\end{block}

#\begin{block}{Pre-processing and calculation of fluxes example}
# In[3]:


fconfig = pymicra.fileConfig('tij_pr.config')
sconfig = pymicra.siteConfig('tij_pr.site')
fname = 'out2/20110224-1340.out'

data, units = pymicra.timeSeries(fname, fconfig, parse_dates=False)
# Prints reports showing which calculations are being done
data = pymicra.micro.preProcess(data, units, solutes=['co2'])
fulldata = data.detrend(units=units, ignore=['p'], join_data=True)

# Prints reports showing which calculations are being done
results = pymicra.micro.eddyCovariance(fulldata, units, 
                                site_config=sconfig, 
                                wpl=True, solutes=['co2'])
print(results.with_units(units))

#\end{block}

#\begin{block}{Visualization is easy}
# In[4]:


data[['u', 'v', 'w']].plot(figsize=(6,4), grid=True)
plt.tight_layout()
plt.savefig('uvw.pdf')


# In[5]:


uw_spectra = pymicra.spectra(fulldata[["u'", "w'"]], 
                        frequency=20, anti_aliasing=True)
axis=uw_spectra.binned(bins_number=100).plot(loglog=True, 
                                        grid=True, figsize=(6,4))
axis.plot(uw_spectra.index, 1e-1*uw_spectra.index**(-5/3), label='-5/3')
axis.legend(); axis.set_ylim(None, 1e4)
plt.tight_layout()
plt.savefig('spec.pdf')
plt.show()

#\end{block}
