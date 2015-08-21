#!/usr/bin/python
import pymicra
#import pandas as pd
from sys import path as spath
spath.append('/home/tomas/Dropbox/sbin')
from plotlab import m2dplot
from matplotlib import pyplot as plt
from sys import exit

fname='/home/tomas/Dropbox/atto/develop/testdata/CSV_9538.fluxo_42m_2012_106_1200.dat'
atto42=pymicra.dataloggerConf(varNames=['h','%Y','%j','%H%M','u','v','w'],
    description="Atto 42m",
    units={'u':'m/s', 'v':'m/s', 'w':'m/s', 'h':'m'}
)
data=pymicra.timeSeries(fname, atto42)

data=data[['u','v','w']]
print data
spec=pymicra.spectrum(data, 'u')
#spec.plot()
#plt.show()

data['detrended']=pymicra.algs.fitByDate(data['u'], degree=2)
data.plot()
plt.show()

exit()
#print data[-3:]
tup=( [spec['frequencies'], spec['spectrum'] ], )

m2dplot(tup, 'spec.png', show=True, subplot_kwargs={'xscale':'log', 'yscale':'log'})

