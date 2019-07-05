import yt
from scipy.optimize import curve_fit
import sys
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from progressbar import *
simname = sys.argv[1]
d = int(sys.argv[2])

def t0(t, gamma, delta):
    return t*delta**(gamma-1)

def manual_mean_density(data):
    mean = 0
    cnt = 1
    for r in data['density']:
        mean += (r-mean)/cnt
        cnt += 1
        bar.update(cnt)
    bar.finish()
    return mean
def _delta(field, data):
    dmean = data['density'].mean()
    return data['density']/dmean
try:
    ds = yt.load('RD%04d/RedshiftOutput%04d'%(d,d))
except:
    ds = yt.load('DD%04d/data%04d'%(d,d))
print('adding field Delta')
ds.add_field(('gas','Delta'), function=_delta, units=None)
print('Added field!')
ad = ds.all_data()
print('making filter/dataset')
low = 0.95
hi = 1.05
filter = np.logical_and(ad['Delta'] > low, ad['Delta'] < hi)
print('made filter')
xdata = ad['Delta'][filter]
ydata = ad['temperature'][filter]
print('input made!')
print('fitting')
popt, pcov = curve_fit(t0, xdata, ydata, p0=[6e3, 1.0005])
print('fitted with')
print (popt)
xpred = np.arange(0.95, 1.05, 0.001)
ypred = [t0(popt[0], popt[1], s) for s in xpred]
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(xdata, ydata, s=2, alpha=0.4)
ax.plot(xpred, ypred, linewidth=2, c='r', \
        label = '$T_0 =$%0.2f\n$\gamma = %0.6e'%(popt[0], popt[1]))
ax.set_title('%s $T_o$ fit z = %0.2f'%(simname, ds.current_redshift))
ax.set_ylabel('$T$[K]')
ax.set_xlabel('$\Delta$')
ax.set_yscale('log')
plt.text(0,0,'$T = %0.3f \Delta ^{%0.6f-1}$')
plt.savefig('r-t_scatter_t0fit.png')
f = open('%s_t0.info'%simname,'w')
f.write('t0\t\tgamma\n%0.3f\t\t%23.16e'%(popt[0], popt[1]))
f.close()
