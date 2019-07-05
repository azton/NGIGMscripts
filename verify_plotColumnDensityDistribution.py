import matplotlib.pyplot as plt
import numpy as np


wd = '/home/azton/simulations/analysis'
simname = '128_uni'

CD = np.loadtxt('%s/%s/columnDensities.txt'%(wd,simname))
for i in sorted(range(len(CD)),reverse=True):
    if CD[i] == 0.0:
        CD = np.delete(CD, i)
nsamp = len(CD)
print("%d good density lines"%nsamp)
np.savetxt('%s/%s/columnDensities_%d.txt'%(wd,simname, nsamp), CD)
plt.figure()
hist, bins = np.histogram(np.log10(CD), bins=100)
hist = [i/nsamp for i in hist ]
plt.plot(bins[:-1], hist)
plt.xlabel('$\log N_c(1/cm^{-2})$')
plt.ylabel('$N/N_{samples}$')
plt.yscale('log')
plt.title('Columns for %d samples'%nsamp)
plt.savefig('%s/%s/CD_%d.png'%(wd,simname, nsamp))
