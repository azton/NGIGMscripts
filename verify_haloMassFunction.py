                                                
"""                                                                                        
    Use Rockstar catalog to generate a                                                     
    halo mass function.                                                                    
    call with                                                                              
    python haloMassFunction.py <simname> <final output #>                                  
"""

import yt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
from yt.mods import *
from yt.analysis_modules.halo_mass_function.api import *

simname = sys.argv[1]
d = int(sys.argv[2])
#workdir = '/media/azton/bigbook/projects/nextGenIGM'
workdir = '/home/azton/simulations'
resultsdir = '/home/azton/simulations/analysis'
rs = yt.load('%s/%s/rockstar_halos/halos_RD%04d.0.bin'%(workdir, simname, d))
ds = yt.load('%s/%s/RD%04d/RedshiftOutput%04d'%(workdir, simname, d,d))
hmf = HaloMassFcn(halos_ds = rs, simulation_ds = ds, log_mass_min=6, log_mass_max=12,\
                    fitting_function=1)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.loglog(hmf.masses_sim.to('Msun/h'), hmf.n_cumulative_sim.to("Mpc**-3*h**3"), linestyle='--', label='data')
ax.loglog(hmf.masses_analytic.to('Msun/h'), hmf.n_cumulative_analytic.to("Mpc**-3*h**3"), \
            linestyle=":", label='Warren')
plt.savefig('%s/%s/HMF_%0.2f_RD%04d.png'%(resultsdir, simname, ds.current_redshift, d))
print(hmf.masses_sim)
print(hmf.n_cumulative_sim)
ad = rs.all_data()
masses = np.array([float(i) for i in ad['particle_mass'].to('Msun/h')])
hist, bins =np.histogram(np.log10(masses), bins=25)
volume = (abs(ds.domain_right_edge-ds.domain_left_edge)[0]).to('Mpc/h')**3
print (volume)
ax.plot(10**bins[:-1], hist/volume/(bins[1]-bins[0]))
plt.savefig('%s/%s/HMF_%0.2f_RD%04d.png'%(resultsdir, simname, ds.current_redshift, d))



