"""
    This generates a bunch of star statistics for the simulations and
    outputs them to the appropriate analysis directory.
    Only works if simulation directories and analysis directory are in the same 
    Call from command line with the simulation name you want to analyze:
    python3 quickcheck.py /path/to/simulation/directory

"""

import yt
import numpy as np
import os
import matplotlib.pyplot as plt
from yt.utilities.cosmology import Cosmology as co
import sys
yt.enable_parallelism()
simname = sys.argv[1]
final = int(sys.argv[2])
def stars (pfilter, data):
    filter = data[(pfilter.filtered_type, 'creation_time')] > 0.0
    return filter

yt.add_particle_filter('stars', function=stars, filtered_type = 'all', \
                       requires=['creation_time'])


DDs = 0
RDs = 1
stars = []
starmass = []
meanZ = []
redshifts = []
theoSFRD = []
for d in range (final+1):
    if DDs == 1:
        fh = '../%s/DD0%03d/data0%03d'%(simname,d,d)
    elif RDs == 1:
        fh = '../%s/RD0%03d/RedshiftOutput0%03d'%(simname,d,d)
    if os.path.exists(fh) == True:
        try:
            ds = yt.load(fh)
        except:
            print ('%s failed to load.\n'%fh)
            continue
    else: continue
    ds.add_particle_filter('stars')
    z = ds.current_redshift
    if os.path.exists('../analysis/%s/phase/%0.3f_profile'\
                %(simname,ds.current_redshift)):
        continue
    ad = ds.all_data()
    vol = ad['cell_volume'].sum().in_units('Mpccm**3')
    if yt.is_root():
        print ('\n\nThere are %d stars at redshift %0.3f'%(len(ad['stars','particle_mass']),ds.current_redshift))
        print ('total star mass: %0.3e'%(ad['stars','particle_mass'].sum().in_units('Msun')))
        print ('in total simulation volume %0.3f Mpccm^3\n\n'%vol)
    for field in ['metallicity','number_density','temperature', 'overdensity']:
        prj = yt.ProjectionPlot(ds,'z',field,weight_field='density')
        fp = '../analysis/%s/projections/%s'%(simname,field)
        if os.path.exists(fp) == False: os.makedirs(fp)
        if field == 'temperature':
            prj.set_zlim('temperature',10**3.5, 10**5)
        prj.save('../analysis/%s/projections/%s/%f'\
                %(simname,field,ds.current_redshift))
        del(prj)
    units = dict(density='g/cm**3', cell_mass='Msun')
    profile = yt.create_profile(ad, ['density','temperature'],\
                    n_bins=[64,64], fields=['cell_mass'],
                    weight_field =None, units=units)
    plot = yt.PhasePlot.from_profile(profile)
    plot.annotate_text([0.1,0.1], 'H+G',  coord_system='axis')
    fp = '../analysis/%s/phase/'%simname
    if yt.is_root():
        if os.path.exists(fp)==False: os.makedirs(fp)
    plot.save('../analysis/%s/phase/%0.3f_profile'\
            %(simname,ds.current_redshift))
    del(plot)
    theoSFR = 0.015 * (1.+z)**(2.7) / (1.+((1.+z)/2.9)**5.6)

    
    
    
    theoSFRD.append(theoSFR)
    stars.append(len(ad['stars','particle_mass'])/vol \
                        if len(ad['stars','particle_mass']) > 0 else 0)
    starmass.append(np.log10(sum(ad['particle_mass'][ad['creation_time'] > 0].in_units('Msun'))\
                    /vol))
    meanZ.append(ad['metallicity'].mean())
    redshifts.append(ds.current_redshift)
    del(ds)
    del(ad)
fig, ax = plt.subplots(1,3, figsize=(20,10))
ax[0].plot(redshifts, stars, label = '$N_{stars}Mpccm^{-3}$')
ax[0].legend()
ax[1].plot(redshifts, np.log10(meanZ), label = '$\log_{10}Z/Z_\odot$')
ax[1].legend()
ax[2].plot(redshifts, starmass, label = '$\log_{10}M_{stars} (M_\odot)Mpccm^{-3}$')
ax[2].legend()
plt.savefig('../analysis/%s/StarStats.png'%simname)
del(fig)
del(ax)
del(stars)
del(meanZ)
del(starmass)

"""
calculate real sfr density from last available dump
"""
ds = yt.load('../%s/RD0%03d/RedshiftOutput0%03d'%(simname, final, final))
ds.add_particle_filter('stars')
ad = ds.all_data()
masses = ad['stars','particle_mass'].in_units('Msun')
formation_time = ad['stars','creation_time'].in_units('yr')
time_range = [0, float(ds.current_time.in_units('yr'))]
n_bins=100
hist, bins = np.histogram(formation_time, bins=n_bins, range=time_range)
inds = np.digitize(formation_time, bins=bins)
time = (bins[:-1]+bins[1:])/2
sfr = np.array([masses[inds==j+1].sum()/(bins[j+1]-bins[j])\
                for j in range(len(time))])
sfr[sfr == 0] = np.nan
sfrd = sfr/vol
time = [co().z_from_t(my_time=ds.quan(t, 'yr')) for t in time]
fig, ax = plt.subplots(1,1, figsize=(10,10))
ax.plot(redshifts, theoSFRD, label="Analytic", c='k',linestyle='--')
ax.plot(time, sfrd, label='SFRD')
ax.set_yscale('log')
ax.set_xlabel('$z$')
ax.set_ylabel('SFRD $M_\odot Mpccm^{-3} yr^{-1}$')
plt.savefig('../analysis/%s/SFRD.png'%simname)



