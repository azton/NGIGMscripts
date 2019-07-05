"""
Although anything but quick, this will iterate all data dumps in the current directory, 
then ouput a series of projections and rho-T phase diagrams.
"""
import yt
import numpy as np
import os
import matplotlib.pyplot as plt
from yt.utilities.cosmology import Cosmology as co


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
for d in range (150):
    if DDs == 1:
        fh = 'DD0%03d/data0%03d'%(d,d)
    elif RDs == 1:
        fh = 'RD0%03d/RedshiftOutput0%03d'%(d,d)
    if os.path.exists(fh) == True:
        try:
            ds = yt.load(fh)
        except:
            print ('%s failed to load.\n'%fh)
            continue
    else: continue
    ds.add_particle_filter('stars')
    z = ds.current_redshift
    ad = ds.all_data()
    vol = sum(ad['cell_volume'].in_units('Mpccm**3'))
    print ('\n\nThere are %d stars at redshift %0.3f'%(len(ad['stars','particle_mass']),ds.current_redshift))
    print ('total star mass: %0.3e'%sum(ad['stars','particle_mass']).in_units('Msun'))
    print ('in total simulation volume %0.3f Mpccm^3\n\n'%vol)
    for field in ['metallicity','number_density','temperature', 'overdensity', 'H_p0_number_density']:
        if field !='H_p0_number_density': prj = yt.ProjectionPlot(ds,'z',field,weight_field='density')
        else: prj = yt.ProjectionPlot(ds,'z',field,weight_field=None)
        fp = 'projections/%s'%field
        if os.path.exists(fp) == False: os.makedirs(fp)
        prj.save('projections/%s/%f'%(field,ds.current_redshift))
    units = dict(density='g/cm**3', cell_mass='Msun')
    profile = yt.create_profile(ad, ['density','temperature'],\
                    n_bins=[64,64], fields=['cell_mass'],
                    weight_field =None, units=units)
    plot = yt.PhasePlot.from_profile(profile)
    plot.annotate_text([0.1,0.1], 'H+G',  coord_system='axis')
    fp = 'phase/'
    if os.path.exists(fp)==False: os.makedirs(fp)
    plot.save('phase/%0.3f_profile'%ds.current_redshift)
    theoSFR = 0.015 * (1.+z)**(2.7) / (1.+((1.+z)/2.9)**5.6)

    
    
    
    theoSFRD.append(theoSFR)
    stars.append(len(ad['stars','particle_mass'])/vol \
                        if len(ad['stars','particle_mass']) > 0 else 0)
    starmass.append(np.log10(sum(ad['particle_mass'][ad['creation_time'] > 0].in_units('Msun'))\
                    /vol))
    meanZ.append(ad['metallicity'].mean())
    redshifts.append(ds.current_redshift)
fig, ax = plt.subplots(1,3, figsize=(20,10))
ax[0].plot(redshifts, stars, label = '$N_{stars}Mpccm^{-3}$')
ax[0].legend()
ax[1].plot(redshifts, np.log10(meanZ), label = '$\log_{10}Z/Z_\odot$')
ax[1].legend()
ax[2].plot(redshifts, starmass, label = '$\log_{10}M_{stars} (M_\odot)Mpccm^{-3}$')
ax[2].legend()
plt.savefig('StarStats.png')


"""
calculate real sfr density from last available dump
"""
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
plt.savefig('SFRD.png')



