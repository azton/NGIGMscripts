"""
    Make projections showing most massive halo and most massive prog.
    make plots showing the evolution of fstar, H density and H fraction
    over time.
    Call:
    python3 massiveHaloTimeSeries.py /path/to/simulation/directory <final redshift dump number>
"""
import yt
import ytree
import sys
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import os

def mass_add (field, data):
    return data['mvir'].to('Msun/h')
def stars (pfilter, data):
    filter = data[(pfilter.filtered_type, 'creation_time')] > 0.0
    return filter

yt.add_particle_filter('stars', function=stars, filtered_type = 'all', \
                       requires=['creation_time'])



final = int(sys.argv[2])
simname = sys.argv[1]
#data = "/media/azton/bigbook/projects/nextGenIGM/%s"%simname
data = '/home/azton/simulations/%s'%simname
results = "/home/azton/simulations/analysis/%s"%simname


tree = ytree.load('%s/rockstar_halos/trees/tree_0_0_0.dat'%data)
tree.add_derived_field('mass',mass_add, units='Msun/h')
mostMassive = tree[np.where(tree['mvir'] == max(tree['mvir']))[0][0]]
zlist = mostMassive['prog','redshift']

mmpInfo = defaultdict(list)

fstar = []
HI_fraction = []
H_fraction = []
star_mass = []
for z in zlist:
    print ('Working z=%0.3f'%z)
    ds = 0.0
    for d in range(100):
        DDfh = '%s/rockstar_halos/halos_DD%04d.0.bin'%(data, d)
        RDfh = '%s/rockstar_halos/halos_RD%04d.0.bin'%(data, d)
        if os.path.exists(DDfh):
            rs = yt.load(DDfh)
            print (rs.current_redshift)
            if not np.isclose(rs.current_redshift, z, atol=1e-3):
                continue
            ds = yt.load('%s/DD%04d/data%04d'%(data, d, d))
        elif os.path.exists(RDfh):
            rs = yt.load(RDfh)
            print(rs.current_redshift)
            if not np.isclose(rs.current_redshift, z, atol = 1e-3):
                continue 
            ds = yt.load('%s/RD%04d/RedshiftOutput%04d'%(data, d, d))
    if ds == 0.0: 
        print ('DS for %0.3f not found!\n\n'%z)
        continue
    ad = rs.all_data()
    ds.add_particle_filter('stars')
    haloid = mostMassive['prog','Orig_halo_ID']\
                [mostMassive['prog','redshift'] == z]
    if haloid not in ad['particle_identifier']: continue
    index = np.where(ad['particle_identifier']==haloid)[0][0]
    print ('Most massive: \nMvir = %0.3e Msun '%ad['particle_mass'][index].to('Msun'))
    center = ad['particle_position'][index].to('unitary')
    radius = ad['virial_radius'][index].to('unitary')
    sp = ds.sphere(center, radius)
    width = 12.0*radius
    prj = yt.ProjectionPlot(ds,'z','temperature',weight_field='density',\
                            width=width, center=center)
    prj.annotate_sphere(center, radius=radius, coord_system='data')
    if os.path.exists('%s/mmp_prj'%results) == False:
        os.makedirs('%s/mmp_prj'%results)
    prj.save('%s/mmp_prj/%0.2f'%(results, ds.current_redshift))
    prj = yt.ProjectionPlot(ds,'z','density',weight_field='density',\
                            width=width, center=center)
    prj.annotate_sphere(center, radius=radius, coord_system='data')
    prj.save('%s/mmp_prj/%0.2f'%(results, ds.current_redshift))
    prj = yt.ProjectionPlot(ds,'z','metallicity',weight_field='density',\
                            width=width, center=center)
    prj.annotate_sphere(center, radius=radius, coord_system='data')
    prj.save('%s/mmp_prj/%0.2f'%(results, ds.current_redshift))
    prj = yt.ProjectionPlot(ds,'z','H_fraction',weight_field='density',\
                            width=width, center=center)
    prj.set_cmap(field='H_fraction', cmap='hot')
    prj.annotate_sphere(center, radius=radius, coord_system='data')
    prj.save('%s/mmp_prj/%0.2f'%(results, ds.current_redshift))
    del(prj)
    HI_fraction.append(sp['gas','H_density'].mean().to('g/cm**3'))
    fstar.append(sp['stars','particle_mass'].sum()\
                /(sp['all','particle_mass'].sum() \
                    + sp['gas','cell_mass'].sum()))
    star_mass.append(sp['stars','particle_mass'].to('Msun').sum())
    H_fraction.append(sp['gas','H_fraction'].mean())
    del(ds)
    del(rs)
    del(ad)
    del(sp)
fstar = np.array(fstar)
fstar[fstar==0] = 1
star_mass = np.array(star_mass)
H_fraction = np.array(H_fraction)
HI_fraction = np.array(HI_fraction)
zlist = zlist[:len(H_fraction)]
fig, ax = plt.subplots(2,2,figsize=(10,10))
ax[0][0].plot(zlist, fstar, label='$log f^*$')
ax[0][1].plot(zlist, np.log10(star_mass), label = "$logM_*$")
ax[1][0].plot(zlist, np.log10(HI_fraction), label = 'Neutral density')
ax[1][1].plot(zlist, np.log10(H_fraction), label = 'Neutral Fraction')
plt.title('f*, M*, H density, H fraction')
plt.savefig('%s/hyd-star-tracking.png'%results)