import yt
import trident
import sys
yt.enable_parallelism()

workdir = '/media/azton/bigbook/projects/nextGenIGM'
simname = sys.argv[1]
final = int(sys.argv[2])
d = final
rs = '%s/%s/rockstar_halos/halos_RD%04d.0.bin'%(workdir, simname, final)
rs = yt.load(rs)
ad = rs.all_data()
ds = yt.load('%s%s/RD%04d/RedshiftOutput%04d'%(workdir, simname,d,d))
#print(ds.field_list)
trident.add_ion_fields(ds, ['C', "O", "O IV", "O VI", "C IV"], ftype='gas')

for field in ds.derived_field_list:
    print (field)

pos = ad['particle_position'][ad['particle_mass'] \
            == ad['particle_mass'].max()][0].to('unitary')
r = 4.0*ad['virial_radius'][ad['particle_mass'] \
            == ad['particle_mass'].max()][0].to('unitary')
sp = ds.sphere(pos, r)
fields = ['H_p1_density','O_p4_density','C_p0_density', 'H_density']
labels = ['HI', "OIV", "C", "H"]


## make a ray that intersects this most massive halo ##
ray = trident.make_simple_ray(ds, start_position=[pos[0], pos[1],0], end_position=[pos[0], pos[1], 1],\
                                lines='H')
width = ds.domain_width[0]
lbase = 1215.67
dz = ds.cosmology.z_from_t(ds.current_time-width/ds.quan(2.998*10**8, 'm/s'))
lmin = lbase*(z-dz+1)
lmax = lbase*(z+1)
sg = trident.SpectrumGenerator(lambda_min = lmin, lambda_max = lmax, dlambda=0.01)
sg.make_spectrum(ray, lines='H')
sg.save_spectrum('/home/azton/simulations/analysis/%s/spectrum_%0.3f.png'%(simname, ds.current_redshift)

for field in range(len(fields)):
    prj = yt.ProfilePlot(sp, "radius",fields[field], label=labels[field]\
                , weight_field='density', n_bins=100)
    prj.annotate_ray(ray, arrow=True)
    prj.save('/home/azton/simulations/analysis/%s/radialProfile'%simname)