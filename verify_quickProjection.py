import yt
import sys

yt.enable_parallelism("MPI_COMM_WORLD")

weight = False
SF = False
outputDir ='/media/azton/bigbook/projects/nextGenIGM/1024_uni+BG'
def stars (pfilter, data):
    filter = data[(pfilter.filtered_type, 'creation_time')] > 0.0
    return filter

yt.add_particle_filter('stars', function=stars, filtered_type = 'all', \
                       requires=['creation_time'])
def _HIColumnDensity(field, data):
    return data['H_p0_number_density']*data['dz']


d = int(sys.argv[1])
field = sys.argv[2]
ds = yt.load('%s/RD00%02d/RedshiftOutput00%02d'%(outputDir,d,d))
ds.add_field(('gas','HI_column_density'), function=_HIColumnDensity, units='1/cm**2')
if SF == True: ds.add_particle_filter('stars')
if weight == False: 
    prj=yt.ProjectionPlot(ds,'z',field, weight_field=None, method='sum')
else: 
    prj = yt.ProjectionPlot(ds, 'z', field, weight_field='density')
if SF == True: star = ds.all_data()['stars','particle_position']
if field == 'temperature':
    prj.set_zlim('temperature',10**3.5, 10**5)
if field == 'HI_column_density':
    prj.set_zlim("HI_column_density",10**14, 10**24)
if SF == True:
    for s in star:
        prj.annotate_marker(s, coord_system='data', marker='.', plot_args={'color':'red', 's':20})
prj.annotate_title('%s'%simname)
prj.save()
ad = ds.all_data()
units = dict(density='g/cm**3', cell_mass='Msun')
profile = yt.create_profile(ad, ['density','temperature'],\
                            n_bins=[64,64], fields=['cell_mass'],\
                            weight_field=None, units=units)
plot = yt.PhasePlot.from_profile(profile)
plot.save('profile%d.png'%ds.current_redshift)
