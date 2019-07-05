import yt
import sys
yt.enable_parallelism()

try:
    title = sys.argv[1]#'128_amr_nf'
    d = int(sys.argv[2])
    focus = sys.argv[3]
except:
    print ('Insufficient arguments:\nUsage: python '+\
               'verify_projectionPlots.py <sim_title> <dump>'+\
               '<focused: max or none>')
    quit()
SF = False
def stars (pfilter, data):
    filter = data[(pfilter.filtered_type, 'creation_time')] > 0.0
    return filter

yt.add_particle_filter('stars', function=stars, filtered_type = 'all', \
                       requires=['creation_time'])
def _DensityCorr(field, data):
    N = data['H_p0_number_density']*data['dz']/(data['H_p0_number_density']*data['dz']).mean()
    O = data['overdensity']
    return N*O/1e10

try:
    ds = yt.load('RD00%02d/RedshiftOutput00%02d'%(d,d))
except:
    ds = yt.load('DD%04d/data%04d'%(d,d))

#ds.add_field(('gas','densityCorr'), function=_DensityCorr, units=None)
mmh = [0.81957571, 0.13369499, 0.740843]
#if SF == True: ds.add_particle_filter('stars')
for field in ['temperature','overdensity','H_p0_number_density']:
    if focus == 'max':
        if field != "H_p0_number_density": 
                prj = yt.ProjectionPlot(ds, 'z', field, weight_field='density',\
                center = mmh, width = 0.2*ds.domain_width[0].to('kpc'))
        else: 
            prj=yt.ProjectionPlot(ds,'z',field, weight_field=None, method='integrate',\
                                      center = mmh, width = 0.2*ds.domain_width[0].to('kpc'))
    elif focus != 'max':
        if field != "H_p0_number_density":
            prj = yt.ProjectionPlot(ds, 'z', field, weight_field='density')
        else: 
            prj=yt.ProjectionPlot(ds,'z',field, weight_field=None, method='integrate')

    if SF == True: star = ds.all_data()['stars','particle_position']
    if field == 'H_p0_number_density':
        prj.set_zlim("H_p0_number_density",10**14, 10**22)
    if field == 'overdensity':
        prj.set_zlim('overdensity', 1e-1, 1e6)
    if field == 'densityCorr':
        prj.set_zlim('densityCorr', 10**-13, 10**2)
    if field == 'temperature':
        prj.set_zlim('temperature', 3e3, 3e5)
    if field == 'density':
        prj.set_zlim('density',1e-30, 1e-21)
    if SF == True:
        for s in star:
            prj.annotate_marker(s, coord_system='data', marker='.', plot_args={'color':'red', 's':20})
    prj.save(suffix='pdf',mpl_kwargs={"bbox_inches":"tight"})

ad = ds.all_data()
#units = dict(cell_mass='Msun')
#profile = yt.create_profile(ad, ['overdensity','temperature'],\
#                            n_bins=[128, 128], fields=['cell_mass'],\
#                            weight_field=None, units=units)
#plot = yt.PhasePlot.from_profile(profile)
#plot.set_zlim('cell_mass',10**5, 10**11)
#plot.annotate_title('%s z = %0.2f'%(title,ds.current_redshift))
#plot.save('r-t-profile%0.2f'%ds.current_redshift,suffix='pdf',mpl_kwargs={"bbox_inches":"tight"})


units = dict(cell_mass='Msun')

#profile = yt.create_profile(ad, ['H_p0_fraction','overdensity'],\
#                            n_bins=[256, 256], fields=['cell_mass'],\
#                            weight_field=None, units=units)
#plot = yt.PhasePlot.from_profile(profile)
plot = yt.PhasePlot(ad,'overdensity','H_fraction', ['cell_mass'], weight_field=None)
plot.set_unit('cell_mass','Msun')
plot.set_zlim('cell_mass',10**3, 10**12)
plot.annotate_title(title+' z = %0.2f'%ds.current_redshift)
plot.save('h-r-profile%0.2f'%ds.current_redshift,suffix='pdf',mpl_kwargs={"bbox_inches":"tight"})
