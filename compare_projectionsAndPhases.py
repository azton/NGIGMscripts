import yt
import numpy as np 
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.colorbar as colorbar
from yt.visualization.base_plot_types import get_multi_plot
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.pyplot as plt
yt.enable_parallelism()
def _Nc(field, data):
        return data['H_p0_number_density']
def _T(field, data):
        return data['temperature']
def _oden(field, data):
        return data['overdensity']
def _logden(field, data):
        return np.log10(data['density'])
def _logT(field,data):
        return np.log10(data['temperature'])
workdir = "/oasis/scratch/comet/azton/temp_project/nextGenIGM"
sims = ["1024_uni+BG","1024_uni_noshield"]
simtitle = ['shielded', 'Unshielded']
dumps = [20,20]
fields = ['overdensity','T','N']
titles = ['Overdensity ($\rho/\bar\rho$)','Temperature [K]','HI Column Density [cm$^-1$]']
ranges = [[1e-1, 1e5], [5e3, 5e5], [1e13, 1e18]]
cmaps = ['bds_highcontrast','hot','Rainbow + black']
oreint = 'horizontal'
mmh = [0.81957571, 0.13369499, 0.740843]
fig = plt.figure(figsize=(20,20))
grid = AxesGrid(fig, (0.05, 0.05, 0.9, 0.9),\
                nrows_ncols=(1,2),\
                axes_pad=0.4,\
                label_mode='L',\
                share_all=True,\
                cbar_location="right",\
                cbar_mode='edge',\
                cbar_size='5%',\
                cbar_pad='0%',\
                direction='column',\
)
sets = [workdir+"/%s/RD%04d/RedshiftOutput%04d"%(sims[s], dumps[s], dumps[s]) for s in range(len(sims))]
cnt = 0
for i, dset in enumerate(sets):
    for j, f in enumerate(fields):
        ds = yt.load("%s"%dset)
        ds.add_field(('gas','N'), function=_Nc, units="cm**-3")
        ds.add_field(('gas','T'), function= _T, units = "K")
        ds.add_field(('gas','Delta'), function=_oden, units=None)
        if f == "N":
            prj = yt.ProjectionPlot(ds, 'z', f, \
                    weight_field=None, center = mmh, \
                    width=0.5*ds.domain_width[0].to('kpc'))
        else:
            prj = yt.ProjectionPlot(ds, 'z', f, \
                weight_field='density', center = mmh, \
                width=0.5*ds.domain_width[0].to('kpc'))
        prj.set_zlim(f, ranges[j][0], ranges[j][1])
        prj.set_cmap(f, cmaps[j])
        plot = prj.plots[f]
        plot.figure = fig
        plot.axes = grid[cnt].axes
        plot.cax = grid.cbar_axes[cnt]
        prj._setup_plots()
        cnt += 1
        del (prj)
plt.rcParams.update({"font.size":8})
plt.savefig('%s-%s_compare.pdf'%(simtitle[0], simtitle[1]),\
		    dpi=200, bbox_inches="tight", pad_inches=0.0)

## phase plots done similarly

for i, dset in enumerate(sets):
        ds = yt.load(dset)
        ds.add_field(('gas','log_Density'), function=_logden, units=None)
        ds.add_field(('gas','log_T'), function=_logT, units=None)
        ad = ds.all_data()
        plot = yt.PhasePlot(ad, 'density','temperature',['cell_mass'],\
                weight_field=None)
        profile = yt.create_profile(ad, ['density','temperature'],\
                n_bins=[128,128], fields=['cell_mass'],\
                weight_field=None, units={'cell_mass':'Msun'})
        plot = yt.PhasePlot.from_profile(profile)
        plot.set_zlim('cell_mass', 1e5, 1e10)
        plot.set_ylim(2e3, 5e6)
        plot.set_xlim(3e-30, 1e-24)
        plot.set_unit('cell_mass','Msun')
        plot.set_cmap('cell_mass',"inferno")
        print(simtitle[i])
        plot.annotate_title('%s'%simtitle[i])
        plot.set_font_size(35)
        plot.title.set_font_size(55)
        if i != 0:
                plot.set_xlabel("")
                plot.set_ylabel("")
        if i != 4:
                plot.set_colorbar_label('cell_mass','')
        plot.save('phase%s.pdf'%(sims[i]), mpl_kwargs={"dpi":200,\
                         "bbox_inches":"tight","pad_inches":0.0})
#os.system('montage -tile 4x1 -geometry +0+0 phase%s.png phase%s.png phase%s.png phase%s.png phases.png'\
#        %(sims[0],sims[1],sims[2],sims[3]))
#os.system('montage -tile 2x1 -geometry +0+0 phases.png Comp_Allplots.png ProdFig.png')
