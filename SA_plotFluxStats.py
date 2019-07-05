"""
Plot up the power and flux that were generated by SA_generatePDFs
and SA_fluxPower.  uses pipeline infofile to figure out where everything is.

Usage:
python -u SA_plotPowerFlux.py <simname> <output#>
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sys, os

f = open('%s_RD%04d_pipeline.info'%(sys.argv[1], int(sys.argv[2])), 'r')
lines = [line for line in f]
simname = lines[0].strip()
d = int(lines[1].strip())
dim = int(lines[2].strip())
dataLoc = sys.argv[3] #lines[3].strip()
nRanks = int(lines[4].strip())
nPerRank = int(lines[5].strip())
if dataLoc == 'BW':
    workdir = '/u/sciteam/wells/scratch/prodRun'
    scriptsdir = '/u/sciteam/wells/scratch/prodRun/scripts'
    dataRepo = scriptsdir
elif dataLoc == 'HD':
    workdir = '/home/azton/simulations'
    dataRepo = '/home/azton/simulations/scripts'#scripts/final_dataRepo'
    scriptsdir = '/home/azton/simulations/scripts'
    analysisDir = '/home/azton/simulations/analysis'
elif dataLoc == 'bigbook':
    workdir = '/media/azton/bigbook/projects/nextGenIGM'
    scriptsdir = '/home/azton/simulations/scripts'
    dataRepo = workdir
    analysisDir = '/home/azton/simulations/analysis'
elif dataLoc == 'comet':
    workdir = '/oasis/scratch/comet/azton/temp_project/nextGenIGM'
    dataRepo = workdir+'/scripts'
    analysisDir = workdir+'/analysis'
else:
    print('Data Location not recognized!\nUse "BW" for bluewaters, "HD" for harddrive, "bigbook" for external!')
    quit()
if not os.path.exists('%s/%s'%(analysisDir, simname)):
    os.makedirs(analysisDir+'/'+simname)
rdir = '%s/rayData/%s/RD%04d'%(dataRepo, simname, d)

fluxhist = np.loadtxt('%s/fluxPDF.txt'%rdir)
fluxbins = np.loadtxt('%s/fluxPDFbins.txt'%rdir)

powerhist = np.loadtxt('%s/fluxPower.txt'%rdir)
wavenums = np.loadtxt('%s/fluxPowerWaveNum.txt'%rdir)

dataFlux=[  np.arange(0, 1.01, 0.05),\
            [5.1e-1, 1.8e-1, 1.2e-1, 1.3e-1, 1.2e-1, 1.23e-1, 1.05e-1, 1.1e-1, 1.2e-1, 1.15e-1, 1.5e-1, 1.7e-1, 0.2, 0.25,0.28,0.38, 0.55, 0.8, 1.7, 5, 9 ]]
dataPower = [[9.5e-3, 5e-3, 1.5e-3, 1.5e-2, 4.5e-2, 9.5e-2, 1.5e-1 ], \
                [18.75, 16.5, 12.5, 6., 1.8,0.14, 2e-2]]

fig, ax = plt.subplots(1,2, figsize=(8,4))
ax[0].plot(fluxbins[:-1], fluxhist, label= simname)
ax[0].set_yscale('log')
ax[0].scatter(dataFlux[0], dataFlux[1], label='est. data', s=2, c = 'r')
#ax[1].set_ylim([1e-2, 5e1])
#ax[1].set_xlim([4e-4, 2.5e-1])
ax[1].loglog(wavenums, powerhist, label = simname)
ax[1].scatter(dataPower[0], dataPower[1], label='est. data', s=2, c='r')
plt.title(simname+'_RD%04d'%d)
plt.savefig('%s/%s/RD%04d_powerFlux.png'%(analysisDir, simname, d))

