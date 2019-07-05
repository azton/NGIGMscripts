import numpy as np
import h5py
import sys
import os
from progressbar import *

"""
This takes in a .h5 file with spectrum, and does statistics to find the 
flux pdf and tau pdf. Uses pipeline info file to get correct parameters.

Usage:
    python generatePDFs.py <simname> <dump>
"""

f = open('lin_%s_RD%04d_pipeline.info'%(sys.argv[1], int(sys.argv[2])), 'r')
lines = [line for line in f]
simname = lines[0].strip()
d = int(lines[1].strip())
dim = int(lines[2].strip())
dataLoc = lines[3].strip()
nRanks = int(lines[4].strip())
nPerRank = int(lines[5].strip())
if dataLoc == 'BW':
    workdir = '/u/sciteam/wells/scratch/prodRun'
    scriptsdir = '/u/sciteam/wells/scratch/prodRun/scripts'
    dataRepo = scriptsdir
elif dataLoc == 'HD':
    workdir = '/home/azton/simulations'
    dataRepo = '/home/azton/simulations/data'#/media/azton/bigbook/projects/nextGenIGM'
    scriptsdir = '/home/azton/simulations/scripts'
elif dataLoc == 'bigbook':
    workdir = '/media/azton/bigbook/projects/nextGenIGM'
    scriptsdir = '/home/azton/simulations/scripts'
    dataRepo = workdir
elif dataLoc == 'comet':
    datarepo = '/oasis/scratch/comet/azton/temp_project/nextGenIGM'
    workdir = '/scratch/%s/%s'%(os.getenv('USER'), os.getenv('SLURM_JOB_ID'))
    scriptsdir = datarepo + '/scripts'
else:
    print('Data Location not recognized!\nUse "BW" for bluewaters, "HD" for harddrive, "bigbook" for external!')
    quit()

if dataLoc == 'comet':
        os.makedirs('%s/rayData/%s/RD%04d'%(workdir, simname, d))

#        os.system('cp %s/scripts/rayData/%s/RD%04d/rays.tar %s/rayData/%s/RD%04d/'\
#                      %(datarepo, simname, d, workdir, simname, d))
#        os.system('cd %s/rayData/%s/RD%04d; tar -xvf rays.tar'\
#                      %(workdir, simname, d))
        os.system('cp %s/scripts/rayData/%s/RD%04d/lin_spectra.tar %s/rayData/%s/RD%04d/'\
                      %(datarepo, simname, d, workdir, simname, d))
        os.system('cd %s/rayData/%s/RD%04d; tar -xvf lin_spectra.tar'\
                      %(workdir, simname, d))
widgets = ['Gathering the fluxes: ', Percentage(), ' ', Bar(marker='%', left='<',\
                right='>'),' ', SimpleProgress(), ' ', Timer()]
pbar = ProgressBar(widgets=widgets, maxval=nRanks*nPerRank+1)
pbar.start()

counts = 0
flux = []
tau = []
for i in range(nPerRank):
        try:
            filepath = '%s/rayData/%s/RD%04d/lin_spectra/spectra_%d.h5'\
                %(workdir, simname, d, i)
            if not os.path.exists(filepath): 
                continue
            f = h5py.File(filepath)
            flux.append(f['flux'][:])
            tau.append(f['tau'][:])
            counts += 1
            pbar.update(counts)
        except: continue
flux = np.array(flux).flatten()
tau = np.array(tau).flatten()
meanflux = np.average(flux)
if not os.path.exists('%s/rayData/%s/RD%04d'%(scriptsdir, simname, d)):
    os.makedirs('%s/rayData/%s/RD%04d'%(scriptsdir, simname, d))
f = open('%s/rayData/%s/RD%04d/lin_meanFluxPerPixel.txt'%(scriptsdir,simname,d), 'w')
f.write('%23.16e\n'%meanflux)
f.close()
print('Mean Flux: %0.5f'%(meanflux))
bins = np.append(0.0, np.linspace(0.025, 0.975, 20))
bins = np.append(bins, 1.0)
normPDFFlux, binsFlux = np.histogram(flux, bins=bins, density=True)
#normPDFTau, binsTau = np.histogram(tau[np.where(tau > 0.0)[0]], bins=100, density = True)
print('PDF: ', normPDFFlux)

# np.savetxt('%s/rayData/%s/RD%04d/allFlux.txt'%(scriptsdir, simname, d), flux)
np.savetxt('%s/rayData/%s/RD%04d/lin_fluxPDFbins.txt'%(scriptsdir,simname, d),binsFlux)
np.savetxt('%s/rayData/%s/RD%04d/lin_fluxPDF.txt'%(scriptsdir,simname, d),normPDFFlux)
#np.savetxt('%s/rayData/%s/RD%04d/tauPDFbins.txt'%(scriptsdir,simname, d), binsTau)
#np.savetxt('%s/rayData/%s/RD%04d/tauPDF.txt'%(scriptsdir,simname, d), normPDFTau)
    
f = open('lin_%s_RD%04d_pipeline.info'%(sys.argv[1], int(sys.argv[2])), 'a')
f.write('\nCompleted: SA_generatePDFs with %d and mean flux %23.16f\n'%(counts, meanflux))
f.close()
