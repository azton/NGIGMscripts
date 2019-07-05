import yt
from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit as GTF
from progressbar import *
import h5py
import sys
import os
from mpi4py import MPI
"""
Generate voigt profile fits using YT.  Run after
LR+spectra.py. Uses pipeline infofile to work
Use MPI parallelism to distribute the work.

usage:
python -u voigtFit_YT.py <simname> <output#>
"""
comm = MPI.COMM_WORLD
comm_rank = comm.Get_rank()
comm_size = comm.Get_size()
f = open('%s_RD%04d_pipeline.info'%(sys.argv[1], int(sys.argv[2])), 'r')
lines = [line for line in f]
lines = comm.bcast(lines)
simname = lines[0].strip()
d = int(lines[1].strip())
dim = int(lines[2].strip())
dataLoc = lines[3].strip()
nRanks = int(lines[4].strip())
nPerRank = int(lines[5].strip())
if dataLoc == 'BW':
    workdir = '/u/sciteam/wells/scratch/prodRun'
    scriptsdir = '/u/sciteam/wells/scratch/prodRun/scripts'
    dataRepo = '/u/sciteam/wells/scratch/prodRun/scripts'
elif dataLoc == 'HD':
    workdir = '/home/azton/simulations'
    dataRepo = '/media/azton/bigbook/projects/nextGenIGM'
    scriptsdir = '/home/azton/simulations/scripts'
elif dataLoc == 'bigbook':
    workdir = '/media/azton/bigbook/projects/nextGenIGM'
    dataRepo = '/media/azton/bigbook/projects/nextGenIGM'
    scriptsdir = '/home/azton/simulations/scripts'
elif dataLoc == 'comet':
    workdir = 'oasis/scratch/comet/azton/temp-project/nextGenIGM'
else:
    print('Data Location <%s> not recognized!'%dataLoc)
    print('Use "BW" for bluewaters, "HD" for harddrive, "bigbook" for external!')
    quit()
localRanks = range(comm_rank, nRanks, comm_size)
resultsdir = '%s/rayData/%s/RD%04d/vpfits/yt'%(dataRepo, simname, d)
if comm_rank == 0:
    if not os.path.exists(resultsdir):
        os.makedirs(resultsdir)
f_value = 4.164e-1
gamma = 6.265e8
lambda0=1215.6700
cLight = 2.99792458e5 #km/s
## define parameters for fitting
HI_params = {'name':'lyz',\
                'field':'H_p0_number_density',\
                'f':[f_value],\
                'Gamma':[gamma],\
                'wavelength':[lambda0],\
                'mass':1.00794,\
                'numLines':1,\
                'maxN':1e22,\
                'minN':1e11,\
                'maxb':300,\
                'minb':1,\
                'init_b':50,\
                'init_N':1e14}
speciesDict={'lya':HI_params}
orderFits = ['lya']
count = 0
if comm_rank == 0:
    widgets = ['Spectra Gen: ', Percentage(), ' ',\
            Bar(marker='0'), ' ', Timer()]
    bar = ProgressBar(widgets=widgets, max_value=len(localRanks)*nPerRank)
    bar.start()
## main loop over spectra
for r in localRanks:
    for i in range(nPerRank):
        sFile ='%s/rayData/%s/RD%04d/spectra/spectrum_%d_%d.h5'\
                    %(dataRepo, simname, d, r, i)
        output = '%s/%d_%d_fitted.h5'%(resultsdir, r, i)
        if os.path.exists(output): 
            count += 1
            if comm_rank == 0: bar.update(count)
            continue
        try:
            f = h5py.File(sFile)
            specFlux = f['flux'][:]
            specLambda = f['wavelength'][:]
            lines, flux = GTF(specLambda, specFlux, \
                        orderFits, speciesDict,
                        output_file=output )
                        
            count += 1
            if comm_rank==0: bar.update(count)
        except:
           print ('bad file')
           continue
counts = comm.reduce(count, root=0, op=MPI.SUM)
if comm_rank==0:
    bar.finish()
    print('Fitted %d spectra'%counts)
