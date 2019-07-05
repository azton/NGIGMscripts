"""
Make spectra after generating all dem light-rays.

usage: 
python -u SA_makeSpectra.py <simname> <dump>
"""

import gc
import yt
import trident
import sys
import numpy as np
import os
from mpi4py import MPI
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
wall = 60*60*48


## some initialization with pipeline variables ##
f = open('pipelines/%s_RD%04d_pipeline.info'%(sys.argv[1], int(sys.argv[2])), 'r')
lines = [line for line in f]
print (lines)
simname = lines[0].strip()
d = int(lines[1].strip())
dim = int(lines[2].strip())
dataLoc = lines[3].strip()
nRanks = int(lines[4].strip())
nPerRank = int(lines[5].strip())
if dataLoc == 'BW':
    workdir = '/u/sciteam/wells/scratch/prodRun'
    scriptsdir = '/u/sciteam/wells/scratch/prodRun/scripts'
    datarepo = '/u/sciteam/wells/scratch/prodRun/scripts'
elif dataLoc == 'HD':
    DataRepo = '/home/azton/research/nextGenIGM'
    OutsDir = '/home/azton/research/nextGenIGM/scripts/rayData'
    scriptsdir = OutsDir
elif dataLoc == 'bigbook':
    workdir = '/media/azton/bigbook/projects/nextGenIGM'
    datarepo = '/media/azton/bigbook/projects/nextGenIGM'
    scriptsdir = '/home/azton/simulations/scripts'
elif dataLoc == 'comet':
    datarepo = '/oasis/scratch/comet/azton/temp_project/nextGenIGM'
    workdir = '/scratch/azton/%s'%os.getenv('SLURM_JOB_ID')
    scriptsdir = '%s/scripts'%datarepo
else:
    print('Data Location <%s> not recognized!'%dataLoc)
    print('Use "BW" for bluewaters, "HD" for harddrive, "bigbook" for external!')
    quit()

## get to work ##
try:
    try:
        ds = yt.load('%s/%s/RD%04d/RedshiftOutput%04d'\
                %(DataRepo, simname, d, d))
    except: 
        ds = yt.load('%s/%s/RD%04d/RD%04d'\
                         %(DataRepo, simname, d,d))
except:
    ds = yt.load('%s/%s/DD%08d/DD%08d'\
                %(DataRepo, simname, d, d))
z = ds.current_redshift
width = ds.domain_width.to("Mpc")[0]
dz = ds.current_redshift - ds.cosmology.z_from_t(\
                            ds.current_time + ds.domain_width[0]\
                                /ds.quan(2.998*10**8, 'm/s'))
lambda_shift = 1215.67
lambda_min = lambda_shift*(z-dz+1)
lambda_max = lambda_shift*(z+1)
count = 0
localRanks = range(rank, nRanks, size)
start = time.time()
if dataLoc == 'comet':
    if rank == 0:
        if not os.path.exists("%s/spectra"%(workdir)):
            os.makedirs('%s/spectra'%(workdir))

        os.system('cd %s; /bin/cp %s/scripts/rayData/%s/RD%04d/rays.tar ./'\
                      %(workdir, datarepo, simname, d))
        os.system('cd %s; tar -xvf rays.tar'\
                      %(workdir))
#        if os.path.exists("%s/rayData/%s/RD%04d/spectra.tar"\
#                              %(scriptsdir, simname, d)):
#            os.system("cd %s; /bin/cp %s/scripts/rayData/%s/RD%04d/spectra.tar ./"\
#                      %(workdir, datarepo, simname, d))
#            os.system("cd %s; tar -xvf spectra.tar"%workdir)

comm.Barrier()
for r in localRanks:
    for i in range(nPerRank):
            startSpec = time.time()
            if dataLoc == "comet":
                rFile = "%s/rays/ray_%d_%d.h5"%(OutsDir, r, i)
                sFile = "%s/spectra/spectra_%d_%d.h5"%(OutsDir, r, i)
            else:
                rFile = '%s/%s/RD%04d/rays/ray_%d_%d.h5'\
                        %(OutsDir,simname, d, r, i)
                sFile ='%s/%s/RD%04d/spectra/spectrum_%d_%d.h5'\
                        %(OutsDir, simname, d, r, i)
            if not os.path.exists(rFile): continue
            try:
                sg = trident.SpectrumGenerator(\
		    lambda_min=lambda_min, \
		    lambda_max=lambda_max, \
		    dlambda=(lambda_max-lambda_min)/(4*dim))
                sg.make_spectrum(rFile, lines=['H I 1216'])
                sg.save_spectrum (sFile)
                count += 1
                print ('[%d] %0.3f seconds to generate'%(rank, time.time()-startSpec))
                del(sg, startSpec)
            except: continue
            if (rank == 0) & (count % 50 == 0) & (dataLoc == "comet"):
                os.system("cd %s; tar -cvf spectra.tar spectra"\
                              %(workdir))
                os.system("cd %s; cp spectra.tar %s/rayData/%s/RD%04d/"\
                              %(workdir, scriptsdir, simname, d))
counts = comm.gather(count, root=0)
if (rank == 0)& (dataLoc == 'comet') & (count >= nPerRank):
                print("tarring and copying files")
                os.system("cd %s; tar -cvf spectra.tar spectra"\
			      %(workdir))
                os.system("cd %s; cp spectra.tar %s/rayData/%s/RD%04d/"\
			      %(workdir, scriptsdir, simname, d))

