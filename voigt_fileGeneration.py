import yt
import trident
import numpy as np
import sys, os
from mpi4py import MPI
import h5py
from progressbar import *
"""
 After running makeLightRays.py and makeSpectra.py, this will 
 take in spectra files and generate the files for input into autoVP
 This version of autoVP uses 4 columns input:
 <wavelength, velocity, flux, noise>

 USAGE:
    python voigtProfileFitting.py <simname> <output #> 
"""

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
lines = None
if rank == 0:
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
    dataRepo = '/oasis/scratch/comet/azton/temp_project/nextGenIGM'
    workdir = '/scratch/azton/%s'%os.getenv('SLURM_JOB_ID')
    scriptsdir = dataRepo+'/scripts'
else:
    print('Data Location <%s> not recognized!'%dataLoc)
    print('Use "BW" for bluewaters, "HD" for harddrive, "bigbook" for external!')
    comm.Abort()

if dataLoc == 'comet':
    resultsdir = '%s/vpinput'%(workdir)
else:
    resultsdir = '%s/rayData/%s/RD%04d/vpfits/vpinput'%(dataRepo, simname, d)
f_value = 4.164e-1
gamma = 6.265e8
lambda0=1215.6700
cLight = 2.99792458e5 #km/s
count = 0
if rank == 0:
    if not os.path.exists(resultsdir):
        os.makedirs(resultsdir)
    os.system('cp /home/azton/autovp/*fit %s'%resultsdir)
localRanks = range(rank, nRanks, size)
if (dataLoc == 'comet') & (rank == 0):
        os.system('cp %s/rayData/%s/RD%04d/spectra.tar %s/'%(scriptsdir, simname, d, workdir))
        os.system('cd %s; tar -xvf spectra.tar'%workdir)
comm.Barrier()
for rank in localRanks:
    for n in range(nPerRank):
        if dataLoc == 'comet':
              sFile = '%s/spectra/spectra_%d_%d.h5'%(workdir, rank, n)
        else:
              sFile ='%s/rayData/%s/RD%04d/spectra/spectra_%d_%d.h5'\
                %(dataRepo, simname, d, rank, n)
        if not os.path.exists(sFile): continue
        try:
            f = h5py.File(sFile)
        except:
            print('Cannot load %s\n or %s'%(sFile, rfile))
            continue
        try:
            flux = f['flux'][:]
            wavelength = f['wavelength'][:]
            f.close()
        except:
            print('could not load data from file \n%s'%sFile)
            continue
        z0 = wavelength[0]/lambda0 - 1
        velocity = (wavelength/lambda0/(1+z0)-1)*cLight
        noise = np.random.uniform(0.025, 0.040, len(flux))
        out = open('%s/autoVPinput_%d_%d.cln'%(resultsdir, rank, n), 'w')
        for l in range(1, len(flux[:-1])):
            out.write('%0.4f %0.4f %0.6f %0.5f\n'\
                        %(wavelength[l], velocity[l], flux[l], noise[l]))#,\
        out.close()
        count += 1
counts = comm.reduce(count, root=0, op=MPI.SUM)
if rank==0:
    print('%d files written to %s'%(counts, resultsdir))    
    f = open('%s_RD%04d_pipeline.info'%(sys.argv[1], int(sys.argv[2])), 'a')
    f.write('Completed: voigtFileGen with %d spectra\n'%counts)
    f.close()
    
    if dataLoc == 'comet':
        os.system('cd %s; tar -cvf vpinput.tar vpinput'%(workdir))
        os.system('cp %s/vpinput.tar %s/rayData/%s/RD%04d/'%(workdir, scriptsdir, simname, d))
exit()

