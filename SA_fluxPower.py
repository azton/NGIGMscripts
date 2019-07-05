import yt
import trident
import sys
import numpy as np
import os
from mpi4py import MPI
import h5py
yt.enable_parallelism()
"""
Open mean flux & spectra to generate the flux power
Use after fullBoxLightRays.py AND generatePDFs.py
    aprun -n nProcs python fluxPower.py <sim> <dump> 

"""

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
f = open('pipelines/%s_RD%04d_pipeline.info'%(sys.argv[1], int(sys.argv[2])), 'r')
lines = [line for line in f]
simname = lines[0].strip()
d = int(lines[1].strip())
dim = int(lines[2].strip())
dataLoc = lines[3].strip()
nRanks = int(lines[4].strip())
nPerRank = int(lines[5].strip())
bins = dim
if dataLoc == 'BW':
    workdir = '/u/sciteam/wells/scratch/prodRun'
    scriptsdir = '/u/sciteam/wells/scratch/prodRun/scripts'
    dataRepo = scriptsdir
elif dataLoc == 'HD':
    DataRepo = '/home/azton/research/nextGenIGM'
    OutsDir = '/home/azton/research/nextGenIGM/scripts/rayData'
    scriptsdir = OutsDir
elif dataLoc == 'bigbook':
    workdir = '/media/azton/bigbook/projects/nextGenIGM'
    scriptsdir = '/home/azton/simulations/scripts'
    dataRepo = workdir
elif dataLoc == 'comet':
    dataRepo = '/oasis/scratch/comet/azton/temp_project/nextGenIGM' 
    scriptsdir = dataRepo+'/scripts'
    workdir = '/scratch/azton/%s'%(os.getenv('SLURM_JOBID'))

else:
    print('Data Location not recognized!\nUse "BW" for bluewaters, "HD" for harddrive, "bigbook" for external!')
    quit()

localRanks = range(rank, nRanks, size)
if dataLoc == 'comet':
    if rank == 0:
        os.system('cd %s;cp  %s/rayData/%s/RD%04d/spectra.tar ./'\
                      %(workdir, scriptsdir, simname, d))
        os.system('cd %s/; tar -xvf spectra.tar'%(workdir))
try:
    try:
        ds = yt.load('%s/%s/RD%04d/RedshiftOutput%04d'\
                         %(dataRepo, simname, d, d))
    except:
        ds = yt.load('%s/%s/RD%04d/RD%04d'\
                         %(dataRepo, simname, d,d))
except:
    ds = yt.load('%s/%s/DD%08d/DD%08d'\
        %(DataRepo, simname, d, d)) 
z = ds.current_redshift
Hz = ds.cosmology.hubble_parameter(z).to('km/s/Mpc')
L_com = ds.domain_width[0].to('Mpccm') #Mpc
fmean = np.loadtxt('%s/%s/RD%04d/meanFluxPerPixel.txt'\
                        %(OutsDir,simname,d))
fmean = np.float(fmean)
if rank == 0: print('using mean flux: %f'%fmean)
L_v = L_com * Hz/(1+z)
dv = L_v/(4*dim) #match the spectra generation lambda binwidth
k_v = 2.0*np.pi/L_v
wavenumber = k_v*np.array(range(1, 4*dim//2))

powerSum = 0
cnt = 0
comm.Barrier()
print('completed init, starting work')
print('Using L_com = %0.3f Mpc, L_v = %0.3f Mpc'%(L_com, L_v))
for r in localRanks:
    for i in range(nPerRank):
        if dataLoc == 'comet':
            file = '%s/spectra/spectra_%d_%d.h5'%(workdir, r,i)
        else:
            file = '%s/%s/RD%04d/spectra/spectrum_%d_%d.h5'\
                %(OutsDir, simname, d, r, i)
        try:
            f = h5py.File(file)
            flux = f['flux'][:]
            f.close()
            flux = flux/fmean-1
            power = np.fft.fft(flux)/L_v*dv
            power = np.abs(power)**2*L_v
            powerSum += power
            cnt += 1
        except:
            continue
powerSums = comm.reduce(powerSum, root=0, op=MPI.SUM)
cnts = comm.reduce(cnt, root=0, op=MPI.SUM)

if rank == 0:
    print('Used %d spectra to generate'%cnts)
    powerSums = np.array(powerSums)/cnts
    powerSums = powerSums[1:4*dim//2]
    np.savetxt('%s/%s/RD%04d/fluxPowerWaveNum.txt'\
                    %(OutsDir, simname, d),\
                    wavenumber)
    np.savetxt('%s/%s/RD%04d/fluxPower.txt'\
                    %(OutsDir, simname, d),\
                    powerSums)
    print('Completed power files')
    f = open('pipelines/%s_RD%04d_pipeline.info'%(sys.argv[1], int(sys.argv[2])), 'a')
    f.write('Completed: SA_fluxPower with %d spectra\n'%cnts)
    f.close()
