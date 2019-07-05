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
    dataRepo = '/oasis/scratch/comet/azton/temp_project/nextGenIGM' 
    scriptsdir = dataRepo+'/scripts'
    workdir = '/scratch/azton/%s'%(os.getenv('SLURM_JOBID'))

else:
    print('Data Location not recognized!\nUse "BW" for bluewaters, "HD" for harddrive, "bigbook" for external!')
    quit()

localRanks = range(rank, nRanks, size)
if dataLoc == 'comet':
    if rank == 0:
        os.system('cd %s;cp  %s/rayData/%s/RD%04d/lin_spectra.tar ./'\
                      %(workdir, scriptsdir, simname, d))
        os.system('cd %s/; tar -xvf lin_spectra.tar'%(workdir))
try:
    ds = yt.load('%s/%s/RD%04d/RedshiftOutput%04d'\
        %(dataRepo, simname, d, d)) 
except:
    ds = yt.load('%s/%s/DD%04d/data%04d'\
        %(dataRepo, simname, d, d)) 
z = ds.current_redshift
Hz = ds.cosmology.hubble_parameter(z).to('km/s/Mpc')
if '80' in simname:
    dx=0.08
else: dx = 0.02
L_com = dx*dim/0.677 #Mpc
print ('Lcom: ',L_com)
fmean = np.loadtxt('%s/rayData/%s/RD%04d/lin_meanFluxPerPixel.txt'\
                        %(scriptsdir,simname,d))
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
for i in range(nPerRank):
        if dataLoc == 'comet':
            file = '%s/lin_spectra/spectra_%d.h5'%(workdir, i)
        else:
            file = '%s/rayData/%s/RD%04d/lin_spectra/spectrum_%d.h5'\
                %(scriptsdir, simname, d, i)
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
    np.savetxt('%s/rayData/%s/RD%04d/lin_fluxPowerWaveNum.txt'\
                    %(scriptsdir, simname, d),\
                    wavenumber)
    np.savetxt('%s/rayData/%s/RD%04d/lin_fluxPower.txt'\
                    %(scriptsdir, simname, d),\
                    powerSums)
    print('Completed power files')
    f = open('lin_%s_RD%04d_pipeline.info'%(sys.argv[1], int(sys.argv[2])), 'a')
    f.write('Completed: SA_fluxPower with %d spectra\n'%cnts)
    f.close()
