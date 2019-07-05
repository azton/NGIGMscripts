import gc
import yt
import trident
import sys
import numpy as np
import os
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
yt.enable_parallelism()


"""
combined spectra and light-ray generator.  Uses make_simple_ray functionality of trident.
    call using:
    aprun -n nProcs python fullBoxLightRays.py <sim> <dump> <root-grid-dimension> <data location> <N to keep>

    on BW do first:
    > module load bwpy/2.0.0-pre0
    > module load bwpy-mpi
    > source /u/sciteam/wells/my_yt_env/bin/activate
"""


dataLoc = sys.argv[4]
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
    print('Data Location not recognized!\nUse "BW" for bluewaters, "HD" for harddrive, "bigbook" for external!')
    quit()
simname = sys.argv[1]
d = int(sys.argv[2])
dim = float(sys.argv[3])
dataPath = workdir+simname
resultsPath = workdir+'/analysis/%s'%simname
if yt.is_root():
    if os.path.exists('%s/rayData/%s/RD%04d/rays'%(dataRepo, simname, d)) == False:
        os.makedirs('%s/rayData/%s/RD%04d/rays'%(dataRepo, simname, d))
    if os.path.exists('%s/rayData/%s/RD%04d/spectra'%(dataRepo, simname, d)) == False:
        os.makedirs('%s/rayData/%s/RD%04d/spectra'%(dataRepo, simname, d))
dx = 1./dim
nSamples = int(dim*dim)
ds = yt.load('%s/%s/RD%04d/RedshiftOutput%04d'\
                %(workdir, simname, d, d))
z = ds.current_redshift
width = ds.domain_width.to("Mpc")[0]
dz = ds.current_redshift - ds.cosmology.z_from_t(\
                            ds.current_time + ds.domain_width[0]/ds.quan(2.998*10**8, 'm/s'))
lambda_shift = 1215.67
lambda_min = lambda_shift*(z-dz+1)
lambda_max = lambda_shift*(z+1)
locals = range(rank, nSamples, size)
LR = trident.LightRay(ds)
LocalID = range(rank, nSamples, size)
count = 0
keepFrac = float(sys.argv[5])/1024**2
keepPerRank = float(sys.argv[5])//size

#
# Make light-rays first and save to disk
#

for i in yt.parallel_objects(range(nSamples)):
    if i % 5000 == 0: print('%d/%d'%(i,nSamples))
    if np.random.uniform(0,1) > keepFrac: continue # want to cut down and NOT make 1e6 spectra...
    sFile ='%s/rayData/%s/RD%04d/spectra/rank%d_%d.h5'\
                    %(dataRepo, simname, d, rank, count)
    rFile = '%s/rayData/%s/RD%04d/rays/rank%d_%d.h5'\
                    %(dataRepo,simname, d, rank, count)
    if os.path.exists(sFile):
        count += 1
        continue
    x = (i//dim+0.5)*dx
    y = (i%dim+0.5)*dx
    start = [x,y,0]
    end = [x,y, 1]
    ray = LR.make_light_ray(start_position = start,\
                        end_position = end,\
                        use_peculiar_velocity=True,\
                        data_filename = rFile,\
                        fields = ['density','metallicity','temperature','H_p0_number_density'])
    sg = trident.SpectrumGenerator(lambda_min=lambda_min, lambda_max=lambda_max, dlambda=0.01)
    sg.make_spectrum(rFile, lines=['H I 1216'])
    sg.save_spectrum (sFile)
    count += 1
    if count > keepPerRank:
        break
    gc.collect()
if rank==0:
    print('generated %d light-rays!'%comm.reduce(count, root=0, op=MPI.SUM))
