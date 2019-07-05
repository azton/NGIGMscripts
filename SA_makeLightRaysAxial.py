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


"""
combined spectra and light-ray generator.  Uses make_simple_ray functionality of trident.
    call using:
    aprun -n nProcs python fullBoxLightRays.py <sim> <dump> <root-grid-dimension> <data location> <N to keep>

    on BW do first:
    > module load bwpy/2.0.0-pre0
    > module load bwpy-mpi
    > source /u/sciteam/wells/my_yt_env/bin/activate
"""
if rank == 0:
    scripttime = time.time()
linear = False # whether lightrays are linear along z or randomly oriented.
dataLoc = sys.argv[4]
if dataLoc == 'BW':
    workdir = '/u/sciteam/wells/scratch/prodRun'
    scriptsdir = '/u/sciteam/wells/scratch/prodRun/scripts'
    dataRepo = '/u/sciteam/wells/scratch/prodRun/scripts'
elif dataLoc == 'HD':
    workdir = '/home/azton/simulations'
    dataRepo = '/home/azton/simulations/data'#/media/azton/bigbook/projects/nextGenIGM'
    scriptsdir = '/home/azton/simulations/scripts'
elif dataLoc == 'bigbook':
    workdir = '/media/azton/bigbook/projects/nextGenIGM'
    dataRepo = '/media/azton/bigbook/projects/nextGenIGM'
    scriptsdir = '/home/azton/simulations/scripts'
elif dataLoc == 'comet':
    workdir = '/oasis/scratch/comet/azton/temp_project/nextGenIGM'
    dataRepo = '/scratch/azton/%s'%os.getenv('SLURM_JOBID')
    scriptsdir = '%s/scripts'%workdir
else:
    print('Data Location not recognized!\nUse "BW" for bluewaters, "HD" for harddrive, "bigbook" for external!')
    comm.Abort()
simname = sys.argv[1]
d = int(sys.argv[2])
dim = float(sys.argv[3])
dataPath = workdir+simname
resultsPath = workdir+'/analysis/%s'%simname

if (rank == 0) & (dataLoc == 'comet'):
    if not os.path.exists('%s/lin_rays'%(dataRepo)):
        os.makedirs('%s/lin_rays'%(dataRepo))
    if not os.path.exists("%s/scripts/rayData/%s/RD%04d/"%(workdir, simname, d)):
        os.makedirs("%s/scripts/rayData/%s/RD%04d/"%(workdir, simname, d))
elif (rank == 0) & (dataLoc != "comet"):
    if os.path.exists('%s/rayData/%s/RD%04d/lin_rays'%(dataRepo, simname, d)) == False:
        os.makedirs('%s/rayData/%s/RD%04d/lin_rays'%(dataRepo, simname, d))
    if os.path.exists('%s/rayData/%s/RD%04d/lin_spectra'%(dataRepo, simname, d)) == False:
        os.makedirs('%s/rayData/%s/RD%04d/lin_spectra'%(dataRepo, simname, d))


dx = 1./dim
nSamples = int(dim*dim)
try:
    ds = yt.load('%s/%s/RD%04d/RedshiftOutput%04d'\
                %(workdir, simname, d, d))
    print ('using redshift output at z = %0.3f'%(ds.current_redshift))
except:
    ds = yt.load('%s/%s/DD%04d/data%04d'\
                %(workdir, simname, d, d))
    print ('using data output at z = %0.3f'%(ds.current_redshift))
z = ds.current_redshift
width = ds.domain_width.to("Mpc")[0]
dz = ds.current_redshift - ds.cosmology.z_from_t(\
    ds.current_time + ds.domain_width[0]/ds.quan(2.998*10**8, 'm/s'))
lambda_shift = 1215.67
lambda_min = lambda_shift*(z-dz+1)
lambda_max = lambda_shift*(z+1)
locals = range(rank, nSamples, size)
count = 0
dx = 1./dim
#
# Make light-rays first and save to disk
#
comm.Barrier()
for i in locals:
    startTime = time.time()
    if dataLoc == "comet":
        rFile = "%s/lin_rays/ray_%d.h5"%(dataRepo,i)
    else:
        rFile = '%s/rayData/%s/RD%04d/lin_rays/ray_%d.h5'\
                    %(dataRepo,simname, d, i)
    row = (i // dim)
    col = (i % dim)
    xstart = (row+0.5) * dx
    ystart = (col+0.5) * dx
    start = np.array([xstart, ystart, 0.0])
    end = start + np.array([5*dx,5*dx,1.0])
    ray = trident.make_simple_ray(ds, start_position=start,\
                end_position=end, data_filename=rFile,\
                lines=['H I 1216'], ftype='gas',\
                fields = ['density','metallicity','temperature','H_p0_number_density'])
    count += 1
    counts = comm.gather(count, root=0)
    
    if rank==0: 
        if count % 50 == 0:
            counts=sum(counts)
            print('%d light-ray objects created!'%(counts)\
                      +'\nIntermediate pipelinefile being created')
            f = open('lin_%s_RD%04d_pipeline.info'%(simname, d), 'w')
            f.write('%s\n%d\n%d\n%s\n%d\n%d\n'\
                        %(simname,d,dim,sys.argv[4], size, counts ))
            f.close()
            os.system('cd %s; tar -cvf lin_rays.tar lin_rays;'%(dataRepo)\
                          +' mv lin_rays.tar %s/rayData/%s/RD%04d/'\
                          %(scriptsdir, simname, d))
counts = comm.gather(count, root=0)
if rank==0: 
    counts=sum(counts)
    print('%d light-ray objects created!'%counts)
    os.system('touch lin_%s_RD%04d_pipeline.info'%(simname, d))
    f = open('lin_%s_RD%04d_pipeline.info'%(simname, d), 'w')
    f.write('%s\n%d\n%d\n%s\n%d\n%d\n'%(simname,d,dim,sys.argv[4], size, counts ))
    f.close()
    print('InfoFile: \n%s\n%d\n%d\n%s\n%d\n%d\n'\
              %(simname,d,dim,sys.argv[4], size, counts ))
    if dataLoc=='comet':
        os.system('cd %s; tar -cvf lin_rays.tar lin_rays'%(dataRepo))
        if not os.path.exists('%s/scripts/rayData/%s/RD%04d/'\
                              %(workdir, simname, d)):
            os.makedirs('%s/scripts/rayData/%s/RD%04d/'\
                              %(workdir, simname, d))
        os.system('cp %s/lin_rays.tar %s/scripts/rayData/%s/RD%04d/'\
                      %(dataRepo, workdir, simname, d))
