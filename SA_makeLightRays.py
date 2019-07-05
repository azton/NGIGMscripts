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
"""
    3 data locations: DataRepo: where the data Outputs are stored
                      OutsDir: where rays, spectra, etc will be stored
                      scriptsdir: where is this script at?
"""
dataLoc = sys.argv[4]
if dataLoc == 'BW':
    workdir = '/u/sciteam/wells/scratch/prodRun'
    scriptsdir = '/u/sciteam/wells/scratch/prodRun/scripts'
    dataRepo = '/u/sciteam/wells/scratch/prodRun/scripts'
elif dataLoc == 'HD':
    DataRepo = '/home/azton/research/nextGenIGM'
    OutsDir = '/home/azton/research/nextGenIGM/scripts/rayData'
    scriptsdir = OutsDir
elif dataLoc == 'bigbook':
    workdir = '/media/azton/bigbook/projects/nextGenIGM'
    dataRepo = '/media/azton/bigbook/projects/nextGenIGM'
    scriptsdir = '/home/azton/simulations/scripts'
elif dataLoc == 'comet':
    workdir = '/oasis/scratch/comet/azton/temp_project/nextGenIGM'
    dataRepo = '/scratch/azton/%s'%os.getenv('SLURM_JOBID')
    scriptsdir = '%s/scripts'%workdir
else:
    print('Read %s as data location'%dataLoc)
    print('Data Location not recognized!\nUse "BW" for bluewaters, "HD" for harddrive, "bigbook" for external!')
    comm.Abort()
simname = sys.argv[1]
path = '' # path for pipeline file to be put into
if '/' in simname:
    prepath = simname.split('/')[:-1]

    for p in prepath:
        path = path + p + '/'
if not os.path.exists('%s/pipelines/%s'%(scriptsdir,path)):
        os.makedirs('%s/pipelines/%s'%(scriptsdir,path))
d = int(sys.argv[2])
dim = float(sys.argv[3])
dataPath = DataRepo+simname
resultsPath = DataRepo+'/analysis/%s'%simname
if (rank == 0) & (dataLoc == 'comet'):
    if not os.path.exists('%s/rays'%(dataRepo)):
        os.makedirs('%s/rays'%(dataRepo))
    if not os.path.exists("%s/scripts/rayData/%s/RD%04d/"%(workdir, simname, d)):
        os.makedirs("%s/scripts/rayData/%s/RD%04d/"%(workdir, simname, d))
#    if os.path.exists('%s/scripts/rayData/%s/RD%04d/rays.tar'%(workdir, simname, d)):
#        os.system('cd %s; cp %s/scripts/rayData/%s/RD%04d/rays.tar ./; tar -xvf rays.tar'%(dataRepo, workdir, simname, d))
elif (rank == 0) & (dataLoc != "comet"):
    if os.path.exists('%s/%s/RD%04d/rays'%(OutsDir, simname, d)) == False:
        os.makedirs('%s/%s/RD%04d/rays'%(OutsDir, simname, d))
    if os.path.exists('%s/%s/RD%04d/spectra'%(OutsDir, simname, d)) == False:
        os.makedirs('%s/%s/RD%04d/spectra'%(OutsDir, simname, d))

dx = 1./dim
nSamples = int(dim*dim)
try:
    try:
        ds = yt.load('%s/%s/RD%04d/RedshiftOutput%04d'\
                         %(DataRepo, simname, d, d))
        print ('using redshift output at z = %0.3f'%(ds.current_redshift))
    except:
        ds = yt.load('%s/%s/RD%04d/RD%04d'\
                     %(DataRepo, simname, d, d))
        print ('using redshift output at z = %0.3f'%(ds.current_redshift))
except:
    ds = yt.load('%s/%s/DD%08d/DD%08d'\
                %(DataRepo, simname, d, d))
    print ('using data output at z = %0.3f'%(ds.current_redshift))


z = ds.current_redshift
width = ds.domain_width.to("Mpc")[0]
dz = ds.current_redshift - ds.cosmology.z_from_t(\
                            ds.current_time + ds.domain_width[0]/ds.quan(2.998*10**8, 'm/s'))
lambda_shift = 1215.67
lambda_min = lambda_shift*(z-dz+1)
lambda_max = lambda_shift*(z+1)
locals = range(rank, nSamples, size)
LocalID = range(rank, nSamples, size)
count = 0
keepPerRank = int(float(sys.argv[5])//size)
#
# Make light-rays first and save to disk
#
comm.Barrier()
for i in range(keepPerRank):
    startTime = time.time()
    if dataLoc == "comet":
        rFile = "%s/rays/ray_%d_%d.h5"%(dataRepo, rank, count)
    else:
        rFile = '%s/%s/RD%04d/rays/ray_%d_%d.h5'\
                    %(OutsDir,simname, d, rank, count)
    rand=np.random.RandomState(seed=[rank+count, count, count*rank])
    start = rand.uniform(0,1,3)
    a = rand.randint(-100,100, 3)
    shift = a/np.sqrt(np.dot(a,a))
    end = start+shift
    ray = trident.make_simple_ray(ds, start_position=start,\
                end_position=end, data_filename=rFile,\
                lines=['H I 1216'], ftype='gas',\
                fields = ['density','metallicity','temperature','H_p0_number_density'])
    count += 1
    counts = comm.gather(count, root=0)
    if rank==0: 
        if count % 50 == 0:
            print('%d light-ray objects created!\nIntermediate pipelinefile being created'%sum(counts))
            f = open('pipelines/%s_RD%04d_pipeline.info'%(simname, d), 'w')
            f.write('%s\n%d\n%d\n%s\n%d\n%d\n'%(simname,d,dim,sys.argv[4], size, max(counts) ))
            f.close()
            os.system('cd %s; tar -cvf rays.tar rays; mv rays.tar %s/rayData/%s/RD%04d/'\
                          %(dataRepo, scriptsdir, simname, d))

if rank==0: 
    print('%d light-ray objects created!'%sum(counts))
    os.system('touch pipelines/%s_RD%04d_pipeline.info'%(simname, d))
    f = open('pipelines/%s_RD%04d_pipeline.info'%(simname, d), 'w')
    f.write('%s\n%d\n%d\n%s\n%d\n%d\n'%(simname,d,dim,sys.argv[4], size, max(counts) ))
    f.close()
    print('InfoFile: \n%s\n%d\n%d\n%s\n%d\n%d\n'%(simname,d,dim,sys.argv[4], size, max(counts) ))
    if dataLoc=='comet':
        os.system('cd %s; tar -cvf rays.tar rays'%(dataRepo))
        if not os.path.exists('%s/scripts/rayData/%s/RD%04d/'\
                              %(workdir, simname, d)):
            os.makedirs('%s/scripts/rayData/%s/RD%04d/'\
                              %(workdir, simname, d))
        os.system('cp %s/rays.tar %s/scripts/rayData/%s/RD%04d/'%(dataRepo, workdir, simname, d))
