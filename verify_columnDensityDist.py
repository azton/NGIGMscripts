import yt
import numpy as np
from mpi4py import MPI
import gc
from progressbar import *
yt.enable_parallelism()

"""
Generate a single file that has column densities along 'z' axis.
Use to quantify how many samples are needed to get a full accounting of 
flux PDF and flux power.
"""

outputsDir = '/media/azton/bigbook/projects/nextGenIGM/1024_uni+BG'
dataFileDest = '/home/azton/simulations/scripts/final_dataRepo/1024_uni+BG'
#outputsDir = '/u/sciteam/wells/scratch/prodRun/1024_uni+BG'
#dataFileDest = '/u/sciteam/wells/scratch/prodRun/scripts/final_dataRepo/1024_uni+BG'
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
nPix = 1024.
ds = yt.load("%s/RD0020/RedshiftOutput0020"%outputsDir)
dx = 1/nPix
samples = int(nPix**2)
local = range(rank, samples, size)
n_kept = 100000.
densities = []
cnt = 0
if rank == 0: 
    widgets = ['Analysis: ', Percentage(), ' ', Bar(marker='%', left='<',\
                right='>'),' ', SimpleProgress(), ' ', Timer()]
    pbar = ProgressBar(widgets=widgets, maxval=n_kept+1)
for i in local:
    if np.random.uniform(0,1) < (size*n_kept+10.)/samples:
        redge = (i%nPix)*dx
        ledge = (i%nPix+1)*dx
        top = i//nPix*dx
        bottom = (i//nPix+1)*dx
        right = [redge, top, 0]
        left = [ledge, bottom, 1]
        box = ds.box(right, left)
        cden = (box['H_p0_number_density'].to('1/cm**3')*box['dz'].to('cm')).sum()
        densities.append(cden)
        # print ('%d/%d'%(i, samples))
        cnt += 1
        if rank==0: pbar.update(cnt)
        if cnt % 100 == 0:
            gc.collect()
        if cnt >= n_kept:
            break
sizes = len(densities)
print('[%d] has size: %d'%(rank, sizes))
cnt = comm.reduce(sizes, root=0, op=MPI.SUM)
densities = np.array(densities)
density = None
if rank == 0:
    density = np.empty(cnt, dtype=np.float64) 
comm.Gather(densities, density, root=0)
if rank==0:
    print (density)
    f = open('%s/columnDensities.txt'%(dataFileDest),'w')
    for line in density:
        f.write('%3.16e\n'%line)
    f.close()
    print('%d densities written to file'%len(density))
    pbar.finish()
