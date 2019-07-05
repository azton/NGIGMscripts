
"""
Final script of parametric stats; use after voigt_autovpFitting.py

usage: python -u voigt_plotAutovpFits.py <simname> <dump>

"""


from mpi4py import MPI
import sys
import os
import numpy as np
import yt
import subprocess as SP
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from progressbar import *

def failed (string):
    print ('command failed:\n%s'%string)
    comm.Abort()


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


f = open('%s_RD%04d_pipeline.info'%(sys.argv[1], int(sys.argv[2])), 'r')
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
    workdir = '/scratch/azton/%s'%(os.getenv("SLURM_JOB_ID"))
    scriptsdir = dataRepo+'/scripts'
else:
    print('Data Location <%s> not recognized!'%dataLoc)
    print('Use "BW" for bluewaters, "HD" for harddrive, "bigbook" for external!')
    comm.Abort()

if (dataLoc == "comet") & (rank == 0):
    cmd = SP.run ('cd %s; cp %s/rayData/%s/RD%04d/fitted.tar ./; tar -xvf fitted.tar'\
                      %(workdir, scriptsdir, simname, d), shell = True)
    if cmd.returncode != 0:
        print('Failed Move/untar')
        comm.Abort()
localRanks = range(rank, nRanks, size)
count = 0
comm.Barrier()

## Begin working ##
colDen = []
width = []
if rank == 0:
    widgets = ['Analysis: ', Percentage(), ' ', Bar(marker='%', left='<',\
                right='>'),' ', SimpleProgress(), ' ', Timer()]
    pbar = ProgressBar(widgets=widgets, maxval=len(localRanks)*nPerRank+1)
    pbar.start()
for r in localRanks:
    for i in range(nPerRank):
        if dataLoc == "comet":
            f = "%s/fitted/%d_%d.vpm"%(workdir, r, i)
            if not os.path.exists(f): continue
            N, b = np.loadtxt(f, skiprows=2, usecols=(1,3), unpack=True)
            colDen.append(N)
            width.append(b)
            count += 1
            if rank == 0: pbar.update(count)
if rank == 0: pbar.finish()
colDen = np.array(colDen)
width = np.array(width)
allN = comm.gather(colDen.to_list(), root = 0)
allB = comm.gather(colDen.to_list(), root = 0)
counts = comm.reduce(count, root = 0, op = MPI.SUM)
    
if rank == 0:
    f = open("%s/vpcoldens_widths.txt"%workdir, 'w')
    for line in range(len(allN)):
        f.write("%16.6e          %16.6e\n"%(allN[line], allB[line]))
    f.close()
    cmd = SP.run("cd %s; mv vpcoldens_widths.txt %s/rayData/%s/RD%04d/"\
               %(workdir, scriptsdir, simname, d), shell=True)
    if cmd.returncode != 0:
        print ('Failed to copy final data file')
        comm.Abort()

## basic plot... the data is saved to plot better locally ##

fig, ax = plt.subplots(2,1,figsize = (8,4))
ax[0].set_title("Fitted Column Density dist.")
ax[1].set_title("Fitted Line Width")

ax[0].hist(allN, bins = 50, alpha = 0.5, density = True)
ax[1].hist(allB, bins = 50, alpha = 0.5, density=True)
plt.savefig("%d/analysis/%s/autovpFitted_RD%04d.png"%(dataRepo, simname))
