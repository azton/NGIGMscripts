

"""
Use after LR+spectra.py and voigtFileGeneration.py
Take in a autovp *.cls file, and generate
a fit of N, b parameters.  autofit and minfit input and output
files are fixed in the fortran code; to avoid race conditions, need to
either only use one thread, or create separate directories for each rank
so that H1216N.%%%, when created, is only used by the process that 
made it.
Usage:
python -u voigtProfileFitting.py <simname> <output #>
"""

from mpi4py import MPI
import sys
import os
import numpy as np
import yt
import subprocess as SP



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
    quit()
if dataLoc == 'comet':
    fittedoutputs = '%s/fitted'%(workdir)
else:
    fittedoutputs = '%s/rayData/%s/RD%04d/vpfits/fitted'%(dataRepo, simname, d)
if rank == 0: 
    if not os.path.exists(fittedoutputs): 
        os.makedirs(fittedoutputs)
if dataLoc == "comet":
    inputfiles = '%s/vpinput'%workdir
    if rank == 0:
        mvfiles = SP.run('cd %s; cp %s/rayData/%s/RD%04d/vpinput.tar ./; tar -xvf vpinput.tar'\
                   %(workdir, scriptsdir, simname, d), shell=True)
        if mvfiles.returncode != 0: failed("copying vpinputs")
else:
    inputfiles = '%s/rayData/%s/RD%04d/vpfits/vpinput'%(dataRepo, simname, d)
autofit = '/home/azton/autovp'
localRanks = range(rank, nRanks, size)
count = 0
## create rank-specific directory to copy files and do work ##
if not os.path.exists('%s/rank%d'%(inputfiles, rank)):
    os.makedirs('%s/rank%d'%(inputfiles, rank))
os.system('/bin/cp %s/*fit %s/rank%d'%(autofit, inputfiles, rank))   
rankworkdir =  '%s/rank%d'%(inputfiles, rank)
comm.Barrier()
for r in localRanks:
    for i in range(nPerRank):
        if os.path.exists('%s/%d_%d.vpm'%(fittedoutputs, r, i)):
            count += 1
            continue
        if not os.path.exists('%s/autoVPinput_%d_%d.cln'%(inputfiles, r, i)):
            continue
                
        command = 'cd %s; /bin/cp ../autoVPinput_%d_%d.cln H1216N.cln;'%(rankworkdir, r, i)\
           +' ./autofit H1216N; ./minfit H1216N.pro;'\
           +' /bin/mv H1216N.res %s/%d_%d.res;'%(fittedoutputs,r,i)\
           +' /bin/mv H1216N.vpm %s/%d_%d.vpm'%(fittedoutputs,r,i)
       
        cpCMD = SP.run(command, shell=True)
#        if cpCMD.returncode != 0: failed("making fit")
        
        #os.system(command)
        count += 1
total = comm.reduce(count, root=0, op=MPI.SUM)
if rank == 0:
    print('%d spectra fitted!'%total)
    f = open('%s_RD%04d_pipeline.info'%(sys.argv[1], int(sys.argv[2])), 'a')
    f.write('Completed:  voigt_autovpFitting with %d spectra\n'%total)
    f.close()
    if dataLoc=="comet":
        os.system('cd %s; tar -cvf fitted.tar fitted; cp fitted.tar %s/rayData/%s/RD%04d/'\
                      %(workdir, scriptsdir, simname, d))
