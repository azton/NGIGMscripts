import numpy as np
import h5py
import sys
import os
from mpi4py import MPI

"""
This takes in a .h5 file with spectrum, and does statistics to find the 
flux pdf and tau pdf. Uses pipeline info file to get correct parameters.

Usage:
    python generatePDFs.py <pipeline file>
attempt to use parallelism via MPI
"""
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

f = open('%s'%(sys.argv[1]), 'r')
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
    DataRepo = '/home/azton/research/nextGenIGM'
    OutsDir = '/home/azton/research/nextGenIGM/scripts/rayData'
    scriptsdir = OutsDir
elif dataLoc == 'bigbook':
    workdir = '/media/azton/bigbook/projects/nextGenIGM'
    scriptsdir = '/home/azton/simulations/scripts'
    dataRepo = workdir
elif dataLoc == 'comet':
    datarepo = '/oasis/scratch/comet/azton/temp_project/nextGenIGM'
    workdir = '/scratch/%s/%s'%(os.getenv('USER'), os.getenv('SLURM_JOB_ID'))
    scriptsdir = datarepo + '/scripts'
else:
    print('Data Location not recognized!\nUse "BW" for bluewaters, "HD" for harddrive, "bigbook" for external!')
    quit()
if dataLoc == 'comet':
    if rank == 0:
        os.makedirs('%s/rayData/%s/RD%04d'%(workdir, simname, d))

#        os.system('cp %s/scripts/rayData/%s/RD%04d/rays.tar %s/rayData/%s/RD%04d/'\
#                      %(datarepo, simname, d, workdir, simname, d))
#        os.system('cd %s/rayData/%s/RD%04d; tar -xvf rays.tar'\
#                      %(workdir, simname, d))
        os.system('cp %s/scripts/rayData/%s/RD%04d/spectra.tar %s/rayData/%s/RD%04d/'\
                      %(datarepo, simname, d, workdir, simname, d))
        os.system('cd %s/rayData/%s/RD%04d; tar -xvf spectra.tar'\
                      %(workdir, simname, d))



localRanks = range(rank, nRanks, size)
sumFlux = 0.0
cntFlux = 0
binsFlux = np.append(0.0, np.linspace(0.025, 0.975, 20))
binsFlux = np.append(binsFlux, 1.0)
pdfFlux = np.zeros(len(binsFlux)-1)
binsTau = np.linspace(-5,3,101)
pdfTau = np.zeros(100)
count = 0
print ('Rank %d has %d locals to analyze!'%(rank, len(localRanks)*nPerRank))
for r in localRanks:
    for i in range(nPerRank):
        filepath = '%s/%s/RD%04d/spectra/spectra_%d_%d.h5'\
            %(OutsDir, simname, d,r, i)
        if not os.path.exists(filepath): 
            continue
        f = h5py.File(filepath)
        try:
            flux = f['flux'][:]

            sumFlux += np.sum(flux)
            cntFlux += len(flux)
        
            hist, x = np.histogram(flux, binsFlux)
            pdfFlux += hist

            tau = f['tau'][:][np.where(f['tau'][:] > 0)[0]]
            hist, x = np.histogram(np.log10(tau), binsTau)
            pdfTau += hist
            f.close()
            count += 1
        except:
            print ('[%d] Somethings bad in file:\n%s'%(rank,filepath))
            continue
        f.close()
sumAllFlux = comm.reduce(sumFlux, root=0, op=MPI.SUM)
cntAllFlux = comm.reduce(cntFlux, root=0, op=MPI.SUM)

allFluxPDF = comm.reduce(pdfFlux, root=0, op=MPI.SUM)
allTauPDF = comm.reduce(pdfTau, root=0, op=MPI.SUM)
counts = comm.reduce(count, root=0, op=MPI.SUM)
if rank==0:print ('Summed!')
if rank==0: print('Of %d spectra files, used %d in statistics'%(nRanks*nPerRank, cntAllFlux/4096))
if rank==0:
    meanFlux = sumAllFlux/float(cntAllFlux)
    if not os.path.exists('%s/%s/RD%04d'%(OutsDir, simname, d)):
        os.makedirs('%s/%s/RD%04d'%(OutsDir, simname, d))
    f = open('%s/%s/RD%04d/meanFluxPerPixel.txt'%(OutsDir,simname,d), 'w')
    f.write('%23.16e'%meanFlux) 
    f.close()
    normPDFFlux = allFluxPDF/np.sum(allFluxPDF)
    normPDFFlux /= 0.05
    np.savetxt('%s/%s/RD%04d/fluxPDFbins.txt'%(OutsDir,simname, d),binsFlux)
    np.savetxt('%s/%s/RD%04d/fluxPDF.txt'%(OutsDir,simname, d),normPDFFlux)

    normPDFTau = allTauPDF/np.sum(allTauPDF)
    np.savetxt('%s/%s/RD%04d/tauPDFbins.txt'%(OutsDir,simname, d), binsTau)
    np.savetxt('%s/%s/RD%04d/tauPDF.txt'%(OutsDir,simname, d), normPDFTau)
    
    f = open('pipelines/%s_RD%04d_pipeline.info'%(sys.argv[1], int(sys.argv[2])), 'a')
    f.write('\nCompleted: SA_generatePDFs with %d and mean flux %23.16f\n'%(counts, meanFlux))
    f.close()

    
