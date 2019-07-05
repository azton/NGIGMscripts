"""
    This will call all SA_* scripts, generating Light rays, spectra, flux pdf and power, and the final plot in 
    one python call.
    Use with the same call as makeLightRays, but with required number of ranks:
    python -u RUN_saPipeline.py <simname> <output> <dim> <datalocation> <N_samples> <nRanks>
"""

import os
import sys

simname = sys.argv[1]
d = int(sys.argv[2])
dim = int(sys.argv[3])
dataloc = sys.argv[4]
nKeep = int(sys.argv[5])
nRanks = int(sys.argv[6])
## Call Light rays ##
os.system('ibrun -np %d python -u SA_makeLightRays.py %s %d %d %s %d'\
            %(nRanks, simname, d, dim, dataloc, nKeep))

## Call Spectra gen ##
os.system('ibrun -np %d python -u SA_makeSpectra.py %s %d'\
            %(nRanks, simname, d))

## generate PDFs ##
os.system('ibrun -np %d python -u SA_generatePDFs_nompi.py %s %d'\
            %(nRanks, simname, d))

## Power spectrum ##
os.system('ibrun -np %d python -u SA_fluxPower.py %s %d'\
            %(nRanks, simname, d))

## Plot Results ##
#os.system('python -u SA_plotPowerFlux.py %s %d %s'\
#            %(simname, d, dataloc))
