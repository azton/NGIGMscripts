"""
Runs all of the scripts to generate voigt fits to the spectra.  
Run after SA_makeSpectra.py, but doesnt depend on PDFs or power.
USAGE:
    python -u RUN_voigtPipeline.py <mpi run command> <simname> <output#>
"""
import os
import sys

procs = int(sys.argv[1])
simname = sys.argv[2]
d = int(sys.argv[3])

# generate files for input to autovp
os.system('ibrun -np %d python -u voigt_fileGeneration.py %s %d'\
            %(procs, simname, d))

# run autovp to fit
os.system('ibrun -np %d python -u voigt_autovpFitting.py %s %d'\
            %(procs, simname, d))

# run the voigt profiler from trident/yt for comparitive purposes
os.system('ibrun -np %d python -u voigt_FitByYT.py %s %d'\
            %(procs, simname, d))
