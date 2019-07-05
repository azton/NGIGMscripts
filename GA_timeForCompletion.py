import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
import numpy as np

sims = ['1024_AMR+BG']
workdir = "/oasis/scratch/comet/azton/temp_project/nextGenIGM/"
fig = plt.figure()
ax = fig.add_subplot(111)
for sim in sims:
    print (sim)
    z = []
    time = []
    z.append(99)
    time.append(0)
#    if sim == '1024_AMR+BG':
#        time.append(3600*9)
    print('starting: ',time)
    for run in range(500):
        print(run,)
        tstart = max(time)
        outlog = '%s/%s/estd_%d.out'%(workdir, sim, run)
        fin = 0
        wall = []
        shift = []
        if not os.path.exists(outlog):
            continue
        f = open(outlog, 'r')
        lines = [line for line in f]
        for line in lines:
            if "z =" in line:
                redshift = float(line.split()[-1])
                if True:
                    wallclock = line.split()[-4]
                    wallclock = float(wallclock.strip('='))
                    shift.append(redshift)
                    wall.append(wallclock+tstart)
            if "Exiting" in line:
                fin = 1
        if fin == 0: continue
        for i in range(len(shift)):
            time.append(wall[i])
            z.append(shift[i])
    #z = 1/(1+np.array(z))
    #time = np.array(time)/3600
zs = zip(z,time)
zs.sort()
zs = np.array(zs)
ax.plot(zs[0,:], zs[1,:], linewidth=1.0, label=sim)
ax.legend()
ax.set_yscale('log')
ax.set_ylabel('$T$[hrs]')
ax.set_xlabel('$(1+z)^{-1}$')
plt.savefig('timeToZ.png') 
