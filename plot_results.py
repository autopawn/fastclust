import numpy as np
import matplotlib.pyplot as plt
import os

RESULT_FOLDER="resultsf"

fnames = os.listdir(RESULT_FOLDER)

optks = set()
dims  = set()
ns    = set()
kps   = set()

def parse_line(lin):
    lin = lin.strip().split()
    optk  = lin[0]
    dim   = int(lin[1])
    n     = int(lin[2])
    kp    = lin[3]
    dists = int(lin[4])
    time  = float(lin[5])#+1e-8
    return optk,dim,n,kp,dists,time


# Read value alternatives
with open(os.path.join(RESULT_FOLDER,fnames[0])) as f:
    for lin in f:
        if len(lin)>0 and lin[0]=="#": continue
        optk,dim,n,kp,_,_ = parse_line(lin)
        optks.add(optk)
        dims.add(dim)
        ns.add(n)
        kps.add(kp)
optks = sorted(list(optks))
dims = sorted(list(dims))
ns = sorted(list(ns))
kps = sorted(list(kps))
print(kps)

ndists = np.zeros((len(optks),len(dims),len(ns),len(kps),len(fnames)),dtype=np.float64)
times  = np.zeros((len(optks),len(dims),len(ns),len(kps),len(fnames)),dtype=np.float64)

ndists[...] = np.nan
times[...] = np.nan

# Read all files
for i,fname in enumerate(fnames):

    # Read files
    with open(os.path.join(RESULT_FOLDER,fname)) as f:
        for lin in f:
            o,d,n,kp,di,t = parse_line(lin)
            if di>=0:
                ndists[optks.index(o),dims.index(d),ns.index(n),kps.index(kp),i] = di
                times[optks.index(o),dims.index(d),ns.index(n),kps.index(kp),i] = t

# Average data along files
ndists = np.mean(ndists,axis=-1)
times = np.mean(times,axis=-1)

# Plot data
for i,dims in enumerate(dims):
    for j,kp in enumerate(kps):
        fig,(ax1,ax2) = plt.subplots(1,2,sharex=True,squeeze=True)
        for k,optk in enumerate(optks):
            xs = ns
            ys = ndists[k,i,:,j]
            ax1.loglog(xs,ys,"o-",label=optk)
        ax1.grid()
        ax1.set_xlabel('$n$')
        ax1.set_ylabel('avg. distances computed')

        for k,optk in enumerate(optks):
            xs = ns
            ys = times[k,i,:,j]
            ax2.loglog(xs,ys,"o-",label=optk)
        ax2.grid()
        ax2.set_xlabel('$n$')
        ax2.set_ylabel('avg. time (s)')

        plt.suptitle("Clustering of points on %d dimensions, $k/n = %s$"%(dims,kp))
        plt.legend()
        plt.show()
