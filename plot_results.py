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
    dists = float(lin[4])
    time  = float(lin[5])#+1e-8
    if '/' not in kp:
        kp = round(np.log2(int(kp)/n))
    return optk,dim,n,kp,dists,time

def style(opname):
    lin = "-"
    if "lo" in opname: lin+= "-"
    if "m2" in opname: lin+= "."
    lw = 2
    if ("_" not in opname) and ("1" in opname): lw = 5
    if ("_" not in opname) and ("2" in opname): lw = 4
    if ("_" in opname): lw = 2
    return (lin,lw)


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
for j,kp in enumerate(kps):
    fig,axs = plt.subplots(len(dims),2,sharex=True,sharey='col',squeeze=True)
    for i,ds in enumerate(dims):
        for k,optk in enumerate(optks):
            xs = ns
            ys = ndists[k,i,:,j]
            lin,lw = style(optk)
            axs[i][0].loglog(xs,ys,lin,lw=lw,label=optk)
        axs[i][0].grid()
        axs[i][0].set_xlabel('$n$')
        axs[i][0].set_ylabel('n. dist. (%d dims)'%ds)


        for k,optk in enumerate(optks):
            xs = ns
            ys = times[k,i,:,j]
            lin,lw = style(optk)
            axs[i][1].loglog(xs,ys,lin,lw=lw,label=optk)
        axs[i][1].grid()
        axs[i][1].set_xlabel('$n$')
        axs[i][1].set_ylabel('time [s] (%d dims)'%ds)

    plt.suptitle("Clustering of points with $k/n = %s$"%(kp,))
    plt.legend()
    plt.show()
