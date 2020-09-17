import numpy as np
import matplotlib.pyplot as plt
import os
import sys

RESULT_FOLDER="results"

TITLES = "-t" in sys.argv

STYLES = {
    "op0" : ("k-", {"lw":2, "label":"original"}),
    "op1" : ("-",  {"lw":5, "label":"op1"}),
    "op2" : ("-",  {"lw":4, "label":"op2"}),
    "op3" : ("-",  {"lw":3, "label":"op3"}),
    "op4" : ("--", {"lw":2, "label":"op4"}),
    }

FIGSIZE = (10,7)

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
    time  = float(lin[5])

    return optk,dim,n,kp,dists,time

def div2tfrac(s):
    s = s.split("/")
    assert(len(s)==2)
    return r"\frac{%s}{%s}\, "%tuple(s)

# Read value alternatives
with open(os.path.join(RESULT_FOLDER,fnames[0])) as f:
    for lin in f:
        lin = lin.strip()
        if len(lin)==0: continue
        if len(lin)>0 and lin[0]=="#": continue
        optk,dim,n,kp,_,_ = parse_line(lin)

        optks.add(optk)
        dims.add(dim)
        ns.add(n)
        kps.add(kp)

optks = sorted(list(optks))
optks = optks[1:] + [optks[0]]
dims = sorted(list(dims))
ns = sorted(list(ns))
kps = sorted(list(kps))
print(optks)
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
            lin = lin.strip()
            if len(lin)==0: continue
            o,d,n,kp,di,t = parse_line(lin)
            if di>=0:
                ndists[optks.index(o),dims.index(d),ns.index(n),kps.index(kp),i] = di
                if t>=0:
                    times[optks.index(o),dims.index(d),ns.index(n),kps.index(kp),i] = t


# Average data along files
ndists = np.mean(ndists,axis=-1)
times = np.mean(times,axis=-1)

# Plot data
for j,kp in enumerate(kps):

    # PLOT 1: loglog
    fig,axs = plt.subplots(len(dims),2,sharex=True,sharey='col',squeeze=True,figsize=FIGSIZE)

    for i,ds in enumerate(dims):
        for k,optk in enumerate(optks):
            xs = ns
            ys = ndists[k,i,:,j]
            axs[i][0].loglog(xs,ys,STYLES[optk][0],**STYLES[optk][1])

        axs[i][0].grid()

        for k,optk in enumerate(optks):
            xs = ns
            ys = times[k,i,:,j]
            axs[i][1].loglog(xs,ys,STYLES[optk][0],**STYLES[optk][1])
        axs[i][1].grid()

        axs[i][0].set_ylabel('n. dist. ($D=%d$)'%ds)
        axs[i][1].set_ylabel('time [s] ($D=%d$)'%ds)

    axs[0][0].set_title('distances computed vs. $n$')
    axs[0][1].set_title('CPU time vs. $n$')

    axs[-1][0].set_xlabel('$n$')
    axs[-1][1].set_xlabel('$n$')

    if TITLES:
        plt.suptitle("Computational cost of clustering points with $k = %s n$"%(div2tfrac(kp),))
    plt.legend()

    fname = "img_%d_loglog"%j+("_t" if TITLES else "")+".png"
    plt.savefig(fname,bbox_inches='tight')
    # plt.show()

    # PLOT 2: logx, comparison with op1
    fig,axs = plt.subplots(len(dims),2,sharex=True,sharey='col',squeeze=True,figsize=FIGSIZE)

    for i,ds in enumerate(dims):
        for k,optk in enumerate(optks):
            if optk=="op0": continue
            xs = ns
            ys = ndists[k,i,:,j]/ndists[0,i,:,j]
            axs[i][0].plot(xs,ys,STYLES[optk][0],**STYLES[optk][1])
        axs[i][0].set_xscale('log')
        axs[i][0].set_ylim(bottom=0, top=1.5)
        axs[i][0].grid()

        for k,optk in enumerate(optks):
            if optk=="op0": continue
            xs = ns
            ys = times[k,i,:,j]/times[0,i,:,j]
            axs[i][1].plot(xs,ys,STYLES[optk][0],**STYLES[optk][1])
        axs[i][1].set_xscale('log')
        # axs[i][1].set_ylim(bottom=0, top=10)
        axs[i][1].grid()


        axs[i][0].set_ylabel(r'$\frac{\mathrm{n.\, dist.}}{\mathrm{n.\, dist.\, op1}}$ ($D=%d$)'%ds)
        axs[i][1].set_ylabel(r'$\frac{\mathrm{time}}{\mathrm{time\, op1}}$ ($D=%d$)'%ds)

    axs[0][0].set_title(r'$\frac{\mathrm{distances\,\, computed}}{\mathrm{distances\,\, computed \,\, by \,\, op1}}$ vs. $n$')
    axs[0][1].set_title(r'$\frac{\mathrm{CPU\,\, time}}{\mathrm{CPU\,\, time\,\, op1}}$ vs. $n$')

    axs[-1][0].set_xlabel('$n$')
    axs[-1][1].set_xlabel('$n$')

    if TITLES:
        plt.suptitle("Computational cost of clustering points with $k = %s n$, relative to optimization 1"%(div2tfrac(kp),))
    plt.legend()

    fname = "img_%d_rel"%j+("_t" if TITLES else "")+".png"
    plt.savefig(fname,bbox_inches='tight')
    # plt.show()
