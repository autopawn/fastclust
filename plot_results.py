import numpy as np
import matplotlib.pyplot as plt
import os

RESULT_FOLDER="resultsf"

OPTK_MAP = {'op2_m1':'op3','op2_lo_m1':'op4'}
STYLES = {"op1":('-',5),"op2":('-',4),"op3":('-',3),"op4":('--',2)}
IGNORE_OPTKS=['op2_s', 'op2_s_lo', 'op2_s_lo_m2', 'op2_s_m2']

FIGSIZE = (10,7)

fnames = os.listdir(RESULT_FOLDER)

optks = set()
dims  = set()
ns    = set()
kps   = set()

def parse_line(lin):
    lin = lin.strip().split()
    optk  = lin[0]
    if optk in OPTK_MAP: optk = OPTK_MAP[optk]
    dim   = int(lin[1])
    n     = int(lin[2])
    kp    = lin[3]
    dists = float(lin[4])
    time  = float(lin[5])
    if '/' not in kp:
        kp = round(np.log2(int(kp)/n))
    return optk,dim,n,kp,dists,time

def style(opname):
    if opname in STYLES:
        return STYLES[opname]
    lin = "-"
    if "lo" in opname: lin+= "-"
    if "m2" in opname: lin+= "."
    if "0" in opname: lin+="k"
    lw = 2
    if ("_" not in opname) and ("1" in opname): lw = 5
    if ("_" not in opname) and ("2" in opname): lw = 4
    if ("_" not in opname) and ("3" in opname): lw = 3
    if ("_" in opname): lw = 2
    return (lin,lw)

def div2tfrac(s):
    s = s.split("/")
    assert(len(s)==2)
    return r"\frac{%s}{%s}\, "%tuple(s)

def div2val(s):
    s = s.split("/")
    assert(len(s)==2)
    return float(s[0])/float(s[1])

# Read value alternatives
with open(os.path.join(RESULT_FOLDER,fnames[0])) as f:
    for lin in f:
        lin = lin.strip()
        if len(lin)==0: continue
        if len(lin)>0 and lin[0]=="#": continue
        optk,dim,n,kp,_,_ = parse_line(lin)
        if optk not in IGNORE_OPTKS:
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
            if o not in IGNORE_OPTKS:
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
            lin,lw = style(optk)
            axs[i][0].loglog(xs,ys,lin,lw=lw,label=optk)
        # op0 prediction
        if False:
            xs = np.array(xs,dtype=np.float64)
            ks = np.round(xs * div2val(kp))
            ys = - (ks - 1) * ks * (2 * ks - 3 * xs - 1) / 6;
            axs[i][0].loglog(xs,ys,"k:",lw=lw,label="op0")

        #
        axs[i][0].grid()


        for k,optk in enumerate(optks):
            xs = ns
            ys = times[k,i,:,j]
            lin,lw = style(optk)
            axs[i][1].loglog(xs,ys,lin,lw=lw,label=optk)
        axs[i][1].grid()

        axs[i][0].set_ylabel('n. dist. (%d dims)'%ds)
        axs[i][1].set_ylabel('time [s] (%d dims)'%ds)

    axs[0][0].set_title('distances computed vs. $n$')
    axs[0][1].set_title('CPU time vs. $n$')

    axs[-1][0].set_xlabel('$n$')
    axs[-1][1].set_xlabel('$n$')

    plt.suptitle("Computational cost of clustering points with $k = %s n$"%(div2tfrac(kp),))
    plt.legend()

    plt.savefig("img_%d_loglog.png"%j)
    # plt.show()

    # PLOT 2: logx, comparison with op1
    fig,axs = plt.subplots(len(dims),2,sharex=True,sharey='col',squeeze=True,figsize=FIGSIZE)

    for i,ds in enumerate(dims):
        for k,optk in enumerate(optks):
            if optk=="op0": continue
            xs = ns
            ys = ndists[k,i,:,j]/ndists[0,i,:,j]
            lin,lw = style(optk)
            axs[i][0].plot(xs,ys,lin,lw=lw,label=optk)
        axs[i][0].set_xscale('log')
        axs[i][0].set_ylim(bottom=0, top=1.5)
        axs[i][0].grid()

        for k,optk in enumerate(optks):
            if optk=="op0": continue
            xs = ns
            ys = times[k,i,:,j]/times[0,i,:,j]
            lin,lw = style(optk)
            axs[i][1].plot(xs,ys,lin,lw=lw,label=optk)
        axs[i][1].set_xscale('log')
        axs[i][1].set_ylim(bottom=0, top=8)
        axs[i][1].grid()


        axs[i][0].set_ylabel(r'$\frac{\mathrm{n.\, dist.}}{\mathrm{n.\, dist.\, op1}}$ (%d dims)'%ds)
        axs[i][1].set_ylabel(r'$\frac{\mathrm{time}}{\mathrm{time\, op1}}$ (%d dims)'%ds)

    axs[0][0].set_title(r'$\frac{\mathrm{distances\,\, computed}}{\mathrm{distances\,\, computed \,\, by \,\, op1}}$ vs. $n$')
    axs[0][1].set_title(r'$\frac{\mathrm{CPU\,\, time}}{\mathrm{CPU\,\, time\,\, op1}}$ vs. $n$')

    axs[-1][0].set_xlabel('$n$')
    axs[-1][1].set_xlabel('$n$')

    plt.suptitle("Computational cost of clustering points with $k = %s n$, relative to optimization 1"%(div2tfrac(kp),))
    plt.legend()

    plt.savefig("img_%d_rel.png"%j)
    # plt.show()
