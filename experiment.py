import numpy as np
import fastclust

N_SIZES = [20,200,2000,20000]
K_SIZES = [10,100,1000,10000]
D_SIZES = [2,3,5,10]
P_SIZES = [1,4]

def distanceset(a,b):
    if a.shape[0]==1:
        # Euclidian
        return np.linalg.norm(a-b,axis=1)
    else:
        # Mean geometric error
        sa = a.shape[0]
        sb = b.shape[0]

        tot = 0
        for i in range(sa):
            vec = a[i]
            tot += np.min(np.linalg.norm(b-vec,axis=1))

        for i in range(sb):
            vec = b[i]
            tot += np.min(np.linalg.norm(a-vec,axis=1))

        return tot/(sa+sb)


for d in D_SIZES:
    for p in P_SIZES:
        for n in N_SIZES:
            for k in [x for x in K_SIZES if x<n]:
                elems = np.random.random((n,p,d))
                last_clusters   = None
                last_assigns    = None
                print("==========")
                for opt in (4,):
                    clusters,assigns,(ndists,qdists) = fastclust.clustering_opt4(elems,k,distanceset)
                    print("d:%3d\tp:%3d\to:%d\tn:%8d\tk:%8d\tdist:%12d\tqdist%12d"%(d,p,opt,n,k,ndists,qdists))
                    assert(last_clusters is None or np.all(last_clusters==clusters))
                    assert(last_assigns is None or np.all(last_assigns==assigns))
                    last_clusters = clusters
                    last_assigns = assigns
                    _,_,ndists,_ = fastclust.clustering(elems,k,opt_mode=4)
                    print(ndists)
