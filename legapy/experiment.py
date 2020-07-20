import numpy as np
import fastclust

import matplotlib.pyplot as plt

N_SIZES = [20,200,2000,4000]
K_SIZES = [10,100,1000,2000]

def distanceset(a,b):
    if a.shape[0]==1 and b.shape[0]==1:
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


def generate_rand_points(n,p,d):
    return np.random.random((n,p,d))

def generate_dependent_points(l,n,d):
    base  = np.random.random((l,1,d))
    reps  = int( np.ceil(n/l) )
    elems = np.repeat(base,reps,axis=0)[:n]
    for i in range(n):
        # side = (1/l)**(1/d) # regular volume
        # side /= 2 # d/2^d of the volume.
        side = 1e-6
        elems[i,0] += (np.random.random(d)-0.5)*side
        # elems[i,0,np.random.randint(0,d)] = np.random.random()
    return elems

def generate_1hot_points(n,d):
    elems = np.ones((n,1,d))*0.5
    for i in range(n):
        elems[i,0,np.random.randint(0,d)] = np.random.random()
    return elems

def generate_subsets_fom_rand(l,n,p,d):
    points = np.random.random((l,d))
    elems = np.zeros((n,p,d))
    for i in range(n):
        indexes = np.random.randint(0,l,size=p)
        elems[i] = points[ indexes ]
    return elems

cases = [
    lambda n,k: generate_1hot_points(n,d=20),
    lambda n,k: generate_rand_points(n,p=1,d=2),
    lambda n,k: generate_rand_points(n,p=1,d=3),
    lambda n,k: generate_rand_points(n,p=1,d=5),
    lambda n,k: generate_rand_points(n,p=1,d=10),
    lambda n,k: generate_rand_points(n,p=1,d=20),
    # lambda n,k: generate_rand_points(n,p=1,d=10),
    # lambda n,k: generate_rand_points(n,p=1,d=20),
    # lambda n,k: generate_rand_points(n,p=1,d=2),
    # lambda n,k: generate_rand_points(n,p=1,d=3),
    # lambda n,k: generate_dependent_points(l=k,n=n,d=20),
    # lambda n,k: generate_dependent_points(l=k,n=n,d=10),
    # lambda n,k: generate_dependent_points(l=k,n=n,d=5),
    # lambda n,k: generate_dependent_points(l=k,n=n,d=2),
    # lambda n,k: generate_rand_points(n,p=1,d=1),
    # lambda n,k: generate_rand_points(n,p=1,d=5),
]


for i,case in enumerate(cases):
    for n in N_SIZES:
        for k in [x for x in K_SIZES if x<n]:
            elems = case(n,k)
            last_centroids   = None
            last_assigns    = None
            print("========== case %d =========="%i)
            for opt in (6,10,11):#-1,6,10):
                if opt==-1:
                    centroids,assigns,(ndists,qdists) = fastclust.clustering_alt(elems,k,distanceset)
                elif opt==4:
                    centroids,assigns,(ndists,qdists) = fastclust.clustering_opt4(elems,k,distanceset)
                elif opt==5:
                    centroids,assigns,(ndists,qdists) = fastclust.clustering_opt5(elems,k,distanceset)
                elif opt==6:
                    centroids,assigns,(ndists,qdists) = fastclust.clustering_opt6(elems,k,distanceset)
                elif opt==7:
                    centroids,assigns,(ndists,qdists) = fastclust.clustering_opt7(elems,k,distanceset)
                elif opt==10:
                    centroids,assigns,(ndists,qdists) = fastclust.clustering_multicent(elems,k,distanceset)
                elif opt==11:
                    centroids,assigns,(ndists,qdists) = fastclust.clustering_optfinal(elems,k,distanceset)


                efficiency = (n*k)/ndists
                print("n:%8d\tk:%8d\to:%d\tq:%12d\tdist:%10d\te:%.3f"%(n,k,opt,qdists,ndists,efficiency))

                assert(last_centroids is None or np.all(last_centroids==centroids))
                assert(last_assigns is None or np.all(last_assigns==assigns))
                last_centroids = centroids
                last_assigns = assigns
                # _,_,ndists,_ = fastclust.clustering(elems,k,opt_mode=4)
                # print(ndists)

            if elems.shape[2]==2 and elems.shape[1]==1:
                for c in range(k):
                    indexes = [i for i in range(n) if assigns[i]==c]
                    pts_x = elems[indexes][:,0,0]
                    pts_y = elems[indexes][:,0,1]
                    plt.plot(pts_x,pts_y,'o')
                # centroids
                indexes = [i for i in range(n) if i in centroids]
                pts_x = elems[indexes][:,0,0]
                pts_y = elems[indexes][:,0,1]
                plt.plot(pts_x,pts_y,'4k')

                plt.show()
