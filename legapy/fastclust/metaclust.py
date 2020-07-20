import numpy as np


def clustering_opt4(X,k,dist,start=0):
    n = len(X)
    # elem id. -> distance to centroid
    prox      = np.zeros(n,dtype=np.float64)
    # elem id. -> cluster id.
    clust     = np.zeros(n,dtype=np.int32)
    # List with elements on each cluster (cluster id. -> [(distance to centroid, elem id.)])
    clusters  = []
    # List with centroids (cluster id. -> elem id.)
    centroids = [start]

    # Current centroid (elem id. -> elem id.)
    cent = lambda i: centroids[clust[i]]

    # Cluster radious (cluster id. -> float64)
    radious = lambda i: clusters[i][-1][0] if len(clusters[i])>0 else 0


    # Count now many times dist is called
    n_dist = 0; q_dist = 0

    # Initialization
    for i in range(n):
        prox[i] = dist(X[start],X[i]); n_dist += 1
        clust[i] = 0
    clusters.append( sorted( [(prox[i],i) for i in range(n) if i!=start] ) )

    for h in range(1,k):
        # Pick new centroid
        Xcand = [clu[-1] for clu in clusters if len(clu)>0]
        _,pick = max(Xcand)
        clusters[clust[pick]].pop() # Remove X[pick] from its cluster

        # Lower bound for cprox
        cprox_lower_bound_right = prox[pick]

        # Closest centroid known on the left
        closest_centroid_left_dist    = prox[pick]
        closest_centroid_left_t       = clust[pick]

        # Set X[pick] cluster to its own
        centroids.append(pick); clust[pick] = h

        # Elements transfered to the new cluster
        new_cluster = []

        # Iterate old centroids
        for t in range(h-1,-1,-1):

            # Old centroid
            ct = centroids[t]

            # Skip cluster if too far away (right lower bound)
            if radious(t) > cprox_lower_bound_right/2:
                # Skip cluster if too far away (left lower bound)
                if t <= closest_centroid_left_t or radious(t) > (prox[ct] - closest_centroid_left_dist)/2:

                    # Distance from old centroid to new
                    cprox = dist(X[ct],X[pick]); n_dist += 1

                    # Update cprox_lower_bound_right
                    cprox_lower_bound_right = max(cprox_lower_bound_right, prox[ct] - cprox)

                    # Find elements so that prox > cprox/2
                    cluster = clusters[t]
                    A = []
                    while len(cluster)>0 and cluster[-1][0] > cprox/2:
                        A.append( cluster.pop() )
                    for _ in range(len(A)):
                        prox1,i   = A.pop()
                        prox2 = dist(X[pick],X[i]); n_dist += 1; q_dist += 1
                        if prox2 < prox1: # move
                            new_cluster.append( (prox2,i) )
                            prox[i]  = prox2
                            clust[i] = h
                        else:             # back to cluster
                            cluster.append( (prox1,i) )
                            assert(prox[i] == prox1)
                            assert(clust[i] == t)

        # Sort new_cluster and add it to the clusters
        new_cluster = sorted(new_cluster)
        clusters.append(new_cluster)

    return centroids,clust,(n_dist,q_dist)

def clustering_opt5(X,k,dist,start=0):
    n = len(X)
    # elem id. -> distance to centroid
    prox      = np.zeros(n,dtype=np.float64)
    # elem id. -> cluster id.
    clust     = np.zeros(n,dtype=np.int32)
    # List with elements on each cluster (cluster id. -> [(distance to centroid, elem id.)])
    clusters  = []
    # List with centroids (cluster id. -> elem id.)
    centroids = [start]
    # Memory of previous clusters and distances to centroid : elem id. -> [(cluster id., distance to centroid)]
    memory    = []

    # Current centroid (elem id. -> elem id.)
    cent = lambda i: centroids[clust[i]]

    # Cluster radious (cluster id. -> float64)
    radious = lambda i: clusters[i][-1][0] if len(clusters[i])>0 else 0


    # Count now many times dist is called
    n_dist = 0; q_dist = 0

    # Initialization
    for i in range(n):
        prox[i] = dist(X[start],X[i]); n_dist += 1
        clust[i] = 0
        memory.append( [(clust[i],prox[i])] )
    clusters.append( sorted( [(prox[i],i) for i in range(n) if i!=start] ) )
    memory[start] = []

    for h in range(1,k):
        # Pick new centroid
        Xcand = [clu[-1] for clu in clusters if len(clu)>0]
        _,pick = max(Xcand)
        clusters[clust[pick]].pop() # Remove X[pick] from its cluster

        # Lower bound for cprox
        cprox_lower_bound_right = prox[pick]

        # Closest centroid known on the left
        closest_centroid_left_t, closest_centroid_left_dist = memory[pick].pop()

        # Set X[pick] cluster to its own
        centroids.append(pick); clust[pick] = h

        # Elements transfered to the new cluster
        new_cluster = []

        # Iterate old centroids
        for t in range(h-1,-1,-1):

            # Old centroid
            ct = centroids[t]

            # Get another closest centroid left
            if t<=closest_centroid_left_t and len(memory[pick])>0:
                closest_centroid_left_t, closest_centroid_left_dist = memory[pick].pop()


            # Skip cluster if too far away (right lower bound)
            if radious(t) > cprox_lower_bound_right/2:
                # Skip cluster if too far away (left lower bound)
                if t <= closest_centroid_left_t or radious(t) > (prox[ct] - closest_centroid_left_dist)/2:

                    # Distance from old centroid to new
                    cprox = dist(X[ct],X[pick]); n_dist += 1

                    # Update cprox_lower_bound_right
                    cprox_lower_bound_right = max(cprox_lower_bound_right, prox[ct] - cprox)

                    # Find elements so that prox > cprox/2
                    cluster = clusters[t]
                    A = []
                    while len(cluster)>0 and cluster[-1][0] > cprox/2:
                        A.append( cluster.pop() )
                    for _ in range(len(A)):
                        prox1,i   = A.pop()
                        prox2 = dist(X[pick],X[i]); n_dist += 1; q_dist += 1
                        if prox2 < prox1: # move
                            new_cluster.append( (prox2,i) )
                            prox[i]  = prox2
                            clust[i] = h
                            memory[i].append( (clust[i],prox[i]) )
                        else:             # back to cluster
                            cluster.append( (prox1,i) )
                            assert(prox[i] == prox1)
                            assert(clust[i] == t)

        # Sort new_cluster and add it to the clusters
        new_cluster = sorted(new_cluster)
        clusters.append(new_cluster)

    return centroids,clust,(n_dist,q_dist)

def clustering_opt6(X,k,dist,start=0):
    n = len(X)
    # elem id. -> distance to centroid
    prox      = np.zeros(n,dtype=np.float64)
    # elem id. -> cluster id.
    clust     = np.zeros(n,dtype=np.int32)
    # List with elements on each cluster (cluster id. -> [(distance to centroid, elem id.)])
    clusters  = []
    # List with centroids (cluster id. -> elem id.)
    centroids = [start]
    # Memory of previous clusters and distances to centroid : elem id. -> [(cluster id., distance to centroid)]
    memory    = []


    # Current centroid (elem id. -> elem id.)
    cent = lambda i: centroids[clust[i]]

    # Cluster radious (cluster id. -> float64)
    radious = lambda i: clusters[i][-1][0] if len(clusters[i])>0 else 0


    # Count now many times dist is called
    n_dist = [0]; q_dist = 0

    mim = {}
    def distmem(i,j):
        a = min(i,j)
        b = max(i,j)
        if (a,b) in mim:
            return mim[(a,b)]
        else:
            d = dist(X[a],X[b])
            mim[(a,b)] = d
            n_dist[0] += 1
            return d

    def indistmem(i,j):
        a = min(i,j)
        b = max(i,j)
        return (a,b) in mim

    # Initialization
    for i in range(n):
        prox[i] = distmem(start,i)
        clust[i] = 0
        memory.append( [(clust[i],prox[i])] )
    clusters.append( sorted( [(prox[i],i) for i in range(n) if i!=start] ) )
    memory[start] = []

    for h in range(1,k):
        # Pick new centroid
        Xcand = [clu[-1] for clu in clusters if len(clu)>0]
        _,pick = max(Xcand)
        clusters[clust[pick]].pop() # Remove X[pick] from its cluster

        # Lower bound for cprox
        cprox_lower_bound_right = prox[pick]

        # Closest centroid known on the left
        closest_centroid_left_t, closest_centroid_left_dist = memory[pick].pop()

        # Set X[pick] cluster to its own
        centroids.append(pick); clust[pick] = h

        # Elements transfered to the new cluster
        new_cluster = []

        # Iterate old centroids
        for t in range(h-1,-1,-1):

            # Old centroid
            ct = centroids[t]

            # Get another closest centroid left
            if t<=closest_centroid_left_t and len(memory[pick])>0:
                closest_centroid_left_t, closest_centroid_left_dist = memory[pick].pop()

            #
            cprox_lo = cprox_lower_bound_right
            if t>closest_centroid_left_t:
                cprox_lo = prox[ct] - closest_centroid_left_dist
            if indistmem(ct,pick):
                cprox_lo = distmem(ct,pick)

            # Skip cluster if too far away (right lower bound)
            if radious(t) > cprox_lo/2:
                # Distance from old centroid to new
                cprox = distmem(ct,pick)

                # Find elements so that prox > cprox/2
                cluster = clusters[t]
                A = []
                while len(cluster)>0 and cluster[-1][0] > cprox/2:
                    A.append( cluster.pop() )
                for _ in range(len(A)):
                    prox1,i   = A.pop()
                    prox2 = distmem(pick,i); q_dist += 1
                    if prox2 < prox1: # move
                        new_cluster.append( (prox2,i) )
                        prox[i]  = prox2
                        clust[i] = h
                        memory[i].append( (clust[i],prox[i]) )
                    else:             # back to cluster
                        cluster.append( (prox1,i) )
                        assert(prox[i] == prox1)
                        assert(clust[i] == t)

            # Update cprox_lower_bound_right
            cprox_lower_bound_right = max(cprox_lower_bound_right, prox[ct] - cprox_lo)

        # Sort new_cluster and add it to the clusters
        new_cluster = sorted(new_cluster)
        clusters.append(new_cluster)

    return centroids,clust,(n_dist[0],q_dist)

def clustering_opt7(X,k,dist,start=0):
    n = len(X)
    # elem id. -> distance to centroid
    prox      = np.zeros(n,dtype=np.float64)
    # elem id. -> cluster id.
    clust     = np.zeros(n,dtype=np.int32)
    # List with elements on each cluster (cluster id. -> [(distance to centroid, elem id.)])
    clusters  = []
    # List with centroids (cluster id. -> elem id.)
    centroids = [start]
    # Memory of previous clusters and distances to centroid : elem id. -> [(cluster id., distance to centroid)]
    memory    = []


    # Current centroid (elem id. -> elem id.)
    cent = lambda i: centroids[clust[i]]

    # Cluster radious (cluster id. -> float64)
    radious = lambda i: clusters[i][-1][0] if len(clusters[i])>0 else 0
    radious2 = lambda i: clusters[i][-2][0] if len(clusters[i])>1 else 0


    # Count now many times dist is called
    n_dist = [0]; q_dist = 0

    mim = {}
    def distmem(i,j):
        a = min(i,j)
        b = max(i,j)
        if (a,b) in mim:
            return mim[(a,b)]
        else:
            d = dist(X[a],X[b])
            mim[(a,b)] = d
            n_dist[0] += 1
            return d


    # Initialization
    for i in range(n):
        prox[i] = distmem(start,i)
        clust[i] = 0
        memory.append( [(clust[i],prox[i])] )
    clusters.append( sorted( [(prox[i],i) for i in range(n) if i!=start] ) )
    memory[start] = []

    for h in range(1,k):
        # Pick new centroid
        Xcand = [clu[-1] for clu in clusters if len(clu)>0]
        _,pick = max(Xcand)
        clusters[clust[pick]].pop() # Remove X[pick] from its cluster

        # Lower bound for cprox
        cprox_lower_bound_right = prox[pick]

        # Closest centroid known on the left
        closest_centroid_left_t, closest_centroid_left_dist = memory[pick].pop()

        # Set X[pick] cluster to its own
        centroids.append(pick); clust[pick] = h

        # Elements transfered to the new cluster
        new_cluster = []

        # Iterate old centroids
        for t in range(h-1,-1,-1):

            # Old centroid
            ct = centroids[t]

            # Get another closest centroid left
            if t<=closest_centroid_left_t and len(memory[pick])>0:
                closest_centroid_left_t, closest_centroid_left_dist = memory[pick].pop()

            cprox_lower_bound = cprox_lower_bound_right
            if t > closest_centroid_left_t:
                cprox_lower_bound = prox[ct] - closest_centroid_left_dist

            # Skip cluster if too far away (right lower bound)
            if radious(t) > cprox_lower_bound/2:

                # Distance from old centroid to new
                #if radious2(t) > cprox_lower_bound/2:
                cprox_lower_bound = distmem(ct,pick)

                # Find elements so that prox > cprox_lower_bound/2
                cluster = clusters[t]
                A = []
                while len(cluster)>0 and cluster[-1][0] > cprox_lower_bound/2:
                    A.append( cluster.pop() )
                for _ in range(len(A)):
                    prox1,i   = A.pop()
                    prox2 = distmem(pick,i); q_dist += 1
                    if prox2 < prox1: # move
                        new_cluster.append( (prox2,i) )
                        prox[i]  = prox2
                        clust[i] = h
                        memory[i].append( (clust[i],prox[i]) )
                    else:             # back to cluster
                        cluster.append( (prox1,i) )
                        assert(prox[i] == prox1)
                        assert(clust[i] == t)

            # Update cprox_lower_bound_right
            cprox_lower_bound_right = max(cprox_lower_bound_right, prox[ct] - cprox_lower_bound)

        # Sort new_cluster and add it to the clusters
        new_cluster = sorted(new_cluster)
        clusters.append(new_cluster)

    return centroids,clust,(n_dist[0],q_dist)

def clustering_alt(X,k,dist,start=0):
    n = len(X)
    # elem id. -> distance to centroid
    prox      = np.zeros(n,dtype=np.float64)
    # elem id. -> cluster id.
    clust     = np.zeros(n,dtype=np.int32)
    # List with elements on each cluster (cluster id. -> [(distance to centroid, elem id.)])
    clusters  = []
    # List with centroids (cluster id. -> elem id.)
    centroids = [start]
    # Memory of previous clusters and distances to centroid : elem id. -> [(cluster id., distance to centroid)]
    memory    = []

    moves = 0

    # Current centroid (elem id. -> elem id.)
    cent = lambda i: centroids[clust[i]]

    # Count now many times dist is called
    n_dist = [0];

    mim = {}
    def distmem(i,j):
        a = min(i,j)
        b = max(i,j)
        if (a,b) in mim:
            return mim[(a,b)]
        else:
            d = dist(X[a],X[b])
            mim[(a,b)] = d
            n_dist[0] += 1
            return d


    # Initialization
    for i in range(n):
        prox[i] = distmem(start,i)
        clust[i] = 0
        memory.append( [(clust[i],prox[i])] )
    clusters.append( sorted( [(prox[i],i) for i in range(n) if i!=start] ) )
    memory[start] = []

    for h in range(1,k):
        # Pick new centroid
        Xcand = [clu[-1] for clu in clusters if len(clu)>0]
        _,pick = max(Xcand)
        clusters[clust[pick]].pop() # Remove X[pick] from its cluster

        # Set X[pick] cluster to its own
        centroids.append(pick); clust[pick] = h

        # Elements transfered to the new cluster
        new_cluster = []

        cprox = []
        # Iterate old centroids
        for t in range(h):
            # Old centroid
            ct = centroids[t]
            # Distance to old centroid
            cprox.append( distmem(ct,pick) )

        # Iterate old centroids
        for t in range(h):
            # Old centroid
            ct = centroids[t]

            cprox_lower_bound = cprox[t]

            # Find elements so that prox > cprox_lower_bound/2
            cluster = clusters[t]
            A = []
            while len(cluster)>0 and cluster[-1][0] > cprox_lower_bound/2:
                A.append( cluster.pop() )
            for _ in range(len(A)):
                prox1,i   = A.pop()
                # Check if i passes all the tests to move to the new cluster
                passes = 1
                for clus,proxp in memory[i][-2::-1]:
                    if proxp < cprox[clus]/2:
                        passes = 0
                        break

                if passes:
                    prox2 = distmem(pick,i);
                    if prox2 < prox1: # move
                        new_cluster.append( (prox2,i) )
                        prox[i]  = prox2
                        clust[i] = h
                        memory[i].append( (clust[i],prox[i]) )
                        moves += 1
                    else:
                        passes = 0

                if not passes:             # back to cluster
                    cluster.append( (prox1,i) )
                    assert(prox[i] == prox1)
                    assert(clust[i] == t)

        # Sort new_cluster and add it to the clusters
        new_cluster = sorted(new_cluster)
        clusters.append(new_cluster)

    return centroids,clust,(n_dist[0],moves)

def clustering_multicent(X,k,dist,start=0):
    n = len(X)
    # elem id. -> distance to centroid
    prox      = np.zeros(n,dtype=np.float64)
    # elem id. -> cluster id.
    clust     = np.zeros(n,dtype=np.int32)
    # List with elements on each cluster (cluster id. -> [(distance to centroid, elem id.)])
    clusters  = []
    # List with centroids (cluster id. -> elem id.)
    centroids = [start]
    # Memory of clusters and distances to centroid that were evaluated : elem id. -> [(cluster id., distance to centroid)]
    memoryX   = []


    # Current centroid (elem id. -> elem id.)
    cent = lambda i: centroids[clust[i]]

    # Cluster radious (cluster id. -> float64)
    radious = lambda i: clusters[i][-1][0] if len(clusters[i])>0 else 0


    # Count now many times dist is called
    n_dist = [0]; q_dist = 0

    mim = {}
    def distmem(i,j):
        a = min(i,j)
        b = max(i,j)
        if (a,b) in mim:
            return mim[(a,b)]
        else:
            d = dist(X[a],X[b])
            mim[(a,b)] = d
            n_dist[0] += 1
            return d

    def indistmem(i,j):
        a = min(i,j)
        b = max(i,j)
        return (a,b) in mim

    # Initialization
    for i in range(n):
        prox[i] = distmem(start,i)
        clust[i] = 0
        memoryX.append( [(clust[i],prox[i])] )
    clusters.append( sorted( [(prox[i],i) for i in range(n) if i!=start] ) )

    for h in range(1,k):
        # Pick new centroid
        Xcand = [clu[-1] for clu in clusters if len(clu)>0]
        _,pick = max(Xcand)
        clusters[clust[pick]].pop() # Remove X[pick] from its cluster

        # Set X[pick] cluster to its own cluster
        centroids.append(pick); clust[pick] = h

        # Elements transfered to the new cluster
        new_cluster = []

        # Iterate old centroids
        cproxes = [distmem(pick,centroids[t]) for t in range(h)]

        for t in range(h):
            # Old centroid
            ct = centroids[t]

            # Get another closest centroid left
            cprox = cproxes[t]

            # Skip cluster if too far away (right lower bound)
            if radious(t) > cprox/2:

                # Find elements so that prox > cprox/2
                cluster = clusters[t]
                A = []
                while len(cluster)>0 and cluster[-1][0] > cprox/2:
                    A.append( cluster.pop() )
                for _ in range(len(A)):
                    prox1,i   = A.pop()

                    passes = True
                    for cc,dd in memoryX[i][-2::-1]: # is [-2::-1] right?
                        bound = abs(cproxes[cc]-dd)
                        if prox1 <= bound:
                            passes = False
                            break

                    if passes:
                        # Now compute distance to be sure
                        prox2 = distmem(pick,i); q_dist += 1

                        # If we have come to compute the distance to the centroid it may be worth saving:
                        memoryX[i].append( (h,prox2) )

                        # Update passes
                        passes = prox2 < prox1

                    if passes:
                        new_cluster.append( (prox2,i) )
                        prox[i]  = prox2
                        clust[i] = h
                    else:
                        # back to cluster
                        cluster.append( (prox1,i) )
                        assert(prox[i] == prox1)
                        assert(clust[i] == t)


        # Sort new_cluster and add it to the clusters
        new_cluster = sorted(new_cluster)
        clusters.append(new_cluster)

    return centroids,clust,(n_dist[0],q_dist)


def clustering_optfinal(X,k,dist,start=0):
    n = len(X)
    # elem id. -> distance to centroid
    prox      = np.zeros(n,dtype=np.float64)
    # elem id. -> cluster id.
    clust     = np.zeros(n,dtype=np.int32)
    # List with elements on each cluster (cluster id. -> [(distance to centroid, elem id.)])
    clusters  = []
    # List with centroids (cluster id. -> elem id.)
    centroids = [start]
    # Memory of clusters and distances to centroid that were evaluated : elem id. -> [(cluster id., distance to centroid)]
    memory   = []


    # Current centroid (elem id. -> elem id.)
    cent = lambda i: centroids[clust[i]]

    # Cluster radious (cluster id. -> float64)
    radious = lambda i: clusters[i][-1][0] if len(clusters[i])>0 else 0


    # Count now many times dist is called
    n_dist = [0]; q_dist = 0

    mim = {}
    def distmem(i,j):
        a = min(i,j)
        b = max(i,j)
        if (a,b) in mim:
            return mim[(a,b)]
        else:
            d = dist(X[a],X[b])
            mim[(a,b)] = d
            n_dist[0] += 1
            return d

    def indistmem(i,j):
        a = min(i,j)
        b = max(i,j)
        return (a,b) in mim

    # Initialization
    for i in range(n):
        prox[i] = distmem(start,i)
        clust[i] = 0
        memory.append( [(clust[i],prox[i])] )

    clusters.append( sorted( [(prox[i],i) for i in range(n) if i!=start] ) )
    memory[start] = []

    for h in range(1,k):
        # Pick new centroid
        Xcand = [clu[-1] for clu in clusters if len(clu)>0]
        _,pick = max(Xcand)
        clusters[clust[pick]].pop() # Remove X[pick] from its cluster

        # Lower bound for cprox
        cprox_lower_bound_right = prox[pick]

        # Closest centroid known on the left
        closest_centroid_left_t, closest_centroid_left_dist = memory[pick].pop()

        # Set X[pick] cluster to its own
        centroids.append(pick); clust[pick] = h

        # Elements transfered to the new cluster
        new_cluster = []

        # Find lower and upper bounds cproxs
        cprox_lo = np.ones(h,dtype=np.float64)
        cprox_hi = np.ones(h,dtype=np.float64)

        # Iterate old centroids
        for t in range(h-1,-1,-1):
            # Old centroid
            ct = centroids[t]

            # Simple hi bound
            cprox_hi[t] = np.inf #prox[ct] # cprox_lower_bound_right

            # Get another closest centroid left
            if t<=closest_centroid_left_t and len(memory[pick])>0:
                closest_centroid_left_t, closest_centroid_left_dist = memory[pick].pop()

            #
            cprox_lo[t] = cprox_lower_bound_right
            if t>closest_centroid_left_t:
                cprox_lo[t] = prox[ct] - closest_centroid_left_dist
            if indistmem(ct,pick):
                cprox_lo[t] = distmem(ct,pick)
                cprox_hi[t] = cprox_lo[t]

            # Skip cluster if too far away (right lower bound)
            if radious(t) > cprox_lo[t]/2:
                # Distance from old centroid to new
                cprox_lo[t] = distmem(ct,pick)
                cprox_hi[t] = cprox_lo[t]

            # Update cprox_lower_bound_right
            cprox_lower_bound_right = max(cprox_lower_bound_right, prox[ct] - cprox_lo[t])

        # Iterate old centroids
        for t in range(h):

            # Old centroid
            ct = centroids[t]

            # Skip cluster if too far away (right lower bound)
            if radious(t) > cprox_lo[t]/2:
                # Find elements so that prox > cprox/2
                cluster = clusters[t]
                A = []
                while len(cluster)>0 and cluster[-1][0] > cprox_lo[t]/2:
                    A.append( cluster.pop() )
                for _ in range(len(A)):
                    prox1,i   = A.pop()

                    # Test against memory
                    passes = True
                    for cc,dd in memory[i][-2::-1]:
                        bound1 = cprox_lo[cc]-dd
                        bound2 = dd-cprox_hi[cc]

                        if prox1 <= bound1 or prox1 <= bound2:
                            passes = False
                            break

                    if passes:
                        # Now compute distance to be sure
                        prox2 = distmem(pick,i); q_dist += 1

                        # If we have come to compute the distance to the centroid it may be worth saving:
                        memory[i].append( (h,prox2) )

                        # Update passes
                        passes = prox2 < prox1

                    if passes:
                        new_cluster.append( (prox2,i) )
                        prox[i]  = prox2
                        clust[i] = h
                    else:
                        # back to cluster
                        cluster.append( (prox1,i) )
                        assert(prox[i] == prox1)
                        assert(clust[i] == t)



        # Sort new_cluster and add it to the clusters
        new_cluster = sorted(new_cluster)
        clusters.append(new_cluster)

    return centroids,clust,(n_dist[0],q_dist)
