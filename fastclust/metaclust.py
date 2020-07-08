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
        closest_centroid_left_cluster = clust[pick]

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
                if t <= closest_centroid_left_cluster or radious(t) > (prox[ct] - closest_centroid_left_dist)/2:

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
