import numpy as np

from numba import jit

# @jit(nopython=True)
def normaxis1(v):
    res = np.zeros(v.shape[0])
    for i in range(v.shape[0]):
        res[i] = np.linalg.norm(v[i,:])
    return res


# @jit(nopython=True)
def dist(a,b):
    assert(len(a.shape)==2)
    if a.shape[0]==1:
        # Euclidian
        return normaxis1(a-b)
    else:
        # Mean geometric error
        sa = a.shape[0]
        sb = b.shape[0]

        tot = 0
        for i in range(sa):
            vec = a[i]
            tot += np.min(normaxis1(b-vec))

        for i in range(sb):
            vec = b[i]
            tot += np.min(normaxis1(a-vec))

        return tot/(sa+sb)



# @jit(nopython=True)
def clustering(elems,n,start_centroid=0,opt_mode=4):
    m = len(elems)

    n_calls = 0

    if n>=m:
        # No need for clustering, every element is a cluster
        centroids = list(range(m)) #np.arange(m,dtype=np.int32)
        cluster   = np.arange(m,dtype=np.int32)
        centroid_d = np.zeros(m,dtype=np.float64)
        return centroids,cluster,0,centroid_d

    # Centroid indexes: cluster -> centroid
    centroids  = [42][:0] #
    # If a given index is centroid: elemindex -> boolean
    is_centroid = np.zeros(m,dtype=np.int32)
    # Current cluster for each elem: elemindex -> cluster
    cluster    = np.zeros(m,dtype=np.int32)
    # Distance to the centroid; or to last centroid if this elem became centroid: elemindex -> distance
    centroid_d = np.zeros(m,dtype=np.float64)


    # If start_centroid not given, pick elemnt furthest to element 0
    if start_centroid is None:
        for i in range(m):
            centroid_d[i] = dist(elems[0],elems[i])
        start_centroid = np.argmax(centroid_d)

    # Use start_centroid
    for i in range(m):
        centroid_d[i] = dist(elems[start_centroid],elems[i])
        n_calls += 1
    centroids.append(start_centroid)
    is_centroid[start_centroid]   = 1
    centroid_d[start_centroid]    = np.inf


    if opt_mode==0:

        # Add a new centroid until n centroids are selected
        for k in range(1,n):
            # ^ k is the number of centroids now and also the index of the next centroid

            # Find new centroid (ki is centroid of cluster k)
            old_ki_centroid_dist,ki = max([(centroid_d[i],i) for i in range(m) if not is_centroid[i]])

            # Mark ki as centroid
            is_centroid[ki] = True
            centroids.append(ki)
            cluster[ki] = k

            # Update centroids
            for e in range(m):
                # No need to update if centroid
                if is_centroid[e]: continue

                distance = dist(elems[e],elems[ki])
                n_calls += 1

                if distance < centroid_d[e]:
                    centroid_d[e] = distance
                    cluster[e]    = k

    else:

        # Cluster radio (distance to furthest index): cluster -> distance
        cluster_radio = np.zeros(1,dtype=np.float64)
        cluster_size = np.zeros(1,dtype=np.float64)
        # Set initial radio
        cluster_radio[0] = np.max(centroid_d)
        cluster_size[0] = m-1

        # Centroid of the centroid before it becomes a centroid
        cluster_old_cluster = [-1]

        # Add a new centroid until n centroids are selected
        saved1 = 0
        saved12 = 0

        for k in range(1,n):
            # ^ k is the number of centroids now and also the index of the next centroid

            # Find new centroid (ki is centroid of cluster k) # NOTE: may replace with priority queue
            old_ki_centroid_dist,ki = max([(centroid_d[i],i) for i in range(m) if not is_centroid[i]])

            # Old centroid of ki
            old_ki_cluster  = cluster[ki]
            old_ki_centroid = centroids[cluster[ki]]

            # Mark ki as centroid
            is_centroid[ki] = True
            centroids.append(ki)
            cluster[ki] = k

            #
            cluster_old_cluster.append(old_ki_cluster)

            # List of affected old centroids
            is_affected = np.zeros(k,dtype=np.int32)
            # Distance from affected cluster's centroids to ki: cluster -> dist
            ki_dist     = np.zeros(k,dtype=np.float64)

            debug_n_dists_computed = 0


            # The lower bound for the distance of ki to the centroids that remains to be checked
            lower_bound_ki_dist = 0

            if opt_mode==4:
                # Find clusters whose elements may be transfered to a new centroid (affected clusters)

                for i in range(k-1,-1,-1): #reversed(range(k)):        # NOTE: may be replaced with reverse index + just iterate affected clusters
                    # Index of cluster i centroid
                    ci = centroids[i]

                    # Centroid is not affected under this condition (if i is older than old_ki_cluster)
                    if i > old_ki_cluster and (centroid_d[ci] - old_ki_centroid_dist)/2 >= cluster_radio[i]:
                        saved1 += 1
                        continue

                    # Can centroid i elements be stolen by new centroid ki?
                    if lower_bound_ki_dist/2 < cluster_radio[i]:

                        # Compute distance to new centroid ki
                        if i!=old_ki_cluster:
                            ki_dist[i] = dist(elems[ci],elems[ki])
                            n_calls += 1
                        else:
                            ki_dist[i] = old_ki_centroid_dist
                        debug_n_dists_computed += 1

                        # Update lower bound if it increases
                        lower_bound_ki_dist = max(lower_bound_ki_dist, centroid_d[ci] - ki_dist[i])

                        # Check now if cluster i is affected by new cluster k
                        if ki_dist[i]/2 < cluster_radio[i]:
                            is_affected[i] = True
                    else:
                        saved12 += 1

            elif opt_mode==5:

                for i in range(k-1,-1,-1):
                    # Index of cluster i centroid
                    ci = centroids[i]


                    # Centroid is not affected under this condition (if i is older than old_ki_cluster)
                    assert(i>=old_ki_cluster)
                    if i == old_ki_cluster or (centroid_d[ci] - old_ki_centroid_dist)/2 < cluster_radio[i]:
                        # Can centroid i elements be stolen by new centroid ki?
                        if lower_bound_ki_dist/2 < cluster_radio[i]:

                            # Compute distance to new centroid ki
                            if i!=old_ki_cluster:
                                ki_dist[i] = dist(elems[ci],elems[ki])
                                n_calls += 1
                            else:
                                ki_dist[i] = old_ki_centroid_dist #or np.isposinf(old_ki_centroid_dist)
                            debug_n_dists_computed += 1

                            # Update lower bound if it increases
                            lower_bound_ki_dist = max(lower_bound_ki_dist, centroid_d[ci] - ki_dist[i])

                            # Check now if cluster i is affected by new cluster k
                            if ki_dist[i]/2 < cluster_radio[i]:
                                is_affected[i] = True
                        else:
                            saved12 += 1
                    else:
                        saved1 += 1


                    # Move close known candidate
                    if i == old_ki_cluster:
                        # Move to previous close known candidate
                        old_ki_cluster       = cluster_old_cluster[old_ki_cluster]
                        assert(old_ki_cluster<i)

                        if old_ki_cluster>=0 and lower_bound_ki_dist/2 < cluster_radio[old_ki_cluster]:
                            old_ki_centroid_dist = dist(elems[centroids[old_ki_cluster]],elems[ki])
                            n_calls += 1
                        else:
                            old_ki_centroid_dist = np.inf

            else:
                assert(False)

            # Update cluster and centroid_d for elements in affected clusters
            for e in range(m): # NOTE: may be replaced with reverse index + just iterate affected clusters

                # No need to update if centroid
                if is_centroid[e]: continue

                # Cluster of element e
                ce = cluster[e]

                # No need to update if cluster was not affected
                if not is_affected[ce]: continue

                # Check if the distance to element e needs to be computed
                if ki_dist[ce]/2 < centroid_d[e]:
                    distance = dist(elems[e],elems[ki])
                    n_calls += 1

                    if distance < centroid_d[e]:
                        centroid_d[e] = distance
                        cluster[e]    = k



            # Update cluster_radio and size
            cluster_radio = np.zeros(k+1,dtype=np.float64)
            cluster_size = np.zeros(k+1,dtype=np.float64)
            for e in range(m):
                # No need to update if centroid
                if is_centroid[e]: continue

                # Cluster of element e
                ce = cluster[e]

                # Check if updates cluster radious
                cluster_size[ce] += 1
                cluster_radio[ce] = max(cluster_radio[ce],centroid_d[e])

        print(saved1,saved12)

    for i in range(m):
        if is_centroid[i]:
            centroid_d[i] = 0
    # Return centroids and cluster assigns
    return centroids,cluster,n_calls,centroid_d

