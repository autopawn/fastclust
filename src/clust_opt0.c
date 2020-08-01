#include <stdio.h>
#include <stdlib.h>

#include "common.h"

lint clust_opt0(elem **elems, int n, int k, int start, int *cents){
    assert(k>0 && k<=n);
    assert(cents!=NULL);

    // total number of distances
    lint d_computed = 0;

    // centroids
    int *is_cent = malloc(sizeof(int)*n);
    for(int j=0;j<n;j++) is_cent[j] = 0;

    // initialize
    cents[0] = start;
    is_cent[start] = 1;

    for(int h=1;h<k;h++){

        double best_dist = -1.0; // any negative value should do.
        int best_j = -1;

        // Search for new centroid
        for(int j=0;j<n;j++){ // j: index of candidate new centroid

            if(is_cent[j]) continue;

            double min_dist = 1/0.0;

            for(int i=0;i<h;i++){
                double dist = distance(elems[cents[i]],elems[j]); d_computed++;
                if(dist<min_dist) min_dist = dist;
            }

            if(min_dist > best_dist){
                best_dist = min_dist;
                best_j = j;
            }

        }
        assert(best_j!=-1);

        // Set new centroid
        cents[h] = best_j;
        is_cent[best_j] = 1;

    }

    free(is_cent);

    // Return distances computed
    return d_computed;
}

