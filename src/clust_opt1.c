#include <stdio.h>
#include <stdlib.h>

#include "common.h"

lint clust_opt1(elem **elems, int n, int k, int start, int *clus, double *prox){
    assert(k>0 && k<=n);
    assert(clus!=NULL);
    assert(prox!=NULL);

    // total number of distances
    lint d_computed = 0;

    // centroids
    int *cents = malloc(sizeof(int)*k);

    // initialize
    cents[0] = start;
    clus[start] = 0;
    prox[start] = 1/0.0;
    for(int i=0;i<n;i++){
        if(i!=start){
            prox[i]  = distance(elems[start],elems[i]); d_computed++;
            clus[i] = 0;
        }
    }

    // Iterate to find new centroids
    for(int h=1;h<k;h++){
        // Pick new centroid
        int pick = -1;
        double pickdist = -1.0; // any negative value should do
        for(int i=0;i<n;i++){
            if(cents[clus[i]] != i){ // if this is not a centroid
                if(prox[i] > pickdist){ // if it replaces the current candidate for new centroid
                    pick = i;
                    pickdist = prox[i];
                }
            }
        }

        cents[h] = pick;
        assert(pick!=-1);

        // Move to its own cluster
        clus[pick] = h;

        // Update proxs
        for(int i=0;i<n;i++){
            if(cents[clus[i]] != i){ // if this is not a centroid
                double prox2 = distance(elems[pick],elems[i]); d_computed++;
                if(prox2 < prox[i]){
                    prox[i] = prox2;
                    clus[i] = h;
                }
            }
        }
    }

    // Fix the distances to centroid for the centroids
    for(int h=0;h<k;h++) prox[cents[h]] = 0;

    // -- free memory
    free(cents);

    // Return distances computed
    return d_computed;
}

