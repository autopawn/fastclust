#include <stdio.h>
#include <stdlib.h>

#include "common.h"

// Clustering using optimization 3, retrieves array of assignments
int *clust_opt1(elem **elems, int n, int k, int start, lint *dists_computed){
    assert(k>0 && k<=n);
    // total number of distances
    lint d_computed = 0;

    // elem id. -> cluster id.
    int *clus = malloc(sizeof(int)*n);
    // elem id. -> distance to centroid
    double *prox = malloc(sizeof(double)*n);
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
        prox[pick] = 0; /* NOTE: line not actually necessary unless we want to modify
        the code to retrieve distances to centroids. This line is not present in further
        optimizations because the previous value of prox is used for the centroids. */

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

    // -- free memory
    free(cents);
    free(prox);
    // Set the distances computed
    if(dists_computed!=NULL) *dists_computed = d_computed;
    // Return assigns
    return clus;
}

