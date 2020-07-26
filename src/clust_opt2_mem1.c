#include <stdio.h>
#include <stdlib.h>

#include "common.h"

lint clust_opt2_mem1(elem **elems, int n, int k, int start, int *clus, double *prox){
    assert(k>0 && k<=n);
    assert(clus!=NULL);
    assert(prox!=NULL);

    // total number of distances
    lint d_computed = 0;

    // centroids
    int *cents = malloc(sizeof(int)*k);

    // elem id. -> dindx // stored computed distances to the cluster (were the element has been) centroids
    dindxvec **clusmov_dmem = malloc(sizeof(dindxvec *)*n);
    for(int i=0;i<n;i++) clusmov_dmem[i] = dindxvec_init(1);

    // initialize
    cents[0] = start;
    clus[start] = 0;
    prox[start] = 1/0.0;
    for(int i=0;i<n;i++){
        if(i!=start){
            prox[i]  = distance(elems[start],elems[i]); d_computed++;
            clus[i] = 0;
            dindxvec_endadd(clusmov_dmem[i],(dindx){.dist=prox[i],.index=0});
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

        // Compute distance to old centroids
        double *cprox = malloc(sizeof(double)*h);
        for(int j=0;j<h;j++){
            cprox[j] = distance(elems[pick],elems[cents[j]]);  d_computed++;
        }

        // Steal pairs
        for(int i=0;i<n;i++){
            // Skip if centroid
            if(cents[clus[i]]==i) continue;

            int passes = 1;

            // Go over all known distances, check if any gives a big upper bound
            int len = clusmov_dmem[i]->len;
            for(int o=len-1;o>=0;o--){
                dindx dp = clusmov_dmem[i]->items[o];

                double bound = cprox[dp.index] - dp.dist;
                if(bound<0) bound = -bound;

                if(bound >= prox[i]){
                    // Found a larger lower bound, pick is too far away
                    passes = 0;
                    break;
                }
            }

            if(passes){ // pick still has a chance of stealing this element

                // Compute distance to new centroid
                double prox2 = distance(elems[pick],elems[i]); d_computed++;

                if(prox2 < prox[i]){ // Move the element to the new centroid
                    prox[i] = prox2;
                    clus[i] = h;
                    // Save the distance that was computed for future use
                    dindxvec_endadd(clusmov_dmem[i],(dindx){.dist=prox2,.index=h});
                }
            }
        }

        free(cprox);
    }

    // Fix the distances to centroid for the centroids
    for(int h=0;h<k;h++) prox[cents[h]] = 0;

    // -- free memory
    for(int i=0;i<n;i++) dindxvec_free(clusmov_dmem[i]);
    free(clusmov_dmem);
    free(cents);

    // Return distances computed
    return d_computed;
}

