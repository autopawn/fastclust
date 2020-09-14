#include <stdio.h>
#include <stdlib.h>

#include "common.h"

lint clust_opt2_locprox_mem1(elem **elems, int n, int k, int start, int *clus, double *prox){
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

        // -- Compute distance to old centroids, when needed
        double *cprox_lo = malloc(sizeof(double)*h); // lower bound for distance to centroid
        double *cprox_hi = malloc(sizeof(double)*h); // upper bound for distance to centroid, can only be inf or the actual value
        for(int j=0;j<h;j++) cprox_lo[j] = 0;
        for(int j=0;j<h;j++) cprox_hi[j] = 1.0/0.0;

        {{
            // Get the radious of each cluster
            double *rad = malloc(sizeof(double)*h);
            for(int j=0;j<h;j++) rad[j] = 0;
            for(int i=0;i<n;i++){
                if(cents[clus[i]]==i) continue; // skip if centroid
                if(rad[clus[i]] < prox[i]) rad[clus[i]] = prox[i];
            }


            // Closest centroid known on the left
            dindx di = dindxvec_endpop(clusmov_dmem[pick]);
            double closest_centroid_left_dist = di.dist;
            int closest_centroid_left_j       = di.index;

            // Right lower bound
            double cprox_lower_bound_right = prox[pick];
            assert(prox[pick]==di.dist);

            // Iterate from right to left
            for(int j=h-1;j>=0;j--){
                // Get another left lower bound, if needed
                if(j==closest_centroid_left_j && clusmov_dmem[pick]->len>0){
                    cprox_lo[j] = closest_centroid_left_dist;
                    cprox_hi[j] = closest_centroid_left_dist;

                    // Update right lower bound
                    double candidate = prox[cents[j]] - cprox_lo[j];
                    if(cprox_lower_bound_right<candidate) cprox_lower_bound_right = candidate;

                    dindx di = dindxvec_endpop(clusmov_dmem[pick]);
                    closest_centroid_left_dist = di.dist;
                    closest_centroid_left_j    = di.index;

                    continue;
                }

                // Update lower bound with right lower bound
                if(cprox_lo[j] < cprox_lower_bound_right) cprox_lo[j] = cprox_lower_bound_right;
                // Update lower bound with left lower bound
                if(j > closest_centroid_left_j){
                    double cprox_lower_bound_left = prox[cents[j]] - closest_centroid_left_dist;
                    if(cprox_lo[j] < cprox_lower_bound_left) cprox_lo[j] = cprox_lower_bound_left;
                }

                // Skip cluster's centroid by lower bound
                if(rad[j] < cprox_lo[j]/2) continue;

                // This will probably be an affected cluster, compute the distance
                double cprox_j = distance(elems[pick],elems[cents[j]]);  d_computed++;
                cprox_lo[j] = cprox_j;
                cprox_hi[j] = cprox_j;

                // Update right lower bound
                double candidate = prox[cents[j]] - cprox_j;
                if(cprox_lower_bound_right<candidate) cprox_lower_bound_right = candidate;
            }

            free(rad);
        }}

        // Steal pairs
        for(int i=0;i<n;i++){
            // Skip if centroid
            if(cents[clus[i]]==i) continue;

            int passes = 1;

            // Go over all known distances, check if any gives a big upper bound
            int len = clusmov_dmem[i]->len;
            for(int o=len-1;o>=0;o--){
                dindx dp = clusmov_dmem[i]->items[o];

                double bound1 = cprox_lo[dp.index] - dp.dist;
                double bound2 = dp.dist - cprox_hi[dp.index];

                if(bound1 >= prox[i] || bound2 >= prox[i]){
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

        free(cprox_lo);
        free(cprox_hi);
    }

    // Fix the distances to centroid for the centroids
    for(int i=0;i<n;i++) dindxvec_free(clusmov_dmem[i]);
    free(clusmov_dmem);
    for(int h=0;h<k;h++) prox[cents[h]] = 0;

    // -- free memory
    free(cents);

    // Return distances computed
    return d_computed;
}

