#include <stdio.h>
#include <stdlib.h>

#include "common.h"

lint clust_opt2_sort_locprox(elem **elems, int n, int k, int start, int *clus, double *prox){
    assert(k>0 && k<=n);
    assert(clus!=NULL);
    assert(prox!=NULL);

    // total number of distances
    lint d_computed = 0;

    // elements on each cluster
    dindxvec **clusters = malloc(sizeof(dindxvec *)*k);
    for(int i=0;i<k;i++) clusters[i] = dindxvec_init(2);
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
            dindxvec_endadd(clusters[0]    ,(dindx){.dist=prox[i],.index=i});
            dindxvec_endadd(clusmov_dmem[i],(dindx){.dist=prox[i],.index=0});
        }
    }
    // Sort cluster elements from nearest to furthest
    qsort(clusters[0]->items,clusters[0]->len,sizeof(dindx),dindx_cmp);

    // Iterate to find new centroids
    for(int h=1;h<k;h++){
        // Pick new centroid
        int pick = -1;
        double pickdist = -1.0; // any negative value should do.
        for(int j=0;j<h;j++){
            if(clusters[j]->len>0){
                dindx rad = dindxvec_endtop(clusters[j]);
                if(rad.dist > pickdist){
                    pickdist = rad.dist;
                    pick = rad.index;
                    assert(clus[pick]==j);
                }
            }
        }
        cents[h] = pick;
        assert(pick!=-1);

        // Remove from current cluster
        dindxvec_endpop(clusters[clus[pick]]);
        clus[pick] = h;

        // Right lower bound
        double cprox_lower_bound_right = prox[pick];

        // Closest centroid known on the left
        dindx di = dindxvec_endpop(clusmov_dmem[pick]);
        double closest_centroid_left_dist = di.dist;
        int closest_centroid_left_j       = di.index;

        // -- Compute distance to old centroids, when needed
        double *cprox = malloc(sizeof(double)*h);
        for(int j=0;j<h;j++) cprox[j] = 1/0.0; // we will mark centroids too far with inf (here, now)

        // Iterate from right to left
        for(int j=h-1;j>=0;j--){
            // Get another left lower bound, if needed
            if(j==closest_centroid_left_j && clusmov_dmem[pick]->len>0){
                cprox[j] = closest_centroid_left_dist;
                dindx di = dindxvec_endpop(clusmov_dmem[pick]);
                closest_centroid_left_dist = di.dist;
                closest_centroid_left_j    = di.index;
                continue;
            }

            // Skip cluster's centroid by right lower bound
            double radious_j = (clusters[j]->len==0)? 0 : dindxvec_endtop(clusters[j]).dist;
            if(radious_j < cprox_lower_bound_right/2) continue;

            // Skip cluster's centroid by left lower bound
            if(j > closest_centroid_left_j && radious_j < (prox[cents[j]] - closest_centroid_left_dist)/2) continue;

            // This will probably be an affected cluster, compute the distance
            cprox[j] = distance(elems[pick],elems[cents[j]]);  d_computed++;

            // Update right lower bound
            double candidate = prox[cents[j]] - cprox[j];
            if(cprox_lower_bound_right<candidate) cprox_lower_bound_right = candidate;
        }

        // Steal pairs
        for(int j=0;j<h;j++){
            // Find elements in cluster j but far from the centroid
            dindxvec *a = dindxvec_init(2);
            while(clusters[j]->len>0 && dindxvec_endtop(clusters[j]).dist > cprox[j]/2){
                dindxvec_endadd(a,dindxvec_endpop(clusters[j]));
            }
            // Iterate them
            while(a->len>0){
                dindx di = dindxvec_endpop(a);
                int i = di.index;

                // Compute distance to new centroid
                double prox2 = distance(elems[pick],elems[i]); d_computed++;

                // Check if element is moved to the new centroid
                if(prox2 < prox[i]){
                    di.dist = prox2;
                    dindxvec_endadd(clusters[h],di);
                    dindxvec_endadd(clusmov_dmem[i],(dindx){.dist=prox2,.index=h});
                    prox[i] = prox2;
                    clus[i] = h;
                }else{
                    dindxvec_endadd(clusters[j],di);
                }
            }
            dindxvec_free(a);
        }

        // Sort cluster elements from nearest to furthest
        qsort(clusters[h]->items,clusters[h]->len,sizeof(dindx),dindx_cmp);

        // Free cprox
        free(cprox);
    }

    // Fix the distances to centroid for the centroids
    for(int h=0;h<k;h++) prox[cents[h]] = 0;

    // -- free memory
    for(int i=0;i<n;i++) dindxvec_free(clusmov_dmem[i]);
    free(clusmov_dmem);
    free(cents);
    for(int i=0;i<k;i++) dindxvec_free(clusters[i]);
    free(clusters);

    // Return assigns
    return d_computed;
}

