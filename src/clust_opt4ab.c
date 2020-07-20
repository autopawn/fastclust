#include <stdio.h>
#include <stdlib.h>

#include "common.h"

// Clustering using optimization 3, retrieves array of assignments
int *clust_opt4ab(elem **elems, int n, int k, int start, lint *dists_computed){
    assert(k>0 && k<=n);
    // total number of distances
    lint d_computed = 0;

    // elem id. -> cluster id.
    int *clus = malloc(sizeof(int)*n);
    // elem id. -> distance to centroid
    double *prox = malloc(sizeof(double)*n);
    // elements on each cluster
    dindxvec **clusters = malloc(sizeof(dindxvec *)*k);
    for(int i=0;i<k;i++) clusters[i] = dindxvec_init(2);
    // centroids
    int *cents = malloc(sizeof(int)*k);
    // elem id. -> dindx // stored computed distances to the cluster (were the element has been) centroids
    dindxvec **clusmov_dmem = malloc(sizeof(dindxvec *)*n);
    for(int i=0;i<n;i++) clusmov_dmem[i] = dindxvec_init(1);
    // elem id. -> dindx // stored all distances to the cluster centroids
    dindxvec **clusall_dmem = malloc(sizeof(dindxvec *)*n);
    for(int i=0;i<n;i++) clusall_dmem[i] = dindxvec_init(1);

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
            dindxvec_endadd(clusall_dmem[i],(dindx){.dist=prox[i],.index=0});
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

        // Compute lower and upper bounds to old centroids
        double *cprox_lo = malloc(sizeof(double)*h);
        double *cprox_hi = malloc(sizeof(double)*h);
        for(int j=0;j<h;j++){
            cprox_lo[j] = 0;
            cprox_hi[j] = 1/0.0; // Actually, it is always inf unless directly computed.
        }
        // Use already computed distances
        while(clusall_dmem[pick]->len>0){
            dindx dp = dindxvec_endpop(clusall_dmem[pick]);
            cprox_lo[dp.index] = dp.dist;
            cprox_hi[dp.index] = dp.dist;
        }

        // Iterate from right to left
        for(int j=h-1;j>=0;j--){
            // Get another left lower bound, if needed
            if(j==closest_centroid_left_j && clusmov_dmem[pick]->len>0){
                cprox_lo[j] = closest_centroid_left_dist;
                cprox_hi[j] = closest_centroid_left_dist;
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
            double radious_j = (clusters[j]->len==0)? 0 : dindxvec_endtop(clusters[j]).dist;
            if(radious_j < cprox_lo[j]/2) continue;

            // This will probably be an affected cluster, compute the distance
            double cprox_j = distance(elems[pick],elems[cents[j]]);  d_computed++;
            cprox_lo[j] = cprox_j;
            cprox_hi[j] = cprox_j;

            // Update right lower bound
            double candidate = prox[cents[j]] - cprox_j;
            if(cprox_lower_bound_right<candidate) cprox_lower_bound_right = candidate;
        }

        // Steal pairs
        for(int j=0;j<h;j++){
            // Find elements in cluster j but far from the centroid
            dindxvec *a = dindxvec_init(2);
            while(clusters[j]->len>0 && dindxvec_endtop(clusters[j]).dist > cprox_lo[j]/2){
                dindxvec_endadd(a,dindxvec_endpop(clusters[j]));
            }
            // Iterate them
            while(a->len>0){
                dindx di = dindxvec_endpop(a);
                int i = di.index;

                int passes = 1;

                // New centroid too far away from current centroid
                if(prox[i] <= cprox_lo[j]/2){
                    passes = 0;
                }

                if(passes){
                    // Go over all known distances, check if any gives a big upper bound
                    int len = clusall_dmem[i]->len;
                    for(int o=len-1;o>=0;o--){
                        dindx dp = clusall_dmem[i]->items[o];

                        double bound1 = cprox_lo[dp.index] - dp.dist;
                        double bound2 = dp.dist - cprox_hi[dp.index];

                        if(bound1 >= prox[i] || bound2 >= prox[i]){
                            // Found a larger lower bound, pick is too far away
                            passes = 0;
                            break;
                        }
                    }
                }

                double prox2;

                if(passes){
                    // Compute distance to new centroid
                    prox2 = distance(elems[pick],elems[i]); d_computed++;

                    // Save the distance that was computed for future use
                    dindxvec_endadd(clusall_dmem[i],(dindx){.dist=prox2,.index=h});

                    // Check if element is moved to the new centroid
                    if(prox2 >= prox[i]) passes = 0;
                }

                if(passes){ // Move the element to the new centroid
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
        free(cprox_hi);
        free(cprox_lo);
    }

    // -- free memory
    for(int i=0;i<n;i++) dindxvec_free(clusall_dmem[i]);
    free(clusall_dmem);

    for(int i=0;i<n;i++) dindxvec_free(clusmov_dmem[i]);
    free(clusmov_dmem);

    free(cents);

    for(int i=0;i<k;i++) dindxvec_free(clusters[i]);
    free(clusters);

    free(prox);
    // Set the distances computed
    if(dists_computed!=NULL) *dists_computed = d_computed;
    // Return assigns
    return clus;
}

