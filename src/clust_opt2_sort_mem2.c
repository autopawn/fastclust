#include <stdio.h>
#include <stdlib.h>

#include "common.h"

lint clust_opt2_sort_mem2(elem **elems, int n, int k, int start, int *clus, double *prox, lint mem_lim){
    assert(k>0 && k<=n);
    assert(clus!=NULL);
    assert(prox!=NULL);

    // total number of distances
    lint d_computed = 0;
    // total number of stored distances
    lint d_stored = 0;

    // elements on each cluster
    dindxvec **clusters = malloc(sizeof(dindxvec *)*k);
    for(int i=0;i<k;i++) clusters[i] = dindxvec_init(2);
    // centroids
    int *cents = malloc(sizeof(int)*k);
    // elem id. -> dindx // stored all distances to the cluster centroids
    dindxvec **clusall_dmem = malloc(sizeof(dindxvec *)*n);
    for(int i=0;i<n;i++) clusall_dmem[i] = dindxvec_init(1);

    // initialize
    cents[0] = start;
    clus[start] = 0;
    prox[start] = 1/0.0;
    for(int i=0;i<n;i++){
        if(i!=start){
            prox[i] = distance(elems[start],elems[i]); d_computed++;
            clus[i] = 0;
            dindxvec_endadd(clusters[0]    ,(dindx){.dist=prox[i],.index=i});
            dindxvec_endadd(clusall_dmem[i],(dindx){.dist=prox[i],.index=0});
        }
    }
    // Sort cluster elements from nearest to furthest
    qsort(clusters[0]->items,clusters[0]->len,sizeof(dindx),dindx_cmp);


    int abort = 0;
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

        // Compute distance to old centroids
        double *cprox = malloc(sizeof(double)*h);
        for(int j=0;j<h;j++) cprox[j] = -1; // Mark as unknown
        // Use already computed distances
        while(clusall_dmem[pick]->len>0){
            dindx dp = dindxvec_endpop(clusall_dmem[pick]);
            cprox[dp.index] = dp.dist;
        }

        for(int j=0;j<h;j++){
            // If distance is marked as uknown, compute it.
            if(cprox[j]<0){
                cprox[j] = distance(elems[pick],elems[cents[j]]);  d_computed++;
            }
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

                int passes = 1;

                // New centroid too far away from current centroid // NOTE: found it is not needed.
                if(prox[i] <= cprox[j]/2){
                    passes = 0;
                }

                if(passes){
                    // Go over all known distances, check if any gives a big upper bound
                    int len = clusall_dmem[i]->len;
                    for(int o=len-1;o>=0;o--){
                        dindx dp = clusall_dmem[i]->items[o];

                        double bound = cprox[dp.index] - dp.dist;
                        if(bound<0) bound = -bound;

                        if(bound >= prox[i]){
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
                    dindxvec_endadd(clusall_dmem[i],(dindx){.dist=prox2,.index=h});  d_stored++;
                    if(d_stored>mem_lim) abort = 1;

                    // Check if element is moved to the new centroid
                    if(prox2 >= prox[i]) passes = 0;
                }

                if(passes){ // Move the element to the new centroid
                    di.dist = prox2;
                    dindxvec_endadd(clusters[h],di);
                    prox[i] = prox2;
                    clus[i] = h;
                }else{
                    dindxvec_endadd(clusters[j],di);
                }

                if(abort) break;
            }
            dindxvec_free(a);

            if(abort) break;
        }

        // Sort cluster elements from nearest to furthest
        qsort(clusters[h]->items,clusters[h]->len,sizeof(dindx),dindx_cmp);

        // Free cprox
        free(cprox);

        // Free clusal_dmem[pick]
        dindxvec_free(clusall_dmem[pick]);
        clusall_dmem[pick] = NULL;

        if(abort) break;
    }

    if(!abort){
        // Fix the distances to centroid for the centroids
        for(int h=0;h<k;h++) prox[cents[h]] = 0;
    }

    // -- free memory
    for(int i=0;i<n;i++){
        if(clusall_dmem[i]!=NULL) dindxvec_free(clusall_dmem[i]);
    }
    free(clusall_dmem);
    free(cents);
    for(int i=0;i<k;i++) dindxvec_free(clusters[i]);
    free(clusters);

    // Return distances computed
    return abort? -1 : d_computed;
}

