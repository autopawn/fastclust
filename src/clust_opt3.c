#include <stdio.h>
#include <stdlib.h>

#include "common.h"

// Clustering using optimization 3, retrieves array of assignments
int *clust_opt3(elem **elems, int n, int k, int start, lint *dists_computed){
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

    // initialize
    cents[0] = start;
    clus[start] = 0;
    prox[start] = 1/0.0;
    for(int i=0;i<n;i++){
        if(i!=start){
            prox[i]  = distance(elems[start],elems[i]); d_computed++;
            clus[i] = 0;
            dindx di;
            di.dist  = prox[i];
            di.index = i;
            dindxvec_endadd(clusters[0],di);
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

        // Compute distance to old centroids
        double *cprox = malloc(sizeof(double)*h);
        for(int j=0;j<h;j++){
            cprox[j] = distance(elems[pick],elems[cents[j]]);  d_computed++;
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

    // -- free memory
    free(cents);

    for(int i=0;i<k;i++) dindxvec_free(clusters[i]);
    free(clusters);

    free(prox);
    // Set the distances computed
    if(dists_computed!=NULL) *dists_computed = d_computed;
    // Return assigns
    return clus;
}

