#include <stdlib.h>
#include <stdio.h>

#include "cadts_vector.h"

#include "sdbs.h"

// Memory is a vector of (double,int) pairs
typedef struct {
    double dist;
    int indx;
} pair;

// Define memory as a vector of pairs, for opt3 and opt4
CADTS_VECTOR(memory,pair)

// Maximum between two doubles
double dmax(double a, double b){
    return (a>b)? a : b;
}
// Absolute value of a double
double dabs(double a){
    return (a<0)? -a : a;
}

result *result_init(int n, int k){
    result *res = malloc(sizeof(result));
    res->n = n;
    res->k = k;

    res->clus = malloc(sizeof(int)*n);
    res->prox = malloc(sizeof(double)*n);
    res->cents = malloc(sizeof(int)*k);
    res->is_cent = malloc(sizeof(int)*n);
    for(int i=0;i<n;i++) res->is_cent[i] = 0;

    res->n_dists = 0;
    return res;
}

void result_free(result *res){
    free(res->is_cent);
    free(res->cents);
    free(res->prox);
    free(res->clus);
    free(res);
}

result *sdbs_original(elem **elems, int n, int k, int start){
    assert(0<=start && start<n);
    assert(k<=n);

    result *res = result_init(n,k);

    // Pick start as first centroid
    res->cents[0] = start;
    res->is_cent[start] = 1;

    // Iterate until finding all centroids
    for(int h=1;h<k;h++){
        int new_cent = -1;
        double max_min_dist = -INFINITY;

        // Find new centroid
        for(int i=0;i<n;i++){
            if(res->is_cent[i]) continue; // skip elements already chosen as centroids.

            // Find minimum distance to current centroids
            double min_dist = INFINITY;
            for(int j=0;j<h;j++){
                int cent = res->cents[j];
                double dist = distance(elems[i],elems[cent]); res->n_dists++;
                if(dist < min_dist) min_dist = dist;
            }

            // Replace current new centroid candidate, if it has the maximum min_dist so far
            if(min_dist > max_min_dist){
                max_min_dist = min_dist;
                new_cent = i;
            }
        }

        // Add new centroid
        assert(new_cent>-1);
        res->cents[h] = new_cent;
        res->is_cent[new_cent] = 1;

    }

    // opt0 does not retrieve clus nor prox
    for(int i=0;i<n;i++){
        res->clus[i] = -1;
        res->prox[i] = INFINITY;
    }

    return res;
}

result *sdbs_opt1(elem **elems, int n, int k, int start){
    assert(0<=start && start<n);
    assert(k<=n);

    result *res = result_init(n,k);

    // Pick start as first centroid
    res->cents[0] = start;
    res->is_cent[start] = 1;
    res->clus[start] = 0;
    res->prox[start] = 0;

    // Compute initial values of prox and move to initial cluster
    for(int i=0;i<n;i++){
        if(res->is_cent[i]) continue;

        res->clus[i] = 0;
        res->prox[i] = distance(elems[i],elems[start]); res->n_dists++;
    }

    // Iterate until finding all centroids
    for(int h=1;h<k;h++){

        // Pick new centroid using prox
        int new_cent = -1;
        double max_prox = -INFINITY;

        for(int i=0;i<n;i++){
            if(res->is_cent[i]) continue; // skip elements already chosen as centroids.

            if(res->prox[i] > max_prox){
                max_prox = res->prox[i];
                new_cent = i;
            }
        }

        // Add new centroid
        assert(new_cent>-1);
        res->cents[h] = new_cent;
        res->is_cent[new_cent] = 1;
        res->clus[new_cent] = h;
        res->prox[new_cent] = 0;

        // Move elements to new centroid, update proximities
        for(int i=0;i<n;i++){
            if(res->is_cent[i]) continue;  // skip elements already chosen as centroids.

            double prox2 = distance(elems[i],elems[new_cent]); res->n_dists++;
            if(prox2 < res->prox[i]){
                res->clus[i] = h;
                res->prox[i] = prox2;
            }
        }
    }

    return res;
}

result *sdbs_opt2(elem **elems, int n, int k, int start){
    assert(0<=start && start<n);
    assert(k<=n);

    result *res = result_init(n,k);

    // Array for centroid-centroid distances
    double *cprox = malloc(sizeof(double)*k);

    // Pick start as first centroid
    res->cents[0] = start;
    res->is_cent[start] = 1;
    res->clus[start] = 0;
    res->prox[start] = 0;

    // Compute initial values of prox and move to initial cluster
    for(int i=0;i<n;i++){
        if(res->is_cent[i]) continue;

        res->clus[i] = 0;
        res->prox[i] = distance(elems[i],elems[start]); res->n_dists++;
    }

    // Iterate until finding all centroids
    for(int h=1;h<k;h++){

        // Pick new centroid using prox
        int new_cent = -1;
        double max_prox = -INFINITY;

        for(int i=0;i<n;i++){
            if(res->is_cent[i]) continue; // skip elements already chosen as centroids.

            if(res->prox[i] > max_prox){
                max_prox = res->prox[i];
                new_cent = i;
            }
        }

        // Add new centroid
        assert(new_cent>-1);
        res->cents[h] = new_cent;
        res->is_cent[new_cent] = 1;
        res->clus[new_cent] = h;
        res->prox[new_cent] = 0;

        // Compute centroid-centroid distances
        for(int j=0;j<h;j++){
            int cent = res->cents[j];
            cprox[j] = distance(elems[cent],elems[new_cent]); res->n_dists++;
        }

        // Move elements to new centroid, update proximities
        for(int i=0;i<n;i++){
            if(res->is_cent[i]) continue;  // skip elements already chosen as centroids.

            if(res->prox[i] > cprox[res->clus[i]]/2){ // element is far enough from its current centroid

                double prox2 = distance(elems[i],elems[new_cent]); res->n_dists++;
                if(prox2 < res->prox[i]){
                    res->clus[i] = h;
                    res->prox[i] = prox2;
                }

            }
        }
    }

    free(cprox);
    return res;
}

result *sdbs_opt3(elem **elems, int n, int k, int start){
    assert(0<=start && start<n);
    assert(k<=n);

    result *res = result_init(n,k);

    // Array for centroid-centroid distances
    double *cprox = malloc(sizeof(double)*k);

    // Memories for each element
    memory **mems = malloc(sizeof(memory *)*n);
    for(int i=0;i<n;i++){
        mems[i] = memory_init(1);
    }

    // Pick start as first centroid
    res->cents[0] = start;
    res->is_cent[start] = 1;
    res->clus[start] = 0;
    res->prox[start] = 0;
    memory_free(mems[start]); // release picked centroid memory.

    // Compute initial values of prox and move to initial cluster
    for(int i=0;i<n;i++){
        if(res->is_cent[i]) continue;

        res->clus[i] = 0;
        res->prox[i] = distance(elems[i],elems[start]); res->n_dists++;

        memory_endadd(mems[i],(pair){.dist=res->prox[i],.indx=0}); // save new prox and clus in memory
    }

    // Iterate until finding all centroids
    for(int h=1;h<k;h++){

        // Pick new centroid using prox
        int new_cent = -1;
        double max_prox = -INFINITY;

        for(int i=0;i<n;i++){
            if(res->is_cent[i]) continue; // skip elements already chosen as centroids.

            if(res->prox[i] > max_prox){
                max_prox = res->prox[i];
                new_cent = i;
            }
        }

        // Add new centroid
        assert(new_cent>-1);
        res->cents[h] = new_cent;
        res->is_cent[new_cent] = 1;
        res->clus[new_cent] = h;
        res->prox[new_cent] = 0;

        // Compute centroid-centroid distances, use memory
        for(int j=0;j<h;j++) cprox[j] = INFINITY;
        int len = memory_len(mems[new_cent]);
        for(int o=0;o<len;o++){
            pair pa = mems[new_cent]->items[o];
            cprox[pa.indx] = pa.dist;
        }
        for(int j=0;j<h;j++){
            if(!isinf(cprox[j])) continue;
            int cent = res->cents[j];
            cprox[j] = distance(elems[cent],elems[new_cent]); res->n_dists++;
        }

        // Move elements to new centroid, update proximities
        for(int i=0;i<n;i++){
            if(res->is_cent[i]) continue;  // skip elements already chosen as centroids.

            int skip = 0; // if this element will be skipped as it won't move to the new cluster

            /* Traverse memory from newest to older, to see if some pair gives a large
            enough upper bound to discard movement */
            int len = memory_len(mems[i]);
            for(int o=len-1;o>=0;o--){
                pair pa = mems[i]->items[o];

                // If the lower bound is larger than the current value of prox, skip this element.
                if(dabs(pa.dist - cprox[pa.indx]) >= res->prox[i]){
                    skip = 1;
                    break;
                }
            }
            if(skip) continue;

            // Compute distance can't be avoided, check if element moves to the new centroid
            double prox2 = distance(elems[i],elems[new_cent]); res->n_dists++;
            if(prox2 < res->prox[i]){
                res->clus[i] = h;
                res->prox[i] = prox2;

                memory_endadd(mems[i],(pair){.dist=prox2,.indx=h}); // save new prox and clus in memory
            }

        }

        // Release picked centroid memory
        memory_free(mems[new_cent]);
    }

    // Free memories
    for(int i=0;i<n;i++){
        if(!res->is_cent[i]) memory_free(mems[i]);
    }
    free(mems);

    free(cprox);
    return res;
}

result *sdbs_opt4(elem **elems, int n, int k, int start){
    assert(0<=start && start<n);
    assert(k<=n);

    result *res = result_init(n,k);

    // Array for centroid-centroid distances
    double *cprox = malloc(sizeof(double)*k);

    // Memories for each element
    memory **mems = malloc(sizeof(memory *)*n);
    for(int i=0;i<n;i++){
        mems[i] = memory_init(1);
    }

    // Array for radii
    double *radii = malloc(sizeof(double)*k);

    // Array for last distance to a centroid
    double *lprox = malloc(sizeof(double)*k);

    // Cprox lower bound array
    double *lower = malloc(sizeof(double)*k);

    // Pick start as first centroid
    res->cents[0] = start;
    res->is_cent[start] = 1;
    res->clus[start] = 0;
    res->prox[start] = 0;
    lprox[0] = INFINITY;
    memory_free(mems[start]); // release picked centroid memory.

    // Compute initial values of prox and move to initial cluster
    for(int i=0;i<n;i++){
        if(res->is_cent[i]) continue;

        res->clus[i] = 0;
        res->prox[i] = distance(elems[i],elems[start]); res->n_dists++;

        memory_endadd(mems[i],(pair){.dist=res->prox[i],.indx=0}); // save new prox and clus in memory
    }

    // Iterate until finding all centroids
    for(int h=1;h<k;h++){

        // Pick new centroid using prox
        int new_cent = -1;
        double max_prox = -INFINITY;

        for(int i=0;i<n;i++){
            if(res->is_cent[i]) continue; // skip elements already chosen as centroids.

            if(res->prox[i] > max_prox){
                max_prox = res->prox[i];
                new_cent = i;
            }
        }

        // Add new centroid
        assert(new_cent>-1);
        res->cents[h] = new_cent;
        res->is_cent[new_cent] = 1;
        res->clus[new_cent] = h;
        lprox[h] = res->prox[new_cent]; // save prox in lprox before setting it to 0
        res->prox[new_cent] = 0;

        { // Compute radii
            for(int j=0;j<=h;j++) radii[j] = 0;
            for(int i=0;i<n;i++){
                // Increase cluster radious if prox is larger
                radii[res->clus[i]] = dmax(res->prox[i],radii[res->clus[i]]);
            }
        }

        { // Compute centroid-centroid distances
            double R = lprox[h];
            pair pa = memory_endpop(mems[new_cent]);
            double L = pa.dist;
            int l = pa.indx;

            for(int j=h-1;j>=0;j--){
                cprox[j] = INFINITY;
                lower[j] = -INFINITY; // shouldn't be used unless cprox[j] stays as inf.

                if(j==l){
                    cprox[j] = L;
                    // Get next valid L from memory
                    assert((memory_len(mems[new_cent])==0) == (j==0));
                    if(memory_len(mems[new_cent])>0){
                        pa = memory_endpop(mems[new_cent]);
                        L = pa.dist;
                        l = pa.indx;
                    }
                }else{
                    // Check if computing cprox[j] can be avoided
                    lower[j] = dmax(R, lprox[j] - L);
                    if(radii[j] <= lower[j]/2) continue;

                    // Compute distance, as cluster may be affected
                    int cent = res->cents[j];
                    cprox[j] = distance(elems[cent],elems[new_cent]); res->n_dists++;
                }

                R = dmax(R,lprox[j] - cprox[j]);
            }
        }

        // Move elements to new centroid, update proximities
        for(int i=0;i<n;i++){
            if(res->is_cent[i]) continue;  // skip elements already chosen as centroids.

            int skip = 0; // if this element will be skipped as it won't move to the new cluster

            /* Traverse memory from newest to older, to see if some pair gives a large
            enough upper bound to discard movement */
            int len = memory_len(mems[i]);
            for(int o=len-1;o>=0;o--){
                pair pa = mems[i]->items[o];

                // Get lower bound
                double bound;
                if(isinf(cprox[pa.indx])){ // if cprox wasn't computed
                    bound = lower[pa.indx] - pa.dist;
                }else{
                    bound = dabs(pa.dist - cprox[pa.indx]);
                }

                // If the lower bound is larger than the current value of prox, skip this element.
                if(bound >= res->prox[i]){
                    skip = 1;
                    break;
                }
            }
            if(skip) continue;

            // Computing distance can't be avoided, check if element moves to the new centroid
            double prox2 = distance(elems[i],elems[new_cent]); res->n_dists++;
            if(prox2 < res->prox[i]){
                res->clus[i] = h;
                res->prox[i] = prox2;

                memory_endadd(mems[i],(pair){.dist=prox2,.indx=h}); // save new prox and clus in memory
            }

        }

        // Release picked centroid memory
        memory_free(mems[new_cent]);
    }

    // Free memory
    free(lower);
    free(lprox);
    free(radii);
    for(int i=0;i<n;i++){
        if(!res->is_cent[i]) memory_free(mems[i]);
    }
    free(mems);

    free(cprox);
    return res;
}

