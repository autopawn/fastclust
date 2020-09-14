#ifndef SDBS_H
#define SDBS_H

#include <math.h>
#include <assert.h>

typedef long long int lint;

// ===========================================
// ==== DEFINITION OF THE ELEMENTS FROM X ====
// ===========================================
// Note: This code can be used for other types by changing the definition of `elem` and `distance`.

// Elements are vectors
typedef struct {
    int dims;
    double *v;
} elem;

// Distance between elements
static inline double distance(const elem *a, const elem *b){
    assert(a->dims == b->dims);
    double total = 0;
    for(int i=0;i<a->dims;i++){
        total += (a->v[i]-b->v[i])*(a->v[i]-b->v[i]);
    }
    return sqrt(total);
}

// =======================
// ==== RESULT STRUCT ====
// =======================

// Struct to hold the results of a clustering
typedef struct {
    int n; // number of elements
    int k; // number of clusters
    int *clus;    // array of size n with the final cluster of each element.
    double *prox; // array of size n with distances to its centroid for each element.
    int *cents;   // array of size k with indexes of chosen centroids.
    int *is_cent; // array of size n with 0 or 1 values to know if an element is a centroid in O(1).
    lint n_dists; // number of distances computed.
} result;

// Allocate a result in memory
result *result_init(int n, int k);
// Free result from memory
void result_free(result *res);

// ===========================================================
// ==== SIMPLE DIVERSITY BASED STARTING METHOD CLUSTERING ====
// ===========================================================

// All of the following clustering alternatives share the following parameters:
// elems : an array of pointers to elements.
// n     : number of elements
// k     : number of target clusters
// start : initial centroid index

// NOTE: sdbs_original's result doesn't retrieve the prox and cent values.

result *sdbs_original(elem **elems, int n, int k, int start);

result *sdbs_opt1(elem **elems, int n, int k, int start);
result *sdbs_opt2(elem **elems, int n, int k, int start);
result *sdbs_opt3(elem **elems, int n, int k, int start);
result *sdbs_opt4(elem **elems, int n, int k, int start);

#endif
