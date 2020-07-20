#ifndef COMMON_H
#define COMMON_H

#include "cadts_vector.h"
#include <math.h>

typedef long long int lint;

// Elements will be arrays
typedef struct {
    double *v;
    int dims;
} elem;

// Distance between elements
static inline double distance(elem *a, elem *b){
    assert(a->dims == b->dims);
    double total = 0;
    for(int i=0;i<a->dims;i++){
        total += (a->v[i]-b->v[i])*(a->v[i]-b->v[i]);
    }
    return sqrt(total);
}

// distance-index pair
typedef struct {
    double dist;
    int index;
} dindx;

// cmp function to compare distance-index pairs
static inline int dindx_cmp(const void *a, const void *b){
    dindx *aa = (dindx *) a;
    dindx *bb = (dindx *) b;
    if(aa->dist < bb->dist) return -1;
    if(aa->dist > bb->dist) return +1;
    return aa->index - bb->index;
}

// Define dindxvec as a vector of dindx
CADTS_VECTOR(dindxvec,dindx)


int *clust_opt1(elem **elems, int n, int k, int start, lint *dists_computed);
int *clust_opt2(elem **elems, int n, int k, int start, lint *dists_computed);
int *clust_opt3(elem **elems, int n, int k, int start, lint *dists_computed);
int *clust_opt4a(elem **elems, int n, int k, int start, lint *dists_computed);
int *clust_opt4b(elem **elems, int n, int k, int start, lint *dists_computed);
int *clust_opt4ab(elem **elems, int n, int k, int start, lint *dists_computed);

#endif
