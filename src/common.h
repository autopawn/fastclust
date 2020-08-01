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
    // Decreasing index order so that results are the same in
    // case of distance ties for the different optimizations.
    return bb->index - aa->index;
}

// Define dindxvec as a vector of dindx
CADTS_VECTOR(dindxvec,dindx)


/* Simple diversity-based clustering initialization method:
elems : an array of pointers to elements.
n     : number of elements
k     : number of clusters
start : initial centroid index
cents : size k array to store centroids
*/
lint clust_opt0(elem **elems, int n, int k, int start, int *cents);

/* All of the following clustering alternatives share the following parameters:
elems : an array of pointers to elements.
n     : number of elements
k     : number of clusters
start : initial centroid index
clus  : size n array to store final cluster of each element
prox  : size n array to store the distance of each element to the centroid
(some) mem_lim : limit for stored distances, if reached; the operation is aborted and -1 is returned.
--
Returns: the number of distances that had to be computed, or a negative value if the computation could not finish.
*/

lint clust_opt1(elem **elems, int n, int k, int start, int *clus, double *prox);
lint clust_opt2(elem **elems, int n, int k, int start, int *clus, double *prox);
lint clust_opt2_mem1(elem **elems, int n, int k, int start, int *clus, double *prox);
lint clust_opt2_locprox_mem1(elem **elems, int n, int k, int start, int *clus, double *prox);

lint clust_opt2_sort(elem **elems, int n, int k, int start, int *clus, double *prox);
lint clust_opt2_sort_locprox(elem **elems, int n, int k, int start, int *clus, double *prox);

lint clust_opt2_sort_mem2(elem **elems, int n, int k, int start, int *clus, double *prox, lint mem_lim);
lint clust_opt2_sort_locprox_mem2(elem **elems, int n, int k, int start, int *clus, double *prox, lint mem_lim);

#endif
