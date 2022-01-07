#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "sdbs.h"

// Number of dimensions for the vectors
#define VEC_DIMS 2

// Number vectors on each element
#define N_VECS_PER_ELEM 8

// Elements of X are sets of vectors
struct elem {
    int index; // Used to break ties
    double v[N_VECS_PER_ELEM][VEC_DIMS];
};

// The distance is the mean geometric error
double distance(const elem *a, const elem *b){
    double total = 0;

    for(int t=0;t<2;t++){
        // For each element of a, find the closest element of b
        for(int i=0;i<N_VECS_PER_ELEM;i++){
            // Vector from a
            const double *a_vec = a->v[i];

            // Closest distance (so far) to an element from b
            double closest = INFINITY;

            for(int j=0;j<N_VECS_PER_ELEM;j++){
                // Vector from a
                const double *b_vec = b->v[j];

                // Compute euclidean distance between vectors
                double euc_distance = 0;
                for(int d=0;d<VEC_DIMS;d++){
                    euc_distance += (a_vec[d]-b_vec[d])*(a_vec[d]-b_vec[d]);
                }
                euc_distance = sqrt(euc_distance);

                // Replace closest distance
                if(euc_distance < closest) closest = euc_distance;
            }

            // Add closest distance to the total distance
            total += closest;
        }

        // Swap a and b for the 2nd iteration
        const elem *c = a;
        a = b;
        b = c;
    }
    return total/(2*N_VECS_PER_ELEM);
}

// Size of an array (if known at compile-time!)
#define ARRSIZE(X) (sizeof(X)/sizeof((X)[0]))

// Predicted number of distances that need to be computed with opt0
// sum_(h=1)^(k - 1) h (n - h) = -1/6 (k - 1) k (2 k - 3 n - 1)
lint predicted_opt0_dists(lint n, lint k){
    return - (k - 1) * k * (2 * k - 3 * n - 1) / 6;
}

int main(int argc, char const *argv[]){
    // Expect second argument as seed.
    assert(argc==2);
    srand(atoi(argv[1]));

    // Strategies that will be tested
    const char *STRATEGIES[] = {"op0","op1","op2","op3","op4"};

    // Values of n
    const int N[]    = {100,316,1000,3162,10000,31622,100000,316227};//,1000000};

    // Numerators and denominators of alpha, where k = alpha * n
    const float KNUM[] = {1,2};//{2,1,1};
    const float KDEN[] = {4,3};//{3,8,32};

    // Limit for the values of K
    const int KLIMIT = 316227;

    // Limit of distances to be evaluated for the original strategy.
    const lint ORIGINALDISTLIMIT = 1000000000;

    for(int i2=0;i2<ARRSIZE(N);i2++){
        int n = N[i2];
        for(int i3=0;i3<ARRSIZE(KNUM);i3++){
            int knum = KNUM[i3];
            int kden = KDEN[i3];
            int k = (int)((double)knum*(double)n/(double)kden+0.5);
            if(k>KLIMIT) continue;

            // Generate n elements at random
            elem **xs = malloc(sizeof(elem*)*n);
            for(int i=0;i<n;i++){
                xs[i] = malloc(sizeof(elem));
                xs[i]->index = i;
                for(int p=0;p<N_VECS_PER_ELEM;p++){
                    for(int d=0;d<VEC_DIMS;d++) xs[i]->v[p][d] = (float)rand()/RAND_MAX;
                }
            }

            // Time counters
            clock_t start;
            double time;

            result *res_old     = NULL;

            for(int i4=0;i4<ARRSIZE(STRATEGIES);i4++){
                const char* strat = STRATEGIES[i4];

                result *res_current = NULL;

                start = clock(); // start counting time.

                if(strcmp("op0",strat)==0){
                    // Don't run optimization 0 if predicted number of distances is larger than OPT0DISTLIMIT
                    lint predicted_ndists = predicted_opt0_dists(n,k);
                    if(predicted_ndists>ORIGINALDISTLIMIT) continue;

                    res_current = sdbs_original(xs,n,k,0);
                    assert(res_current->n_dists==predicted_ndists);

                }else if(strcmp("op1",strat)==0){
                    res_current = sdbs_opt1(xs,n,k,0);
                }else if(strcmp("op2",strat)==0){
                    res_current = sdbs_opt2(xs,n,k,0);
                }else if(strcmp("op3",strat)==0){
                    res_current = sdbs_opt3(xs,n,k,0);
                }else if(strcmp("op4",strat)==0){
                    res_current = sdbs_opt4(xs,n,k,0);
                }else{
                    assert(!"Valid strategy name in STRATEGIES.");
                }

                // Count time
                time = (double)(clock() - start) / (double)CLOCKS_PER_SEC;

                // Print results
                assert(res_current!=NULL);
                printf("%s\t%d\t%d/%d\t%lld\t%.9lf\n",strat,n,knum,kden,res_current->n_dists,time);

                // Ensure that results are the same for each strategy
                if(res_old!=NULL && res_current!=NULL){
                    for(int j=0;j<k;j++){
                        if(res_current->cents[j] != res_old->cents[j]){
                            printf("Warning: %s inconsistent on cents[%d].\n",strat,j);
                            break;
                        }
                    }
                    for(int i=0;i<n;i++){
                        if(res_current->prox[i] != res_old->prox[i]){
                            printf("Warning: %s inconsistent on prox[%d].\n",strat,i);
                            break;
                        }
                    }
                    for(int i=0;i<n;i++){
                        if(res_current->clus[i] != res_old->clus[i]){
                            printf("Warning: %s inconsistent on clus[%d].\n",strat,i);
                            break;
                        }
                    }
                }


                // Free res_old, and swap results
                if(res_old!=NULL) result_free(res_old);
                if(strcmp(strat,"op0")==0){ // just free if strat is original
                    result_free(res_current);
                }else{
                    res_old = res_current;
                }
                res_current = NULL;

            }

            // Free res_old
            if(res_old!=NULL) result_free(res_old);

            // Free elements memory
            for(int i=0;i<n;i++) free(xs[i]);
            free(xs);

        }
    }
}
