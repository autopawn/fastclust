#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "sdbs.h"

#define ARRSIZE(X) (sizeof(X)/sizeof((X)[0]))

// Predicted number of distances that need to be computed with opt0
// https://www.wolframalpha.com/input/?i=%5Csum_%7Bh%3D1%7D%5E%7Bk-1%7D+%28n-h%29+h
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

    // Dimensions
    const int DIMS[] = {2,8,16};

    // Values of n
    const int N[]    = {100,316,1000,3162,10000,31622,100000,316227,1000000};

    // Numerators and denominators of alpha, where k = alpha * n
    const float KNUM[] = {2,1,1};
    const float KDEN[] = {3,8,32};

    // Limit for the values of K
    const int KLIMIT = 316227;

    // Limit of distances to be evaluated for the original strategy.
    const lint ORIGINALDISTLIMIT = 100000000000;

    for(int i1=0;i1<ARRSIZE(DIMS);i1++){
        int d = DIMS[i1];
        for(int i2=0;i2<ARRSIZE(N);i2++){
            int n = N[i2];
            for(int i3=0;i3<ARRSIZE(KNUM);i3++){
                int knum = KNUM[i3];
                int kden = KDEN[i3];
                int k = (int)((double)knum*(double)n/(double)kden+0.5);
                if(k>KLIMIT) continue;

                // Generate n vectors at random
                elem **xs = malloc(sizeof(elem*)*n);
                for(int i=0;i<n;i++){
                    xs[i] = malloc(sizeof(elem));
                    xs[i]->dims = d;
                    xs[i]->v = malloc(sizeof(xs[i]->v[0])*d);
                    for(int j=0;j<d;j++) xs[i]->v[j] = rand()/(float)RAND_MAX;
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
                    printf("%s\t%d\t%d\t%d/%d\t%lld\t%.9lf\n",strat,d,n,knum,kden,res_current->n_dists,time);

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
                for(int i=0;i<n;i++){
                    free(xs[i]->v);
                    free(xs[i]);
                }
                free(xs);

            }
        }
    }
}
