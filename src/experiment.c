#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "common.h"

#define ARRSIZE(X) (sizeof(X)/sizeof((X)[0]))

// Predicted number of distances that need to be computed with opt0
// https://www.wolframalpha.com/input/?i=%5Csum_%7Bh%3D1%7D%5E%7Bk-1%7D+%28n-h%29+h
// sum_(h=1)^(k - 1) h (n - h) = -1/6 (k - 1) k (2 k - 3 n - 1)
lint predicted_opt0_dists(lint n, lint k){
    return - (k - 1) * k * (2 * k - 3 * n - 1) / 6;
}

int main(int argc, char const *argv[]){
    assert(argc==2);

    srand(atoi(argv[1]));

    // const char *STRATEGIES[] = {"op0","op1","op2","op2_m1","op2_lo_m1","op2_s","op2_s_lo","op2_s_m2","op2_s_lo_m2"};
    const char *STRATEGIES[] = {"op0"};

    const int D[]    = {2,8,16};

    const int N[]    = {100,316,1000,3162,10000,31622,100000,316227,1000000};
    // const int N[]    = {100,316,1000};

    const float KNUM[] = {2,1,1};
    const float KDEN[] = {3,8,32};


    const lint MEMDISTLIMIT  = 100000000; // less than 3GB
    const lint OPT0DISTLIMIT = 1000000000;

    const int KLIMIT = 316227;

    for(int i1=0;i1<ARRSIZE(D);i1++){
        int d = D[i1];
        for(int i2=0;i2<ARRSIZE(N);i2++){
            int n = N[i2];
            for(int i3=0;i3<ARRSIZE(KNUM);i3++){
                int knum = KNUM[i3];
                int kden = KDEN[i3];
                int k = (int)((double)knum*(double)n/(double)kden+0.5);
                if(k>KLIMIT) continue;

                // Generate vectors at random
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

                // Previous results
                int *clus1    = NULL;
                double *prox1 = NULL;

                for(int i4=0;i4<ARRSIZE(STRATEGIES);i4++){
                    const char* strat = STRATEGIES[i4];

                    int *clus    = malloc(sizeof(int)*n);
                    double *prox = malloc(sizeof(double)*n);

                    start = clock();

                    lint ndists;
                    if(strcmp("op0",strat)==0){
                        // Run optimization 0, if predicted number of distances isn't larger than OPT0DISTLIMIT
                        lint predicted_ndists = predicted_opt0_dists(n,k);

                        if(predicted_ndists<=OPT0DISTLIMIT){
                            int *cents = malloc(sizeof(int)*k);
                            ndists = clust_opt0(xs,n,k,0,cents);
                            assert(ndists==predicted_ndists);
                            free(cents);
                        }else{
                            ndists = -1;
                        }

                    }else if(strcmp("op1",strat)==0){
                        ndists = clust_opt1(xs,n,k,0,clus,prox);

                    }else if(strcmp("op2",strat)==0){
                        ndists = clust_opt2(xs,n,k,0,clus,prox);

                    }else if(strcmp("op2_m1",strat)==0){
                        ndists = clust_opt2_mem1(xs,n,k,0,clus,prox);

                    }else if(strcmp("op2_lo_m1",strat)==0){
                        ndists = clust_opt2_locprox_mem1(xs,n,k,0,clus,prox);

                    }else if(strcmp("op2_s",strat)==0){
                        ndists = clust_opt2_sort(xs,n,k,0,clus,prox);

                    }else if(strcmp("op2_s_lo",strat)==0){
                        ndists = clust_opt2_sort_locprox(xs,n,k,0,clus,prox);

                    }else if(strcmp("op2_s_m2",strat)==0){
                        ndists = clust_opt2_sort_mem2(xs,n,k,0,clus,prox,MEMDISTLIMIT);

                    }else if(strcmp("op2_s_lo_m2",strat)==0){
                        ndists = clust_opt2_sort_locprox_mem2(xs,n,k,0,clus,prox,MEMDISTLIMIT);

                    }else{
                        assert(!"Valid strategy name in STRATEGIES.");

                    }

                    time = (double)(clock() - start) / (double)CLOCKS_PER_SEC;

                    printf("%s\t%d\t%d\t%d/%d\t%lld\t%.9lf\n",strat,d,n,knum,kden,ndists,time);

                    if(prox1!=NULL){
                        // Ensure same results than previous optimization.
                        for(int i=0;i<n;i++) assert(clus[i]==clus1[i]);
                        free(prox1);
                        free(clus1);
                    }

                    // Swap results
                    prox1 = prox;
                    clus1 = clus;
                }

                // Free previous result memory
                if(prox1!=NULL) free(prox1);
                if(clus1!=NULL) free(clus1);

                // Free elements memory
                for(int i=0;i<n;i++){
                    free(xs[i]->v);
                    free(xs[i]);
                }
                free(xs);

                printf("\n");
            }

        }
    }

    return 0;
}
