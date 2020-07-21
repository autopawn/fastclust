#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "common.h"

#define ARRSIZE(X) (sizeof(X)/sizeof((X)[0]))


int main(int argc, char const *argv[]){
    assert(argc==2);

    srand(atoi(argv[1]));

    const int D[]    = {16,8,2};
    const int N[]    = {1000000,316227,100000,31622,10000,3162,1000,316,100};
    const float KNUM[] = {1,1,1};
    const float KDEN[] = {2,8,32};
    const lint MEMDISTLIMIT = 100000000; // less than 3GB

    for(int i1=0;i1<ARRSIZE(D);i1++){
        int d = D[i1];
        for(int i2=0;i2<ARRSIZE(N);i2++){
            int n = N[i2];
            for(int i3=0;i3<ARRSIZE(KNUM);i3++){
                int knum = KNUM[i3];
                int kden = KDEN[i3];
                int k = (int)((double)knum*(double)n/(double)kden+0.5);

                // Generate vectors at random
                elem **xs = malloc(sizeof(elem*)*n);
                for(int i=0;i<n;i++){
                    xs[i] = malloc(sizeof(elem));
                    xs[i]->dims = d;
                    xs[i]->v = malloc(sizeof(xs[i]->v[0])*d);
                    for(int j=0;j<d;j++) xs[i]->v[j] = rand()/(float)RAND_MAX;
                }

                // Start counting time (cpu and elapsed)
                clock_t start;
                double time;

                int *clus4b    = malloc(sizeof(int)*n);
                double *prox4b = malloc(sizeof(double)*n);
                start = clock();
                lint ndists4b = clust_opt4b(xs,n,k,0,clus4b,prox4b,MEMDISTLIMIT);
                time = (double)(clock() - start) / (double)CLOCKS_PER_SEC;
                printf("op4b\t%d\t%d\t%d/%d\t%lld\t%.9lf\n",d,n,knum,kden,ndists4b,time);

                int *clus4ab    = malloc(sizeof(int)*n);
                double *prox4ab = malloc(sizeof(double)*n);
                start = clock();
                lint ndists4ab = clust_opt4ab(xs,n,k,0,clus4ab,prox4ab,MEMDISTLIMIT);
                time = (double)(clock() - start) / (double)CLOCKS_PER_SEC;
                printf("op4ab\t%d\t%d\t%d/%d\t%lld\t%.9lf\n",d,n,knum,kden,ndists4ab,time);

                int *clus1    = malloc(sizeof(int)*n);
                double *prox1 = malloc(sizeof(double)*n);
                start = clock();
                lint ndists1 = clust_opt1(xs,n,k,0,clus1,prox1);
                time = (double)(clock() - start) / (double)CLOCKS_PER_SEC;
                printf("op1\t%d\t%d\t%d/%d\t%lld\t%.9lf\n",d,n,knum,kden,ndists1,time);

                int *clus2    = malloc(sizeof(int)*n);
                double *prox2 = malloc(sizeof(double)*n);
                start = clock();
                lint ndists2 = clust_opt2(xs,n,k,0,clus2,prox2);
                time = (double)(clock() - start) / (double)CLOCKS_PER_SEC;
                printf("op2\t%d\t%d\t%d/%d\t%lld\t%.9lf\n",d,n,knum,kden,ndists2,time);

                int *clus3    = malloc(sizeof(int)*n);
                double *prox3 = malloc(sizeof(double)*n);
                start = clock();
                lint ndists3 = clust_opt3(xs,n,k,0,clus3,prox3);
                time = (double)(clock() - start) / (double)CLOCKS_PER_SEC;
                printf("op3\t%d\t%d\t%d/%d\t%lld\t%.9lf\n",d,n,knum,kden,ndists3,time);

                int *clus4a    = malloc(sizeof(int)*n);
                double *prox4a = malloc(sizeof(double)*n);
                start = clock();
                lint ndists4a = clust_opt4a(xs,n,k,0,clus4a,prox4a);
                time = (double)(clock() - start) / (double)CLOCKS_PER_SEC;
                printf("op4a\t%d\t%d\t%d/%d\t%lld\t%.9lf\n",d,n,knum,kden,ndists4a,time);


                for(int i=0;i<n;i++){
                    assert(clus2[i]==clus1[i]);
                    assert(clus3[i]==clus1[i]);
                    assert(clus4a[i]==clus1[i]);
                    assert(clus4b[i]==clus1[i] || ndists4b==-1);
                    assert(clus4ab[i]==clus1[i] || ndists4ab==-1);
                }

                free(prox4a);
                free(clus4a);
                free(prox3);
                free(clus3);
                free(prox2);
                free(clus2);
                free(prox1);
                free(clus1);
                free(prox4ab);
                free(clus4ab);
                free(prox4b);
                free(clus4b);

                // Free memory
                for(int i=0;i<n;i++){
                    free(xs[i]->v);
                    free(xs[i]);
                }
                free(xs);
            }

        }
    }

    return 0;
}
