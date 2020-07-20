#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "common.h"

#define ARRSIZE(X) (sizeof(X)/sizeof((X)[0]))


int main(int argc, char const *argv[]){
    assert(argc==2);

    srand(atoi(argv[1]));

    const int D[]    = {2,3,5,10,20};
    const int N[]    = {50,500,5000,50000,500000};
    const float KP[] = {1./2,1./4,1./8,1./16,1./32};

    for(int i1=0;i1<ARRSIZE(D);i1++){
        int d = D[i1];
        for(int i2=0;i2<ARRSIZE(N);i2++){
            int n = N[i2];
            for(int i3=0;i3<ARRSIZE(KP);i3++){
                float kp = KP[i3];
                int k = (int)(kp*n+0.5);
                if(k>10000) continue;

                // Generate vectors at random
                elem **xs = malloc(sizeof(elem*)*n);
                for(int i=0;i<n;i++){
                    xs[i] = malloc(sizeof(elem));
                    xs[i]->dims = d;
                    xs[i]->v = malloc(sizeof(xs[i]->v[0])*d);
                    for(int j=0;j<d;j++) xs[i]->v[j] = rand()/(float)RAND_MAX;
                }

                // Start counting time (cpu and elapsed)
                lint ndists;
                clock_t start;
                double time;

                start = clock();
                int *assigns1 = clust_opt1(xs,n,k,0,&ndists);
                time = (double)(clock() - start) / (double)CLOCKS_PER_SEC;
                printf("op1\t%d\t%d\t%d\t%lld\t%lf\n",d,n,k,ndists,time);

                start = clock();
                int *assigns2 = clust_opt2(xs,n,k,0,&ndists);
                time = (double)(clock() - start) / (double)CLOCKS_PER_SEC;
                printf("op2\t%d\t%d\t%d\t%lld\t%lf\n",d,n,k,ndists,time);

                start = clock();
                int *assigns3 = clust_opt3(xs,n,k,0,&ndists);
                time = (double)(clock() - start) / (double)CLOCKS_PER_SEC;
                printf("op3\t%d\t%d\t%d\t%lld\t%lf\n",d,n,k,ndists,time);

                start = clock();
                int *assigns4a = clust_opt4a(xs,n,k,0,&ndists);
                time = (double)(clock() - start) / (double)CLOCKS_PER_SEC;
                printf("op4a\t%d\t%d\t%d\t%lld\t%lf\n",d,n,k,ndists,time);

                start = clock();
                int *assigns4b = clust_opt4b(xs,n,k,0,&ndists);
                time = (double)(clock() - start) / (double)CLOCKS_PER_SEC;
                printf("op4b\t%d\t%d\t%d\t%lld\t%lf\n",d,n,k,ndists,time);

                start = clock();
                int *assigns4ab = clust_opt4ab(xs,n,k,0,&ndists);
                time = (double)(clock() - start) / (double)CLOCKS_PER_SEC;
                printf("op4ab\t%d\t%d\t%d\t%lld\t%lf\n",d,n,k,ndists,time);

                for(int i=0;i<n;i++){
                    assert(assigns1[i]==assigns2[i]);
                    assert(assigns2[i]==assigns3[i]);
                    assert(assigns3[i]==assigns4a[i]);
                    assert(assigns4a[i]==assigns4b[i]);
                    assert(assigns4b[i]==assigns4ab[i]);
                }

                free(assigns4ab);
                free(assigns4b);
                free(assigns4a);
                free(assigns3);
                free(assigns1);

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
