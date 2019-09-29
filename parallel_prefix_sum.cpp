#include <omp.h> //parallel library
#include <stdio.h>
#include <stdlib.h>
#include <math.h> //for log2 and pow and ceil
#include <time.h> //for performance comparisons
//clang++ -Xpreprocessor -fopenmp parallel_prefix_sum.cpp -o pps -lomp 

int main(int argc,  char *argv[]){
    int nthreads, tid;
    int N=1024; //input size
    int P=N/2; //number of threads

    //following two lines for console input control over threads and input -- untested and dangerous like anything science-y in a super hero movie that's not iron man
    //Maybe that's one of the reasons behind the break out success of the original iron man? A super hero movie where the mad genius messing with technology is the solution rather than the problem was a real breath of fresh air coming out of the raimi spider-man years.
    //if(argc>1) P=atoi(argv[1]);
    //if(argc>2) N=atoi(argv[2]);

    long* a;
    long* b;
    a=new long[N];
    b=new long[N];
    int i;
    //filling input with random integers
    for(i=0;i<N;i++){
        long rand_num=(rand()*clock())%50;
        a[i]=rand_num;
        b[i]=rand_num;
    }
    //Printing the input
    printf("%s","Input Array\n");
    for(i=0;i<N;i++){
        printf("%d",(int)a[i]);
        printf("%s","\n");
    }
    //serial prefix sum
    for(i=0;i<7;i++){
        printf("%d",(int)b[i]);
        printf("%s","\n");
    }
    for(i=1;i<7;i++){
        b[i]=b[i-1]+b[i];
    }
    //serial output
    printf("%s","Output Array for serial\n");
    for(i=0;i<7;i++){
        printf("%d",(int)b[i]);
        printf("%s","\n");
    }
    //parallel prefix sum
    omp_set_num_threads(P);
    for(int j=0;j<ceil(log2(N));j++){
        for(i=1;i<(N/j+1)/P;i++){
            #pragma omp parallel private(nthreads, tid)
            {
                tid=omp_get_thread_num();
                    if((i*P+tid-pow(2,j))>=0){
                        a[i*P+tid]=a[i*P+tid]+a[int(i*P+tid-pow(2,j))];
                    }
            }
        }
    }
    //printing the output
    printf("%s","Output Array for parallel\n");
    for(i=0;i<N;i++){
        printf("%d",(int)a[i]);
        printf("%s","\n");
    }
return 0;
}