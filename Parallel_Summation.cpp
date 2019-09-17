#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc,  char *argv[]){
    int nthreads, tid;
    int P=4; //number of threads
    int N=1024; //input size
    if(argc>1) P=atoi(argv[1]);
    if(argc>2) N=atoi(argv[2]);
    long* input;
    long* output;

    input=new long[N];
    output=new long(P);

    for(int i=0;i<N;i++) input[i]=randn();
    omp_set_num_threads(P);
    #pragma omp parallel private(nthreads, tid){
        int start=0,stop=0,size=0;
        tid=omp_get_thread_num();
        size=N/P;
        start=tid*size;
        stop=start+size;

        int total=0;
        for(int i=start;i<stop;i++){
            total+=input[i];
        }
        output[tid]=total;
        printf("subtotal = %d, thread = %d",total, tid);
    }
    int t=omp_get_num_threads();
    printf("Number of threads = %d\n",t);
    int total=0;
    for(int i=0;i<P;i++) total+=output[i];
    
    //delete input;
    //delete output;
    #pragma omp parallel private(tid,nthreads){
        tid=omp_get_thread_num();

        //outNext=
        int next=outout[tid];
        for(int i=0;i<lg(p);i++)
        next+=output[tid+2^i];
        #pragma omp barrier
    }
}
