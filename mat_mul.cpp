//Russell Hornbuckle
//2019
//matrix multiplication of 2 same sized N by N arrays in parallel
#include <omp.h> //parallel library
#include <stdio.h>
#include <stdlib.h>
#include <math.h> //for log2 and pow and ceil
#include <time.h> //for performance comparisons
//clang++ -Xpreprocessor -fopenmp mat_mul.cpp -o mat_mul -lomp 
//^macOS
//assumes same size matrices
int main(int argc,  char *argv[]){
    int nthreads, tid;
    int N;
    if(argc>1){
        N = atoi(argv[1]);
    }else{
        N = 16;
    }//input size
    int P=N*N; //number of threads
    int *C,*A,*B;
    

    A=new int[N*N];
    B=new int[N*N];
    C=new int[N*N];

    for(int i=0;i<N*N;i++){
        A[i]=(rand()*clock())%50;
        B[i]=(rand()*clock())%50;
        C[i]=(rand()*clock())%50;
        }

    omp_set_num_threads(P);
    #pragma omp parallel
    {
        int sum=0;
        int pid=omp_get_thread_num();
        int ci=pid/N,cj=pid%N;
        for(int i=0;i<N;i++){
            sum+=A[ci*N+i]*B[cj+N*i];
        }
        C[ci*N+cj]=sum;
    }
    for(int i=0;i<N*N;i++) {
        printf("c: %d",C[i]);
        if(i%N==0){
            printf("\n");
        }
    }
}
int alt(int A,int B,int C){
    int P=16;

}