#include <omp.h> //parallel library
#include <stdio.h>
#include <stdlib.h>
#include <math.h> //for log2 and pow and ceil
#include <time.h> //for performance comparisons
#include <iostream>
using namespace std;
//clang++ -Xpreprocessor -fopenmp sieve_eratosthenes.cpp -o sieve_e -lomp 
//^macOS

int main(int argc, char* argv[]){
    int nthreads,tid;
    int N;
    int P;
    if(argc==3){
        P = atoi(argv[1]);
        N = atoi(argv[2]);
    }else{
        return 0;
    }
    bool *list_primes = new bool[N];
    omp_set_num_threads(P);
    //sets the entire array to true
    memset(list_primes,true,sizeof(list_primes));
    for(int i = 2; i*i<N; i++){
        if(list_primes[i]){
            for(int j = i*i; j < (N-i*i)/P;j+=i*P){
                //gonna probably swap over to omp for to make this easier on me
                #pragma omp parallel private(nthreads,tid)
                {
                    tid = omp_get_thread_num();
                    list_primes[j+tid] = false;
                }
            }
            omp_set_num_threads((N-i*i)%P);
            #pragma omp parallel private(nthreads,tid)
            {
                tid = omp_get_thread_num();

            }
        }
    }

