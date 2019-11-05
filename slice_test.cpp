#include <stdio.h>
#include <stdlib.h>
#include <math.h> //for log2 and pow and ceil
#include <time.h> //for performance comparisons
//clang++ -Xpreprocessor -fopenmp parallel_prefix_sum.cpp -o pps -lomp 
//^macOS

int main(int argc,  char *argv[]){
    int N,P;
    float st_macro,pt_macro;
    if(argc>1){
        N=atoi(argv[1]);
    }else{
        N=16; //input size
    }
    if(argc>2){
        P=atoi(argv[2]);
    }else{
        P=N/2; //number of threads
    }
    int local_size=floor(N/P);
    //b for serial, a and c for parallel
    //runs 30 times to get a decent average running time for each input size
    long *a;
    a=new long[N];
    //filling input with random integers
    for(int i=0;i<N;i++){
        long rand_num=(rand()/(clock())%50);
        a[i]=rand_num;
    }
    for(int i = 0; i < P; i++){
        long *array_slice;
        array_slice = new long[N/P];
        for(int j = 0;j<N/P; j++){
            array_slice[j] = a[j+(N/P)*i];
        }
        //send the arrays off
      
    }
    if(N%P!=0){
        
    }
}