#include <omp.h> //parallel library
#include <stdio.h>
#include <stdlib.h>
#include <math.h> //for log2 and pow and ceil
#include <time.h> //for performance comparisons
//clang++ -Xpreprocessor -fopenmp parallel_prefix_sum.cpp -o pps -lomp 
//^macOS

int main(int argc,  char *argv[]){
    int nthreads, tid;
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
    for(int q=0;q<30;q++){
        clock_t serial_time,parallel_time;
        long *a,*b,*c;
        a=new long[N];
        b=new long[N];
        c=new long[N];
        //filling input with random integers
        for(int i=0;i<N;i++){
            long rand_num=(rand()/(clock())%50);
            a[i]=rand_num;
            b[i]=rand_num;
            c[i]=rand_num;
        }
        //Printing the input
        //printf("%s","Input Array\n");
        //for(i=0;i<N;i++)printf("%d : %ld\n",i,a[i]);
        //serial prefix sum
        serial_time=clock();
        omp_set_num_threads(1);
        #pragma omp parallel private(nthreads,tid)
        {
            for(int i=1;i<N;i++) b[i]=b[i-1]+b[i];
        }
        serial_time=clock()-serial_time;
        //parallel prefix sum
        parallel_time=clock();
        omp_set_num_threads(P);
        int n=N-1;
        for(int j=0;j<(int)ceil(log2(N));j++){
            int j_offset=pow(2,j);
            int ops=N-j_offset;
            #pragma omp parallel private(nthreads, tid)
                {  
                tid=omp_get_thread_num();
                int ops_per_run=floor(ops/P);
                for(int i=0;i<ops_per_run;i++){
                //ops here 
                int pos=n-(i*P+tid);
                //printf("C[%d](%ld)+A[%d](%ld)=%ld \n",pos,c[pos],pos-j_offset,a[pos-j_offset],c[pos]+a[pos-j_offset]);
                    c[pos]=c[pos]+a[pos-j_offset];
                }
                if(ops%P>tid){
                    int pos=n-(ops_per_run*P+tid);
                    //printf("C[%d](%ld)+A[%d](%ld)=%ld \n",pos,c[pos],pos-j_offset,a[pos-j_offset],c[pos]+a[pos-j_offset]);
                    c[pos]=c[pos]+a[pos-j_offset];
                    //cycle extra ops here
                }
                //updating reference array in parallel (very slow)
                //#pragma omp barrier
                //for(int i=0;i<local_size;i++) a[i*P+tid]=c[i*P+tid];
                //if(N%P>tid) a[local_size*P+tid]=c[local_size*P+tid];
            }
            //updating reference array in serial
            for(int i=0;i<N;i++) a[i]=c[i];
        }
        parallel_time=clock()-parallel_time;
        if(q==0){
            pt_macro=(float)parallel_time;
            st_macro=(float)serial_time;
        }else{
            pt_macro+=(float)parallel_time;
            st_macro+=(float)serial_time;
        }
        //printing the output
        //printf("%s","Output Array for serial\n");
        //for(i=0;i<N;i++)printf("%d\n",(int)b[i]);
        //printf("%s","Output Array for parallel\n");
        //for(i=0;i<N;i++) printf("%d %ld\n",(int)a[i],a[i]-b[i]);
        //printf("clocks per second: %d\n",CLOCKS_PER_SEC);
        //printf("serial time in seconds: %f\n",(float)(serial_time)/CLOCKS_PER_SEC);
        //printf("parallel time in seconds: %f\n",(float)(parallel_time)/CLOCKS_PER_SEC);
        //printf("easy copy: %f\n",((float)(serial_time)/CLOCKS_PER_SEC)/((float)(parallel_time)/CLOCKS_PER_SEC));
    }
    pt_macro=(float)pt_macro/30;
    st_macro=(float)st_macro/30;
    printf("serial time: %f\n",st_macro);
    printf("parallel time: %f\n",pt_macro);
    printf("easy copy: %f\n",(float)st_macro/(float)pt_macro);
return 0;
}