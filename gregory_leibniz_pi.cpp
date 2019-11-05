//Russell Hornbuckle
//2019
//pi calculation via the gregory & leibniz series
#include <omp.h> //parallel library
#include <math.h> //for log2 and pow and ceil
#include <iostream>
//clang++ -Xpreprocessor -fopenmp gregory_leibniz_pi.cpp -o gregory_leibniz_pi -lomp 
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
    double pi_sum = 0;
    omp_set_num_threads(P);
    for(int i = 0;i<=floor(N/P);i++){
        #pragma omp parallel private(nthreads,tid)
        {   
            tid=omp_get_thread_num();
            double index_p = i*P+tid;
            double calculation_top = pow(-1.0,index_p);
            double calculation_bottom = (2.0*(index_p)+1.0);
            double calculation = calculation_top/calculation_bottom;
            #pragma omp atomic
            pi_sum+=calculation;
        }
    }
    int remainder = N%P;
    if(remainder > 0){
        omp_set_num_threads(remainder);
        for(int i = 0; i<N%P; i++){
            #pragma omp parallel private(nthreads,tid)
            {
                tid=omp_get_thread_num();
                double index_p = (N-remainder)+tid+1;
                double calculation_top = pow(-1.0,index_p);
                double calculation_bottom = (2.0*(index_p)+1.0);
                double calculation = calculation_top/calculation_bottom;
                #pragma omp atomic
                pi_sum+=calculation;
            }
        }
    }
    pi_sum = pi_sum*4;
    std::cout << pi_sum;
    std::cout << '\n';

}