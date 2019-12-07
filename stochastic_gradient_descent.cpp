//Russell Hornbuckle
//2019
//logistic regression with serial and parallel descent
#include <omp.h> //parallel library
#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <time.h> //for performance comparisons
#include <iostream>
#include <limits>
using namespace std;
typedef std::numeric_limits< double > dbl;
//clang++ -Xpreprocessor -fopenmp stochastic_gradient_descent.cpp -o sgd -lomp 
//^macOS
//scp ./stochastic_gradient_descent.cpp rhornbuckle3@hpclogin:/home/users/rhornbuckle3/PDC
//idev cpu 10 1 1000 qCDER EDU1002
//idev 
//qCDER27
//cat /proc/cpu/info
double const LEARNING_RATE = 0.0001;
double vec_mul(double *mat_one,double *mat_two,int S){
    double result = 0;
    for(int i = 0; i < S; i++) result += mat_one[i]*mat_two[i];

    return result;
}
double* fill_array(int S, int limiter = 10, bool float_rand = true){
    double *array;
    array = new double[S];
    if(float_rand){
        for(int i=0;i<S;i++)array[i]=(double)((rand()*(clock()%100))%limiter) + 1.0/((double)((rand()*clock())%10000+1));
    }else{
        for(int i=0;i<S;i++)array[i]=(double)((rand()*(clock()%100))%limiter);
    }
    return array;
}
//this function is begging for formal parameters
double ** matrix_constructor(int X , int Y, bool rando, double *weights = new double[0], int limiter = 10){
    double **array;
    array = new double*[X];
    for(int i = 1; i < X; i++){
        array[i] = new double[Y];
        if(rando){
        for(int j = 0; j < Y; j++) array[i][j] = (double)((rand()*(clock()%100))%limiter) + 1.0/((double)((rand()*(clock()%100))%10000+1));
        }else{
            for(int j = 0; j < Y; j++) array[i][j] = 0.0;
        }
    }
    array[0] = new double[Y];
    for(int i = 0; i < Y; i++) array[0][i] = weights[i];
    return array;
}
double* subtract_array(double *array,double decrease, int size){
    for(int i = 0; i < size; i++) array[i] = array[i] - decrease;
    return array;
}
double* add_arrays(double *array,double *increase, int size, bool subtract){
    if(subtract){
        for(int i = 0; i < size; i++) array[i] = array[i] - increase[i];
    }else{
        for(int i = 0; i < size; i++) array[i] = array[i] + increase[i];
    }
    return array;
}
double* average_weights(double *array_one, double *array_two, int S){
    for(int i =0; i < S; i++){
        array_one[i] = (array_one[i]+array_two[i])/2;
    }
    return array_one;
}
void print_weights(double* array, int size){
    for(int i = 0; i < size; i++){
        //printf("%lf \n",array[i]);
        std::cout << array[i] << "\n";
    }
}
double* multiply_in(double gradient_scalar,double *array, int S){
    for(int i = 0; i < S; i++) array[i] = array[i] * gradient_scalar;
    return array;
} 
double sigmoid(double *mat_one,double *mat_two, int S){
    double sigmoid_input=vec_mul(mat_one,mat_two,S)*-1;
    double sigmoid_output = 1/(1+exp(sigmoid_input));
    return sigmoid_output;
}

double* logistic_regression(int N, int S,double **X, double *Y, double *beta){
    int nthreads,tid;
    double* gradient = new double[S];
    double cost_old = 0.0;
    double* return_values = new double[3];
    for(int i = 0; i < S; i++) gradient[i] = 0.0;
    int iterator = 1;
    double time_seconds = 0.0;
    double cost = 0.0;
    omp_set_num_threads(1);
    #pragma omp parallel private(nthreads,tid)
    {
        printf("Starting descent\n");
        clock_t serial_time=clock();
        while(true){
            //gradient calculation
            for(int i = 0; i < N; i++ ){
                double sigmoid_in =  sigmoid(beta,X[i],S) - Y[i];
                double *gradient_in = multiply_in(sigmoid_in,X[i],S);
                gradient = add_arrays(gradient,gradient_in,S,false);            
            }  
            beta = add_arrays(beta,multiply_in(LEARNING_RATE,gradient,S),S,true);
            //cost calculation
            cost = 0.0;
            for(int i = 0; i < N; i++){
                double sigmoid_in = sigmoid(beta,X[i],S);
                double cost_one = Y[i]*log10(sigmoid_in);
                double cost_two = (1.0 - Y[i])*log10(1-sigmoid_in);
                cost += cost_one + cost_two;
            }
            cost = (1.0/(double)N)*cost*-1.0;
            //std::cout << "Cost: " << cost << endl;
            if(cost_old == cost){
                printf("Descent iteratiion: %d\n",iterator);
                serial_time = clock() - serial_time;
                return_values[0] = (double)serial_time/CLOCKS_PER_SEC;
                printf("Time: %lf\n",return_values[0]);
                break;
            }
            iterator++;
            cost_old = cost;
        }
        return_values[1] = cost;
        return_values[2] = iterator;
    }
    return return_values;
}
//Parallel Descent Logistic Regression
//will initialize random weights for each thread, perform descent once and then compare new values 
//whichever node has new lowest function value will then serve as new weights prime and every thread will average their weights with
//weights prime, process will repeat until minima is found
double* parallel_regression(int N, int S, int P, double **X, double *Y, double **beta_grand){
    int nthreads, tid;
    double *cost_array = new double[P];
    int best;
    //time and cost
    double* return_values = new double[3];
    double cost_old = 0.0;
    bool cont = true;
    int iterator = 1;
    double **gradient = matrix_constructor(P,S,false);
    double **gradient_in = matrix_constructor(P,S,false);
    double *sigmoid_in = new double[P];
    omp_set_num_threads(P);
    clock_t parallel_time;
    #pragma omp parallel private(nthreads,tid)
    {   
        tid=omp_get_thread_num();
        #pragma omp master
        {
            printf("Starting Descent\n");
            parallel_time = clock();
        }
        while(cont){
            for(int i = 0; i < S; i++) gradient[tid][i] = 0.0;
            //printf("WE MADE IT HERE %d\n", tid);
            //descent begin
            for(int i = 0; i < N; i++ ){
                sigmoid_in[tid] =  sigmoid(beta_grand[tid],X[i],S) - Y[i];
                gradient_in[tid] = multiply_in(sigmoid_in[tid],X[i],S);
                gradient[tid] = add_arrays(gradient[tid],gradient_in[tid],S,false);            
            }  
            beta_grand[tid] = add_arrays(beta_grand[tid],multiply_in(LEARNING_RATE,gradient[tid],S),S,true);
            //cost calculation
            cost_array[tid] = 0.0;
            for(int i = 0; i < N; i++){
                sigmoid_in[tid] = sigmoid(beta_grand[tid],X[i],S);
                cost_array[tid] += Y[i]*log10(sigmoid_in[tid]) + (1.0 - Y[i])*log10(1-sigmoid_in[tid]);
            }
            cost_array[tid] = (1.0/(double)N)*cost_array[tid]*-1.0;
            //descent end
            //sync section
            #pragma omp barrier
            #pragma omp master 
            {
                double comparator = cost_array[0];
                best = 0;
                for(int i = 1; i < P; i++){
                    if(cost_array[i]<comparator){
                        best = i;
                        comparator = cost_array[i];
                    }    
                }
                //std::cout << "Best Cost: " << cost_array[best] << endl;
                if(cost_old == cost_array[best]){
                    printf("Descent iteration: %d\n",iterator);
                    cont = false;
                    parallel_time = clock() - parallel_time;
                    return_values[0] = (double)parallel_time/CLOCKS_PER_SEC;
                    printf("Time: %lf\n",return_values[0]);
                }else{
                    iterator++;
                    cost_old = cost_array[best];
                }
            }
            //average weights
            #pragma omp barrier
            if(cont){
                beta_grand[tid] = average_weights(beta_grand[tid],beta_grand[best],S);
            }else{
                break;    
            }
        }   
    }
    return_values[1] = cost_array[best];
    return_values[2] = iterator;
    return return_values;
}
double* nosync_parallel_regression(int N, int S, int P, double **X, double *Y, double **beta_grand){
    int nthreads, tid;
    double *cost_array = new double[P];
    //time and cost
    double *return_values = new double[3];
    double *cost_old  = new double[P];
    for(int i = 0; i < P; i++) cost_old[i] = 0.0;
    int iterator = 1;
    double **gradient = matrix_constructor(P,S,false);
    double **gradient_in = matrix_constructor(P,S,false);
    double *sigmoid_in = new double[P];
    bool *cont = new bool[P];
    double *time = new double[P];
    for(int i =0; i < P;i++) cont[i] = true;
    omp_set_num_threads(P);
    clock_t parallel_time;
    printf("Starting Descent\n");
    #pragma omp parallel private(nthreads,tid)
    {   
        tid=omp_get_thread_num();
        time[tid] = clock();
        while(cont[tid]){
            for(int i = 0; i < S; i++) gradient[tid][i] = 0.0;
            //printf("WE MADE IT HERE %d\n", tid);
            //descent begin
            for(int i = 0; i < N; i++ ){
                sigmoid_in[tid] =  sigmoid(beta_grand[tid],X[i],S) - Y[i];
                gradient_in[tid] = multiply_in(sigmoid_in[tid],X[i],S);
                gradient[tid] = add_arrays(gradient[tid],gradient_in[tid],S,false);            
            }  
            beta_grand[tid] = add_arrays(beta_grand[tid],multiply_in(LEARNING_RATE,gradient[tid],S),S,true);
            //cost calculation
            cost_array[tid] = 0.0;
            for(int i = 0; i < N; i++){
                sigmoid_in[tid] = sigmoid(beta_grand[tid],X[i],S);
                cost_array[tid] += Y[i]*log10(sigmoid_in[tid]) + (1.0 - Y[i])*log10(1-sigmoid_in[tid]);
            }
            cost_array[tid] = (1.0/(double)N)*cost_array[tid]*-1.0;
            //descent end
            //sync section
                //std::cout << "Best Cost: " << cost_array[best] << endl;
            if(cost_old[tid] == cost_array[tid]){
                printf("Descent iteration: %d\n",iterator);
                cont[tid] = false;
                time[tid] = clock() - time[tid];
                printf("Time: %lf\n",time[tid]);
            }else{
                iterator++;
                cost_old[tid] = cost_array[tid];
            }
            //average weights
        }
        #pragma omp barrier   
    }
    int best = 0;
    for(int i = 0;i<P;i++){
        if(time[best]<time[i]){
            best = i;
        }
    }
    return_values[0] = (double)time[best]/CLOCKS_PER_SEC;
    return_values[1] = cost_array[best];
    return_values[2] = iterator;
    return return_values;
}
int main(int argc,  char *argv[]){
    //serial test
    //S = Features, N = Samples, P = Threads
    int I, S, N, P;
    if(argc>1){I = atoi(argv[1]);}else{I = 1;}
    if(argc>2){S = atoi(argv[2]);}else{S = 10;}
    if(argc>3){N = atoi(argv[3]);}else{N = 1000;}
    if(argc>4){P = atoi(argv[4]);}else{ P = omp_get_max_threads();}
    printf("Number of Features = %d\n",S);
    printf("Number of Samples = %d\n",N);
    printf("Number of Threads = %d\n",P);
    double **parallel, **serial, **nosync_parallel;
    serial = new double*[I];
    parallel = new double*[I];
    nosync_parallel = new double*[I];
    for(int i = 0; i < I; i++){
    serial[i] = new double[3];
    parallel[i] = new double[3];
    nosync_parallel[i] = new double[3];
    }
    for(int i = 0; i < I; i++){
        double **X = matrix_constructor(N,S,true,fill_array(S));
        double *Y = fill_array(N,2,false);
        double *beta = fill_array(S,2);
        serial[i] = logistic_regression(N,S,X,Y,beta);
    }
    for(int i = 0; i < I; i++){
        double **X = matrix_constructor(N,S,true,fill_array(S));
        double *Y = fill_array(N,2,false);
        double **beta_grand = matrix_constructor(P,S,true,fill_array(S,2),2);
        parallel[i] = parallel_regression(N,S,P,X,Y,beta_grand);
    }
    for(int i = 0; i < I; i++){
        double **X = matrix_constructor(N,S,true,fill_array(S));
        double *Y = fill_array(N,2,false);
        double **beta_grand = matrix_constructor(P,S,true,fill_array(S,2),2);
        nosync_parallel[i] = nosync_parallel_regression(N,S,P,X,Y,beta_grand);
    }
    for(int i = 0; i < I; i++){
        std::cout << "Serial time (in seconds) for run " << i+1 <<": " << serial[i][0] << " and % error: " << serial[i][1] <<" which took "<< serial[i][2]<< " iterations." << endl;
        std::cout << "Parallel time (in seconds) for run " << i+1 <<": " << parallel[i][0] << " and % error: " << parallel[i][1] <<" which took "<< parallel[i][2]<< " iterations." << endl;
        std::cout << "No synchronization Parallel time (in seconds) for run " << i+1 <<": " << nosync_parallel[i][0] << " and % error: " << nosync_parallel[i][1] <<" which took "<< nosync_parallel[i][2]<< " iterations." << endl;
    }
    return 0;
}