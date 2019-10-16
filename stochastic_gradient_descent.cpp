#include <omp.h> //parallel library
#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <time.h> //for performance comparisons
#include <iostream>
using namespace std;
//clang++ -Xpreprocessor -fopenmp stochastic_gradient_descent.cpp -o pps -lomp 
//^macOS


double vec_mul(double *mat_one,double *mat_two){
    double result = 0;
    for(int i = 0; i < 8; i++) result += mat_one[i]*mat_two[i];
    cout << result;
    cout << '\n';
    return result;
}
double sigmoid(double *mat_one,double *mat_two){
    double sigmoid_input=vec_mul(mat_one,mat_two);
    double sigmoid_output = 1/(1+exp(sigmoid_input)*-1);
    return sigmoid_output;
}
//serial sigmoid on 8 inputs
double* fill_array(int N){
    double *array;
    array = new double[N];
    for(int i=0;i<N;i++)array[i]=(rand()/(clock())%5);
    return array;
}
double* serial_descent(int N){
    double *argument,*weights;
    argument = fill_array(N);
    weights = fill_array(N);
    
    while(true){
        //calculate gradient here
        double sigmoid_serial = sigmoid(argument, weights);
        cout << sigmoid_serial;
        cout << '\n';
        double cost = 1/2*pow((1-sigmoid_serial),2);
        if(cost == 0){
            break;
        }
    }
    return weights;
}

int main(int argc,  char *argv[]){
    //serial test
    int N;
    if(argc>1){
        N = atoi(argv[1]);
    }else{
        N = 8;
    }
    clock_t serial_time=clock();
    serial_descent(N);
    serial_time = clock() - serial_time;
    return 0;   
}

