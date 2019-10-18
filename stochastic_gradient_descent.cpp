//Russell Hornbuckle
//2019
//gradient descent on sigmoid in serial, parallel, and CUDA - currently implementing in parallel and testing serial implementaion
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
//serial sigmoid on N inputs
double* fill_array(int N){
    double *array;
    array = new double[N];
    for(int i=0;i<N;i++)array[i]=(rand()/(clock())%5);
    return array;
}
void subtract_array(double *array,double decrease, int size){
    for(int i = 0; i < size; i++){
        array[i] = array[i] - decrease;
    }
    return;
}
void print_weights(double* array, int size){
    for(int i = 0; i < size; i++){
        cout << array[i];
        cout << ' ';
    }
    cout << '\n';
}
double* serial_descent(int N, float learning_step){
    double *argument,*weights;
    argument = fill_array(N);
    weights = fill_array(N);
    while(true){
        double sigmoid_serial = sigmoid(argument, weights);
        double gradient =1.0 - sigmoid_serial*(1.0-sigmoid_serial);   
        double decrease = learning_step*gradient;
        subtract_array(weights,decrease,N);   
        print_weights(weights, N);
        double cost = 0.5*pow((1-sigmoid(argument, weights)),2);
        if(cost == 0){
            break;
        }
        //likely will need another method for breaking the descent
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
    serial_descent(N, 0.01);
    serial_time = clock() - serial_time;
    return 0;   
}

