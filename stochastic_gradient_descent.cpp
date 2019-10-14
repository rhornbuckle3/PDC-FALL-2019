#include <omp.h> //parallel library
#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <time.h> //for performance comparisons
//clang++ -Xpreprocessor -fopenmp stochastic_gradient_descent.cpp -o pps -lomp 
//^macOS

//serial sigmoid on 8 inputs
double mat_mul(double *mat_one,double *mat_two){
    //fill in here

    double result;
    return result;
}

int main(){
    double *argument;
    argument = new double[8];
    double *weights;
    weights = new double[8];
    for(int i=0;i<8;i++)weights[i]=(rand()/(clock())%50);
    for(int i=0;i<8;i++)argument[i]=(rand()/(clock())%50);
    double sigmoid_input=mat_mul(argument,weights);
    double sigmoid = 1/(1+exp(sigmoid_input));
}