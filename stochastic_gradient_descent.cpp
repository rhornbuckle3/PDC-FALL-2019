#include <omp.h> //parallel library
#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <time.h> //for performance comparisons
//clang++ -Xpreprocessor -fopenmp stochastic_gradient_descent.cpp -o pps -lomp 
//^macOS


double vec_mul(double *mat_one,double *mat_two){
    double result = 0;
    for(int i = 0; i < 8; i++) result += mat_one[i]*mat_two[i];
    return result;
}
double sigmoid(double *mat_one,double *mat_two){
    double sigmoid_input=vec_mul(mat_one,mat_two);
    double sigmoid_output = 1/(1+exp(sigmoid_input)*-1);
    return sigmoid_output;
}
//serial sigmoid on 8 inputs
double* fill_inputs(){
    double *argument,*weights;
    argument = new double[8];
    weights = new double[8];
    for(int i=0;i<8;i++)weights[i]=(rand()/(clock())%50);
    for(int i=0;i<8;i++)argument[i]=(rand()/(clock())%50);
    return argument,weights;
}
double* serial_descent(){
    double *argument,*weights;
    argument,weights = fill_inputs();
    double sigmoid_serial = sigmoid(argument, weights);
    
    return weights;
}

int main(int argc,  char *argv[]){
    
    
    
    return 0;   
}

