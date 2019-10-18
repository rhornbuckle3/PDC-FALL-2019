//Russell Hornbuckle
//2019
//Decaying sort
/*
N samples, need to sort.
Can only compare 2 at a time.
Values decay after each comparison at a fixed but unknown rate unique to each sample. 
Rate of decay after a comparison is a function of the other sample's value. 
Cannot directly view each sample's value.

Obviously, it's impossible to get a perfect order. How would we go about getting a best fit?

-Assign a placeholder value to each member and recompute as a function of their previous comparisons at each round
-
*/

#include <math.h> 
#include <time.h> //for performance comparisons
#include <iostream>


double* fill_array(int N){
    double *array;
    array = new double[N];
    for(int i = 0; i < N; i++){
        array[i] = (double)(rand()*clock())/7.3;
        array[i] = (rand()*clock()%30)+fmod(array[i],3);
    }
    return array;
}
int main(int argc,  char *argv[]){
    int N;
    if(argc>1){
        N = atoi(argv[1]);
    }else{
        N = 32;
    }
    double *array = fill_array(N);
    for(int i = 0; i<N;i++){
        std::cout << array[i];
        std::cout << '\n';
    }

}