//Russell Hornbuckle
//2019
//logistic regression with serial and parallel descent
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <iostream>
using namespace std;
//clang++ -Xpreprocessor -fopenmp stochastic_gradient_descent.cpp -o sgd -lomp 
//^macOS
//scp ./stochastic_gradient_descent.cpp rhornbuckle3@hpclogin:/home/users/rhornbuckle3/PDC
//module load Compilers/mvapich2_ACoRE
//mpic
//qCDER27
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
double ** array_constructor(int rows , int cols, bool fill){
    double **array;
    array = new double*[rows];
    for(int i = 0; i < rows; i++){
        array[i] = new double[cols];
        if(fill){
            for(int j = 0; j < cols; j++){
                array[i][j] = (double)((rand()*(clock()%100))%10) + 1.0/((double)((rand()*(clock()%100))%10000+1));
            }
        }
    }
    return array;
}
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

main(int argc, char *argv[]){
    MPI_Init(&argc, &argv); //starting message passing interface
    int P, id;
    MPI_Comm_size(MPI_COMM_WORLD, &P); //getting size
    MPI_Comm_rank(MPI_COMM_WORLD, &id); //getting id/rank
    printf("declaration 1 %d\n", id);
    if(id>0 && id < P){
        
        double const LEARNING_RATE = 0.0001;
        int I, S, N;
        I = 1;
        S = 10;
        N = 500;
        //if(argc>1){I = atoi(argv[1]);}else{I = 1;}
        //if(argc>2){S = atoi(argv[2]);}else{S = 10;}
        //if(argc>3){N = atoi(argv[3]);}else{N = 1000;}
        double **X;
        double *Y;
        double time;
        
        //printf("Number of Features = %d\n",S);
        //printf("Number of Samples = %d\n",N);
    
        if(id == 0){
            printf("declaration 1.1 %d\n", id);
            X = array_constructor(N,S,true);
            Y = fill_array(N,2,false);
            printf("declaration 2.1 %d\n", id);
            MPI_Bcast(X,N*S,MPI_DOUBLE,1,MPI_COMM_WORLD);
        }else{
            printf("declaration 2 %d\n",id);
            MPI_Bcast(X,N*S,MPI_DOUBLE,1,MPI_COMM_WORLD);
        }
        printf("declaration 3 %d\n",id);
        MPI_Barrier(MPI_COMM_WORLD);
        if(id==1){
            MPI_Bcast(Y,N,MPI_DOUBLE,1,MPI_COMM_WORLD);
        }else{
            MPI_Bcast(Y,N,MPI_DOUBLE,1,MPI_COMM_WORLD);
        }
        
        double *beta = fill_array(S,2);
        double* gradient = new double[S];
        double cost_old = 0.0;
        for(int i = 0; i < S; i++) gradient[i] = 0.0;
        int iterator = 1;
        double time_seconds = 0.0;
        double *cost;
        *cost = 0.0;
        int *best; 
        *best = 0;
        //time record start
        time = MPI_Wtime();
        while(true){
            
            //gradient calculation
            for(int i = 0; i < N; i++ ){
                double sigmoid_in =  sigmoid(beta,X[i],S) - Y[i];
                double *gradient_in = multiply_in(sigmoid_in,X[i],S);
                gradient = add_arrays(gradient,gradient_in,S,false);            
            }  
            beta = add_arrays(beta,multiply_in(LEARNING_RATE,gradient,S),S,true);
            //cost calculation
            *cost = 0.0;
            for(int i = 0; i < N; i++){
                double sigmoid_in = sigmoid(beta,X[i],S);
                double cost_one = Y[i]*log10(sigmoid_in);
                double cost_two = (1.0 - Y[i])*log10(1-sigmoid_in);
                *cost += cost_one + cost_two;
            }
            
            *cost = (1.0/(double)N)*(*cost)*-1.0;
            if(id != 1 || id != 0){
                MPI_Send(&cost,1,MPI_DOUBLE,1,id,MPI_COMM_WORLD);
                //MPI_Recv();
            }
            if(id==1){
                double *cost_array = new double[P-1];
                for(int i = 2; i < P;i++){
                MPI_Recv(&cost_array[i],1,MPI_DOUBLE,i,i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                }
                double comparator = cost_array[0];
                for(int i = 1; i < P; i++){
                    if(cost_array[i]<comparator){
                        *best = i;
                        comparator = cost_array[i];
                    }    
                }
            }
            double *beta_best;
            MPI_Bcast(&best,1,MPI_INT,1,MPI_COMM_WORLD);
            if(id == *best){
                beta_best = beta;
            }
            MPI_Bcast(&beta_best,S,MPI_DOUBLE,*best,MPI_COMM_WORLD);
            beta = average_weights(beta,beta_best,S);
            //std::cout << "Cost: " << cost << endl;
            if(cost_old == *cost){
                //printf("Descent iteratiion: %d\n",iterator);
                break;
            }
            iterator++;
            cost_old = *cost;
        }
        //sync and report
        if(id == *best){
            time = time - MPI_Wtime();
            printf("Time: %lu\n",time);
            printf("Iterations: %d\n",iterator);
        }
        MPI_Finalize();
    }
}