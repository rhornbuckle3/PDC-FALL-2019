//Russell Hornbuckle
//2019
//cannon's algorithm on 9 nodes of 2 3k by 3k arrays
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <time.h>
using namespace std;

double ** arrayConstructor(int rows , int cols, bool fill){
    double **array;
    array = new double*[rows];
    for(int i = 0; i < rows; i++){
        array[i] = new double[cols];
        if(fill){
            for(int j = 0; j < cols; j++){
                array[i][j] = (rand()/(clock())%50);
            }
        }
    }
    return array;
}
main(){
    //9 processors
    int P, id;
    MPI_Init(&argc, &argv); //starting message passing interface
    MPI_Comm_size(MPI_COMM_WORLD, &P); //getting size
    MPI_Comm_rank(MPI_COMM_WORLD, &id); //getting id/rank

    //head node => generate and distribute starting arrays to workers
    if(id == 0){
            //creating our arrays
            double **arrayA = arrayConstructor(3000,3000,true);
            double **arrayB = arrayConstructor(3000,3000,true);
            double **workingArrayA;
            double **workingArrayB;
            double **resultArray = arrayConstructor(1000,1000,false);
            //need to send skewed versions to respective processors
            //each processor calculates results for a 1000x1000 array
            //send A slice
            for(int i = 0; i < (P**0.5)-1; i++){
                for(int j = 0; j < (P**0.5)-1; j++){
                    double **arraySlice;
                    arraySlice = arrayConstructor(1000,1000,false);
                    //fill with skew
                    for(int k = 0; k < 1000; k++){
                        for(int l = 0; l < 1000; l++){
                            if(((j*1000)+l-k)<0){
                                //overflow
                                arraySlice[k][l] = arrayA[(i*1000)+k][l-k+3000];    
                            }else{
                            arraySlice[k][l] = arrayA[(i*1000)+k][(j*1000)+l-k];
                            }
                        }
                    }
                    //send the arrays off
                    if(i==0){
                        workingArrayA = arraySlice;
                    }else{
                    MPI_Send(&(arraySlice[0][0]),1000*1000,MPI_DOUBLE,i,0,MPI_COMM_WORLD);
                    }
                }
            }
            delete arrayA;
            //send B slice
            for(int i = 0; i < (P**0.5)-1; i++){
                for(int j = 0; j < (P**0.5)-1; j++){
                    double **arraySlice;
                    arraySlice = arrayConstructor(1000,1000,false);
                    //fill with skew
                    for(int k = 0; k <1000; k++){
                        for(int l = 0; l <1000; l++){
                            if(((j*1000)+k-l)<0){
                                //overflow
                                arraySlice[k][l] = arrayA[k-l+3000][(j*1000)+l];    
                            } else{ arraySlice[k][l] = arrayA[(i*1000)+k-l][(j*1000)+l];}
                        }
                    }
                    //send the arrays off
                    if(i==0){ workingArrayB = arraySlice;} else{
                    MPI_Send(&(arraySlice[0][0]),1000*1000,MPI_DOUBLE,i,0,MPI_COMM_WORLD);
                    }
                }
            }
            delete arrayB;
        }
        if(id!=0){
            //calculate
            double **workingArrayA;
            double **workingArrayB;
            double **resultArray = arrayConstructor(1000,1000,false);
            MPI_Recv(&workingArrayA,1000*1000,MPI_DOUBLE,1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            MPI_Recv(&workingArrayB,1000*1000,MPI_DOUBLE,1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
            int row = (id+1)/3-1;
            int col = id%3;  
            int colSend,colRec,rowSend,rowRec;
            if(col - 1 <0){ colSend = 2;} else{ colSend = col - 1;}
            if(col+1>2){ colRec = 0;} else{ colRec = col + 1;}
            if(row - 1 < 0){ rowSend = 2;} else{ rowSend = row - 1;}
            if(row+1>2){ rowRec = 0;} else{ rowRec = row + 1;}
    
            for(int k = 0; k<1000; k++){
                for(int i = 0; i<1000; i++){
                    for(int j = 0; j<1000; j++) resultArray[i][j] = workingArrayA[i][j]*workingArrayB[i][j];
                }

                //communication
                //row communication (array A)
                double *arraySlice =new double[1000];
                for(int q = 0; q<1000; q++) arraySlice[q] = workingArrayA[q][0];
                MPI_Send(&arraySlice,1000,MPI_DOUBLE,row*3+colSend,0,MPI_COMM_WORLD);
                MPI_Recv(&arraySlice,1000,MPI_DOUBLE,row*3+colRec,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

                //array shift (cols)
                for(int i = 0; i<1000; i++){
                    for(int j = 0; j<1000; j++) workingArrayA[j][i] = workingArrayA[j][i+1];
                }
                for(int i = 0; i<1000; i++) workingArrayA[i][999] = arraySlice[i];

                //row communication (array b)
                for(int q = 0; q<1000; q++){
                    arraySlice[q] = workingArrayB[0][q]];
                }
                MPI_Send(&arraySlice,1000,MPI_DOUBLE,rowSend*3+col,0,MPI_COMM_WORLD);
                MPI_Recv(&arraySlice,1000,MPI_DOUBLE,rowRec*3+col,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

                //array shift (rows)
                for(int i = 0; i<1000; i++){
                    for(int j = 0; j<1000; j++) workingArrayB[i][j] = workingArrayB[i+1][j];
                }
                for(int i = 0; i<1000; i++) workingArrayB[999][i] = arraySlice[i];
            }
        }
    }