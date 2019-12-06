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
                array[i][j] = (rand()/(clock())%50) + 1/((rand()*clock())%10000);
            }
        }
    }
    return array;
}
main(int argc, char *argv[]){
    //9 processors
    int P, id;
    MPI_Init(&argc, &argv); //starting message passing interface
    MPI_Comm_size(MPI_COMM_WORLD, &P); //getting size
    MPI_Comm_rank(MPI_COMM_WORLD, &id); //getting id/rank
    //printf("Init success, # of nodes: %d\n",P);
    if(id>0 && id <10){
        double **workingArrayA = arrayConstructor(1000,1000,false);
        double **arrayGenerator = arrayConstructor(1000,1000,true);
        
        //printf("Array init success on node: %d\n", id);
        int row = (int)ceil(float(id)/3);
        int col = id%3;  
        if(col == 0) col = 3;
        if(row == 0) row = 1;
        //initial shift A
        int iter = 0;
        for(int i = 0; i < 1000; i++){
            for(int j = 0; j <1000; j++){
                if(j+iter<0){
                    workingArray[i][j] = arrayGenerator[i][1000+iter];
                }else{
                    workingArray[i][j] = arrayGenerator[i][j+iter];
                }
            }
            iter--;
        }
        delete arrayGenerator;
        double **arrayGenerator = arrayConstructor(1000,1000,true);
        double **workingArrayB = arrayConstructor(1000,1000,false);
        //initial shift B
        iter = 0;
        for(int j = 0; j <1000; j++){
            for(int i = 0; i < 1000; i++){
                if(i+iter<0){
                    workingArray[i][j] = arrayGenerator[1000+iter][j];
                }else{
                    workingArray[i][j] = arrayGenerator[i+iter][j];
                }
            }
            iter--;
        }
        delete iter;
        delete arrayGenerator;
        double **resultArray = arrayConstructor(1000,1000,false);
        int colSend,colRec,rowSend,rowRec;
        //printf("Col: %d\n",col);
        //printf("Row: %d\n",row);
        if(col - 1==0){ 
            colSend = 3; 
        }else{ 
            colSend = col - 1;
        }
        if(col+1>3){ 
            colRec = 1; 
        }else{ 
            colRec = col + 1;
        }
        if(row - 1 == 0){ 
            rowSend = 3; 
        }else{ 
            rowSend = row - 1;
        }
        if(row+1>3){ 
            rowRec = 1; 
        }else{ 
            rowRec = row + 1;
        }
        
        for(int k = 0; k<999; k++){
            for(int i = 0; i<1000; i++){
                for(int j = 0; j<1000; j++) resultArray[i][j] = workingArrayA[i][j]*workingArrayB[i][j];
            }
            //printf("%d calc success\n", k);
            //communication
            //row communication (array A)
            //double  *arraySlice = new double[1000];
            double arraySlice[1000];
            for(int q = 0; q<1000; q++) arraySlice[q] = workingArrayA[q][0];
            //printf("arraySlice created successfully\n");
            int send = row*3-3+colSend;
            int recv = row*3-3+colRec;
            //printf("%d\n",colSend);
            //printf("Sending to %d from %d\n",send,id);
            MPI_Send(arraySlice,1000,MPI_DOUBLE,send,0,MPI_COMM_WORLD);
            //printf("MPISEND success\n");
            MPI_Recv(arraySlice,1000,MPI_DOUBLE,recv,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            //array shift (cols)
            //it was at this point that I greatly regretted not using std::vector
            for(int i = 0; i<999; i++){
                for(int j = 0; j<1000; j++) workingArrayA[j][i] = workingArrayA[j][i+1];
            }
            for(int i = 0; i<1000; i++) workingArrayA[i][999] = arraySlice[i];

            //row communication (array b)
            for(int q = 0; q<1000; q++){
                arraySlice[q] = workingArrayB[0][q];
            }
            send = rowSend*3-3+col;
            recv = rowRec*3-3+col;
            //printf("Sending to %d from %d\n",send,id);
            MPI_Send(arraySlice,1000,MPI_DOUBLE,send,0,MPI_COMM_WORLD);
            //printf("MPISEND success\n");
            MPI_Recv(arraySlice,1000,MPI_DOUBLE,recv,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            //printf("%d communication success\n", k);
            //array shift (rows)
            for(int i = 0; i<999; i++){
                for(int j = 0; j<1000; j++) workingArrayB[i][j] = workingArrayB[i+1][j];
            }
            for(int i = 0; i<1000; i++) workingArrayB[999][i] = arraySlice[i];
        }
        printf("%d",id);
        printf("\n");
        for(int i = 0; i < 1000; i++){
            for(int j = 0; j <1000; j++){
                printf("%f",resultArray[i][j]);
                printf(" ");
            }
            printf("\n");
        }
    //write code for sending and recieving arrays to be printed by the head node
    MPI_Finalize();
    }
}