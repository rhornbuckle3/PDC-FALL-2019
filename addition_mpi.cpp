#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <time.h>
using namespace std;
main(int argc, char *argv[])
{   
    int P, id;
    MPI_Init(&argc, &argv); //starting message passing interface
    MPI_Comm_size(MPI_COMM_WORLD, &P); //getting size
    MPI_Comm_rank(MPI_COMM_WORLD, &id); //getting id/rank
    //head node commands
    printf("Hello");
    if(id == 0){
        long *sum = 0;
        int N;
        //building our array
        if(argc==2){
            N=atoi(argv[1]);
        }else{
            N=16; //input size
        }
        long *a;
        a=new long[N];
        for(int i=0;i<N;i++){
            long rand_num=(rand()/(clock())%50);
            a[i]=rand_num;
        }
        //now we need to split the array and send
        //there is most certainly a way to do this by interacting directly with the memory that is much more efficient -- will look into later
        for(int i = 0; i < P; i++){
            long *array_slice;
            array_slice = new long[N/P];
            for(int j = 0;j<N/P; j++){
                array_slice[j] = a[j+(N/P)*i];
            }
            //send the arrays off
            MPI_Send(&array_slice,int(N/P),MPI_LONG,i,0,MPI_COMM_WORLD);
      
        }
        //handling the remainder
        if(N%P!=0){
            for(int i = 0;i<N%P;i++){
                sum+=a[int(N/P)*P+i];
            }
        }
        //recieve other sums
        for(int i = 0; i < P; i++){
            long *sum_from;
            MPI_Recv(&sum_from,1,MPI_LONG,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            sum += *sum_from;
        }
        long sum_two=0;
        for(int i = 0; i < N; i++){
            sum_two+=a[i];
        }
        printf("%d",sum);
        fflush(stdout);
    }else{
        //worker nodes
        long *sum = 0;
        int N;
        //building our array
        if(argc==2){
            N=atoi(argv[1]);
        }else{
            N=16; //input size
        }
        //recieve array
        long *array = new long[N/P];
        MPI_Recv(&array,int(N/P),MPI_LONG,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        for(int i = 0; i < N/P; i++){
            sum+=array[i];
        }
        MPI_Send(&sum,1,MPI_LONG,0,0,MPI_COMM_WORLD);
    }
    
    MPI_Finalize();
}