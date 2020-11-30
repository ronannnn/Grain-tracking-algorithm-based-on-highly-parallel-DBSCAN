#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <math.h>

int main(int argc, char *argv[]){
    
    int rank, size;
    
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    for(int step=0; step<(int)log2(size); step++){
        int senderRank, recverRank;
        if(step==0){
            MPI_Status status;
            int n = size/2;
            for(int i=0; i<n; i++){
                senderRank = 2*i;
                recverRank = 2*i+1;
                if(rank == senderRank){
                    printf("Im %d as sender and me %d as recver in step %d\n", senderRank, recverRank, step);
            }
        }else{
            int i = pow(2, (step+1)) - 1;
            for(int j=0; j<size/(2 * pow(2, step)); j++){
                recverRank = i;
                senderRank = i - (int) pow(2, step);
                i += pow(2, step+1);
                if(rank == senderRank){
                    printf("Im %d as sender and me %d as recver in step %d\n", senderRank, recverRank, step);
                }
            }
        }
    }
    MPI_Finalize();
    return 0;
}