#include <mpi.h>

#include "prob2.h"

void group_exchange() {
    int nR;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nR);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int* x;
    int* data;
    x = &rank; //to get the rank
    //forcing the order from the start
    if (rank == 0) {
        MPI_Send(x, 1, MPI_INT, (rank+1)%nR, 0, MPI_COMM_WORLD);
    } else {
        MPI_Recv(data, 1, MPI_INT, (rank-1)%nR, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Send(x, 1, MPI_INT, (rank+1)%nR, 0, MPI_COMM_WORLD);
    }   
    //The final rank to receive will be rank 0 
    if(rank == 0) {
        MPI_Recv(data, 1, MPI_INT, (rank-1)%nR, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }  
}