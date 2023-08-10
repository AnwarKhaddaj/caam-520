#include <mpi.h>

#include "halo_exchange.h"

// writes from memory to buffer
void pack(int* memory, int* buffer, int nx_local, int ny_local, int halo_radius) {
    int count = 0;
    int jump = 0;
    for (int j = 0; j < ny_local; j++){
        for (int i = 0; i < halo_radius; i++) {
            buffer[count] = memory[jump];
            count++;
            jump++;
        }
        jump = jump + nx_local+halo_radius; //this is to jump to the right memory
    }
}

// writes from buffer to memory
void unpack(int* memory, int* buffer, int nx_local, int ny_local, int halo_radius) {
    int count = 0;
    int jump = 0;
    for (int j = 0; j < ny_local; j++){
        for (int i = 0; i < halo_radius; i++) {
            memory[jump] = buffer[count];
            count++;
            jump++;
        }
        jump = jump + nx_local+halo_radius; //this is to jump to the right memory
    } 
}

void process_block(int* memory, int nx_local, int ny_local, int Nx, int Ny, int halo_radius) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    int nx_total_local = nx_local + 2*halo_radius;
    int ny_total_local = ny_local + 2*halo_radius;

    //Creating requests for all sends and receives
    MPI_Request request1 = MPI_REQUEST_NULL;
    MPI_Request request2 = MPI_REQUEST_NULL;
    MPI_Request request3 = MPI_REQUEST_NULL;
    MPI_Request request4 = MPI_REQUEST_NULL;
    MPI_Request request5 = MPI_REQUEST_NULL;
    MPI_Request request6 = MPI_REQUEST_NULL;
    MPI_Request request7 = MPI_REQUEST_NULL;
    MPI_Request request8 = MPI_REQUEST_NULL;

    //Top halo receiving from bottom 
    if (rank < Nx)
        MPI_Irecv(&memory[cartesian2flat(halo_radius,0,nx_total_local)],nx_total_local*halo_radius-2*halo_radius,MPI_INT, Nx*Ny-(Nx-rank), 0,MPI_COMM_WORLD,&request5);
    else
        MPI_Irecv(&memory[cartesian2flat(halo_radius,0,nx_total_local)],nx_total_local*halo_radius-2*halo_radius,MPI_INT, rank-Nx, 0, MPI_COMM_WORLD, &request5);

    //Bottom halo receiving from top
    if (Nx*Ny-rank <= Nx)
        MPI_Irecv(&memory[cartesian2flat(halo_radius, ny_local+halo_radius,nx_total_local)],nx_total_local*halo_radius-2*halo_radius,MPI_INT, rank % Nx , 0,MPI_COMM_WORLD,&request6);
    else
        MPI_Irecv(&memory[cartesian2flat(halo_radius, ny_local+halo_radius,nx_total_local)],nx_total_local*halo_radius-2*halo_radius,MPI_INT, rank+Nx, 0, MPI_COMM_WORLD, &request6);
    
    //Left halo receiving from right
    int* buffer3 = new int[halo_radius*ny_local];
    if (rank %Nx==0)
        MPI_Irecv(buffer3,ny_local*halo_radius,MPI_INT, rank + Nx - 1, 0,MPI_COMM_WORLD,&request7);
    else
        MPI_Irecv(buffer3,ny_local*halo_radius,MPI_INT, rank-1, 0, MPI_COMM_WORLD, &request7);
    
    //Right halo receiving from leftmost
    int* buffer4 = new int[halo_radius*ny_local];
    if (rank % Nx == Nx-1)
        MPI_Irecv(buffer4,ny_local*halo_radius,MPI_INT, rank - Nx + 1, 0,MPI_COMM_WORLD,&request8);
    else
        MPI_Irecv(buffer4,ny_local*halo_radius,MPI_INT, rank+1, 0, MPI_COMM_WORLD, &request8);

    // Computing top inner halo
    for (int i=halo_radius;i<nx_local+halo_radius;i++){
        for (int j=halo_radius;j<2*halo_radius;j++){
            memory[cartesian2flat(i,j,nx_total_local)] = rank;  
        }
    }
    //Sending the top inner halo
    if (rank < Nx) 
        MPI_Isend(&memory[cartesian2flat(halo_radius,halo_radius,nx_total_local)], nx_total_local*halo_radius-2*halo_radius, MPI_INT, Nx*Ny-(Nx-rank), 0, MPI_COMM_WORLD, &request3);
    else
        MPI_Isend(&memory[cartesian2flat(halo_radius,halo_radius,nx_total_local)], nx_total_local*halo_radius-2*halo_radius, MPI_INT, rank-Nx, 0, MPI_COMM_WORLD, &request3);
    
    //Computing bottom inner halo
    for (int i=halo_radius;i<nx_local+halo_radius;i++){
        for (int j=ny_local;j<ny_local+halo_radius;j++){
            memory[cartesian2flat(i,j,nx_total_local)] = rank;  
        }
    }
    // Sending the bottom inner halo
    if (Nx*Ny -rank <= Nx)
        MPI_Isend(&memory[cartesian2flat(halo_radius,ny_local,nx_total_local)], nx_total_local*halo_radius-2*halo_radius, MPI_INT, rank % Nx, 0, MPI_COMM_WORLD, &request4);
    else
        MPI_Isend(&memory[cartesian2flat(halo_radius,ny_local,nx_total_local)], nx_total_local*halo_radius-2*halo_radius, MPI_INT, rank+Nx, 0, MPI_COMM_WORLD, &request4);
    
    // Computing left inner halo
    for (int i=halo_radius;i<2*halo_radius;i++){
        for (int j=2*halo_radius;j<ny_local;j++){
            memory[cartesian2flat(i,j,nx_total_local)] = rank; 
        }
    }
    //packing left inner halo
    int* buffer1 = new int[halo_radius*ny_local];
    pack(&memory[cartesian2flat(halo_radius,halo_radius,nx_total_local)], buffer1, nx_local, ny_local, halo_radius);

    //Sending left inner halo
    if (rank % Nx == 0) 
        MPI_Isend(buffer1, halo_radius*ny_local, MPI_INT, rank + Nx - 1, 0, MPI_COMM_WORLD, &request1);
    else
        MPI_Isend(buffer1, halo_radius*ny_local, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &request1);
    
    // Computing right inner halo
    for (int i=nx_local;i<nx_local+halo_radius;i++){
        for (int j=2*halo_radius;j<ny_local;j++){
            memory[cartesian2flat(i,j,nx_total_local)] = rank; 
        }
    }

    // Packing right inner halo
    int* buffer2 = new int[halo_radius*ny_local];
    pack(&memory[cartesian2flat(nx_local,halo_radius,nx_total_local)], buffer2, nx_local, ny_local, halo_radius);

    //Sending right inner halo
    if (rank % Nx == Nx-1) 
        MPI_Isend(buffer2, halo_radius*ny_local, MPI_INT, rank - Nx + 1, 0, MPI_COMM_WORLD, &request2);
    else
        MPI_Isend(buffer2, halo_radius*ny_local, MPI_INT, rank + 1, 0, MPI_COMM_WORLD, &request2);
    
    //Computing innermost region (overlapping communication and computation)
    for (int i=2*halo_radius;i<nx_local;i++){
        for (int j=2*halo_radius;j<ny_local;j++){
            memory[cartesian2flat(i,j,nx_total_local)] = rank;
        }
    }

    // Waiting for the sends 
    MPI_Wait(&request1,MPI_STATUS_IGNORE);
    MPI_Wait(&request2,MPI_STATUS_IGNORE);
    MPI_Wait(&request3,MPI_STATUS_IGNORE);
    MPI_Wait(&request4,MPI_STATUS_IGNORE);
    
    // Waiting for the receives
    MPI_Wait(&request5,MPI_STATUS_IGNORE);
    MPI_Wait(&request6,MPI_STATUS_IGNORE);
    MPI_Wait(&request7,MPI_STATUS_IGNORE);
    MPI_Wait(&request8,MPI_STATUS_IGNORE);
    
    // Unpacking the right and left halos after the receive is over
    unpack(&memory[cartesian2flat(0,halo_radius,nx_total_local)], buffer3,nx_local, ny_local, halo_radius);
    unpack(&memory[cartesian2flat(nx_local+halo_radius,halo_radius,nx_total_local)], buffer4,nx_local, ny_local, halo_radius);

    // Freeing up memory
    delete[] buffer1;
    delete[] buffer2;
    delete[] buffer3;
    delete[] buffer4;  
}