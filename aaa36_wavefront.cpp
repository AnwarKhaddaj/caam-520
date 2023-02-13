#include <iostream>
#include <cstdio>
#include "wavefront.h"
#include <cmath>
#include <omp.h>



void wavefront420(double* data, int nx, int ny, int Nx);

void process_block(double* data, int I, int J, int nx, int ny, int Nx, int Ny){
    // Calculate the start and end indices for the nodes in the current block
    // with nx and ny being the size of the FD grid
    // Nx and Ny are the number of blocks in each dimension
    int i_start = I*(nx/Nx); 
    int i_end = (I+1)*(nx/Nx) - 1;
    if (I == Nx - 1){
        i_end = nx - 1;
    }

    int j_start = J*(ny/Ny);
    int j_end = (J+1)*(ny/Ny) - 1;
    if (J == Ny - 1){
        j_end = ny - 1;
    }

    //Loop over all the nodes in current block with implementation of cartesian2flat funciton
    // Using ++i/++j to return the value after incrementation
    for (int i = i_start; i <=i_end; ++i){
        
        for (int j = j_start; j <= j_end; ++j){

            int curr = cartesian2flat(i,j,ny);
            int up = cartesian2flat(i - 1, j, ny);
            int down = cartesian2flat(i + 1 ,j ,ny);
            int right = (i, j+1, ny);
            int left = (i, j - 1, ny);

            //updating node that is currently being processed
            data[curr] = data[curr] + Cx*(data[right] + data[left] - 2*data[curr]) + Cy*(data[up] + data[down] - 2*data[curr]);
        }
    


    }
}
        


void wavefront520(double* data, int nx, int ny, int Nx, int Ny){
    // Determining the number of blocks in each direction
    int nb_x = nx/Nx;
    int nb_y = ny/Ny;

    // initializations:
    int I, J, i, j;
    #pragma omp parallel for private(I, J, i, j) schedule (dynamic)
    for (I = 0; I < Nx; I++){
        for (J = 0; J < Ny; J++){
            
            process_block(data, I, J, nx, ny, Nx, Ny);
        }
    }

    //Implement wrap around for when the number of blocks in each dimension does not equal the number of threads
    if (nx % Nx != 0 || ny % Ny != 0){
        process_block(data, Nx,Ny,nx,ny,Nx,Ny);
    }


}
