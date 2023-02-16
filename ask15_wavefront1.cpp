#include <algorithm> 
#include <cmath>
#include <cstdio>

#include "wavefront.h"
#include "omp.h"

using std::min;
using std::max;
using std::printf;

int elements(int wave, int I, int J, int waves);
void update_node(double* data, int i, int j, int ny);

const double PI = M_PI;

void process_block(double* data, int I, int J, int nx, int ny, int Nx, int Ny) {
    int x; 
    int y;

    //determine number of nodes in this specific block in x direction
    //base case of 1 block in x direction
    if((nx + 1) % Nx == 0) {
        x = (nx + 1) / Nx;
    }
    else 
        x = (nx+1)/Nx+1;
    

    //determine number of nodes in this specific block in y direction
    //base case of 1 block in y direction
    if((ny + 1) % Ny == 0) {
        y = (ny + 1) / Ny;
    }
    else 
        y = (ny+1)/Ny+1;
    
    //calculate nodes using backward finite difference method
    for(int i = I * x; i < min((I + 1) * x,nx+1); i++) {
        for(int j = J * y; j < min((J + 1) * y,ny+1); j++) {
            if(i == 0) {
                data[cartesian2flat(i, j, ny + 1)] = sin(2 * M_PI * j / ny);
                //printf("%f\n", data[cartesian2flat(i, j, ny + 1)]);
            }
            else if(j == 0) {
                data[cartesian2flat(i, j, ny + 1)] = sin(2 * M_PI * i / nx);
                // printf("%f\n", data[cartesian2flat(i, j, ny + 1)]);
            }
            else {
                data[cartesian2flat(i, j, ny + 1)] = Cx * data[cartesian2flat(i - 1, j, ny + 1)] + Cy * data[cartesian2flat(i, j - 1, ny + 1)];
                //printf("%f\n", data[cartesian2flat(i, j, ny + 1)]);
            }
        }
    }
}

/*
void process_block(double* data, int I, int J, int nx, int ny, int Nx, int Ny) {
    // calculate size of blocks in each dimension
    int stride_I;

    // taking "ceiling"
    if((nx + 1) % Nx == 0)
        stride_I = (nx + 1) / Nx; // if num blocks divides evenly, then divide
    else
        stride_I = (nx + 1) / Nx + 1; // if not, then will have extra remainder block

    int stride_J;

    // taking "ceiling"
    if((ny + 1) % Ny == 0)
        stride_J = (ny + 1) / Ny; // if num blocks divides evenly, then divide
    else
        stride_J = (ny + 1) / Ny + 1; // if not, then will have extra remainder block
    
    // compute i and j ranges of block (I, J)
    int i_low = I * stride_I;
    int j_low = J * stride_J;

    // accounting for uneven blocks at ends of array
    int i_high = min(i_low + stride_I, nx + 1);
    int j_high = min(j_low + stride_J, ny + 1);

    // block size
    int nx_local = i_high - i_low;
    int ny_local = j_high - j_low;

    int num_waves = nx_local + ny_local - 1;

    // process block in column-major order
    for(int i = i_low; i < i_high; i++)
        for(int j = j_low; j < j_high; j++)
            update_node(data, i, j, ny);
}

void update_node(double* data, int i, int j, int ny) {
    if(i == 0)
        data[cartesian2flat(i, j, ny + 1)] = sin(2.0 * PI * j / ny);
    else if(j == 0)
        data[cartesian2flat(i, j, ny + 1)] = sin(2.0 * PI * i / nx);
    else
        data[cartesian2flat(i, j, ny + 1)] = Cx * data[cartesian2flat(i-1, j, ny+1)] + Cy * data[cartesian2flat(i, j-1, ny+1)];
}
*/

void wavefront520(double* data, int nx, int ny, int Nx, int Ny) {
    int waves = Nx + Ny - 1;
    int counter;
    int J; 
    int* I_track = new int[Nx*Ny];
    int* J_track = new int[Nx*Ny];
    int count = 0;
    for(int w = 0; w < waves; w++) {
        counter = 0;
        //parallelize using OpenMP
        //#pragma omp parallel for
        //loop over the blocks in a wavefront diagonal fashion
        
        for(int I = std::max(0, 1 + w - Nx); I <= std::min(Nx, w); I++) {
            //printf("I = %d\n",1+w-Nx);
            //for(int J = w-I; J >= min(Ny - 1, w - I)-elements(w,I,J,waves); J--) {
            //for(int J = min(Ny - 1, w - I); J >= min(Ny - 1, w - I)-elements(w,I,J,waves); J--) {
            J = w - I;
            //process the domain in blocks
            if (I+J == w && I <= Nx-1 && J <= Ny -1) {
            //process_block(data, I, J, nx, ny, Nx, Ny);
            //printf("%d %d %d\n",w,I, J);
            I_track[count] = I;
            J_track[count] = J;
            count++;
            }
            //counter++;
            
        }
        
        //printf("\n");

    }
    // for (int i = 0; i < Nx*Ny; i++) {
    //     printf("%d\n",I_track[i]);
    // }
    // for (int i = 0; i < Nx*Ny; i++) {
    //     printf("%d\n",J_track[i]);
    // }
    int Nt;
    #pragma omp parallel 
    {
    Nt = omp_get_num_threads();
    }
    //spinup
    for (int w = 0; w < Nt-1; w++) {
        #pragma omp parallel for
        for (int I = 0; I <= w; I++ ){
            J = w-I;
            process_block(data, I, J, nx, ny, Nx, Ny);
            printf("%d %d\n",I,J);
        }
        //printf("\n");
    }
    // printf("\n");
    // printf("\n");
    // printf("\n");

    int total_parallel_blocks = Nx*Ny-Nt*(Nt-1);
    int num_cycles;
    int leftover=0;
    if (total_parallel_blocks % Nt == 0)
        num_cycles = total_parallel_blocks/Nt;
    else {
        num_cycles = total_parallel_blocks/Nt;
        leftover = total_parallel_blocks%Nt;
    }
    //printf("%d\n",Nt);

    //fully parallel
    for (int k = 0; k < num_cycles; k++) {
    #pragma omp parallel for
    for (int i = 0; i < Nt; i++) {
        process_block(data, I_track[Nt*(Nt-1)/2+k*Nt+i],J_track[Nt*(Nt-1)/2+k*Nt+i],nx,ny,Nx,Ny);
        //printf("%d %d\n",I_track[Nt*(Nt-1)/2+k*Nt+i],J_track[Nt*(Nt-1)/2+k*Nt+i]);
    }
    }
    // printf("\n");
    // printf("\n");
    // printf("\n");

    if (leftover != 0) {
        #pragma omp parallel for
        for (int i = 0; i < leftover; i++) {
            process_block(data, I_track[Nt*(Nt-1)/2+(num_cycles-1)*Nt+(Nt-1)+i+1],J_track[Nt*(Nt-1)/2+(num_cycles-1)*Nt+(Nt-1)+i+1],nx,ny,Nx,Ny);
            //printf("%d %d\n",I_track[Nt*(Nt-1)/2+(num_cycles-1)*Nt+(Nt-1)+i+1],J_track[Nt*(Nt-1)/2+(num_cycles-1)*Nt+(Nt-1)+i+1]);
        }
    }
    //printf("\n");
    //printf("\n");
    //printf("\n");

    //spindown
    int current_wave;
    //for (int w = Nt-1-1; w >= 0; w--) {
        int w;
        int wave;
    for ( w = Nt-1-1; w >= 0; w--) {
        wave = Nx+Ny-2-w;
        #pragma omp parallel for
        //for (int I = Nx-Nt+1; I <= Nx-Nt+1+w; I++ ){
        for (int I = Nx-1-w; I <= Nx-1; I++ ){
            //current_wave = (Nx+Ny-Nt-1);
            //J = current_wave -I;
            J = wave-I;
            process_block(data, I, J, nx, ny, Nx, Ny);
            //printf("% d %d %d\n",w,I,J);
        }
        //printf("\n");
        //w++;
    }
}

int elements(int wave, int I, int J, int waves) {
    if (wave <= min(I,J)) 
        return wave+1;
    else if (wave > min(I,J) && wave <= waves - min(I,J))
        return min(I,J);
    else
        return waves-wave;
}

void wavefront420(double* data, int nx, int ny, int Nx) {
    
}
