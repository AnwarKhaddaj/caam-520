#include <algorithm>
#include <cmath>
#include <cstdio>
#include <math.h>

#include "omp.h"

#include "wavefront.h"

using std::min;
using std::max;
using std::printf;

const double PI = M_PI;

void update_node(double* data, int i, int j, int ny) {
    data[cartesian2flat(i, j, ny)] = Cx * data[cartesian2flat(i-1, j, ny)] + Cy * data[cartesian2flat(i, j-1, ny)];
}

int cartesian2diagonal(int& Id, int& Jd, int I, int J, int nx, int ny) {
    int diagonal = I + J;

    return diagonal * (diagonal + 1) / 2 + I;
}

void initialize_boundary(double* data, int nx, int ny) {
    #pragma omp parallel
    {
        for(int i=0; i < ny + 1; i++) {
            data[cartesian2flat(i, 0, ny)] = sin(2.0 * PI * i / ny);
        }

        for(int j=0; j < nx + 1; j++) {
            data[cartesian2flat(0, j, ny)] = sin(2.0 * PI * j / nx);
        }
    }
}


void process_block(double* data, int I, int J, int nx, int ny, int Nx, int Ny){
    int counterx;
    int countery;
    //Finding the number of nodes in a selected block
    if((nx+1)%Nx==0){
        counterx=(nx+1)/Nx;
    }
    else{
        if(Nx==1){ //base case
        counterx=nx+1;  
     }
    else{
        if(I<Nx-1){
            counterx=(nx+1)/(Nx-1);
        }
        else{
            counterx=(nx+1)%(Nx-1);
        }  
    }

    }
    
    if((ny+1)%Ny==0){
        countery=(ny+1)/Ny;
    }
    else{
        if(Ny==1){ //base case
        countery=ny+1;  
    }
    else{
        if(J<Ny-1){
            countery=(ny+1)/(Ny-1);
        }
        else{
            countery=(ny+1)%(Ny-1);
        }  
    } 
        
    }
       
    for(int i=I*counterx;i<I*counterx+counterx;i++){
        for(int j=J*countery;j<J*countery+countery;j++){
            if(i > 0 && j > 0)
                printf("Processing node (%d, %d)...\n", i, j);
            if(i==0){
                data[cartesian2flat(i,j,ny+1)]=sin(2*M_PI*j/ny);//initial boundary conditions
            }
            else if(j==0){
                data[cartesian2flat(i,j,ny+1)]=sin(2*M_PI*i/nx);//initial boundary conditions
            }
            else{
                data[cartesian2flat(i,j,ny+1)]=Cx*data[cartesian2flat(i-1,j,ny+1)]+Cy*data[cartesian2flat(i,j-1,ny+1)];//3-point stencil
            }
        }
    }
}



// void process_block(double* data, int I, int J, int nx, int ny, int Nx, int Ny) {
//     int i_low = I * Ny;
//     int j_low = J * Nx;

//     int i_high = min(i_low + Ny, ny+1);
//     int j_high = min(j_low + Nx, nx+1);

//     int nx_local = j_high - j_low;
//     int ny_local = i_high - i_low;

//     int num_diags = ny_local + nx_local - 1;

//     for(int d = 0; d < num_diags; d++){    
//         int j = j_low + min(d, nx_local);
//         int i = i_low + max(0, d - nx_local);

//         while(i <= i_high && j >= j_low) {
//             printf("Updating entry (%d, %d)...\n", i, j);
//             if(i > 0 && j > 0)
//                 update_node(data, i, j, ny);
//             i++;
//             j--;
//         }
//     }
// }

void wavefront420(double* data, int nx, int ny, int Nx) { // Ny will = the number of threads
    return;
}

void wavefront520(double* data, int nx, int ny, int Nx, int Ny){
    int Nw = Nx + Ny - 1;

    int Nt = omp_get_num_threads();

    int spin_up_down_num_blocks = Nt * (Nt - 1) / 2;
    int spin_up_num_waves = Nt - 1;
    int end_parallel_region = Nx * Ny - spin_up_down_num_blocks;

    // spin-up phase
    for(int w = 0; w < spin_up_num_waves; w++){
        #pragma omp parallel for
        for(int J = 0; J <= w; J++) {
            int I = w - J;

            process_block(data, J, I, nx, ny, Nx, Ny);
        }
    }

    // fully parallelized region
    int wave_num = spin_up_num_waves;
    int wave_size = spin_up_num_waves;
    int wave_track = 0;
    int counter_I = 0;
    int counter_J = min(wave_num - counter_I, Ny - counter_I);
    int N = spin_up_down_num_blocks;

    for(; N < end_parallel_region; N += Nt) {
        int chunk_size = Nt; // assuming num_threads <= num blocks per wave

        // compute I's and J's before looping
        int* I_chunk = new int[chunk_size];
        int* J_chunk = new int[chunk_size];

        for(int b = 0; b < chunk_size; b++) {
            I_chunk[b] = counter_I;
            J_chunk[b] = counter_J;

            counter_I++;
            counter_J--;

            if(counter_J < 0 || counter_I > min(wave_num, Nx)) {
                wave_num++;
                counter_J = min(wave_num - counter_I, Ny - counter_I);
                counter_I = 0;
            }
        }

        #pragma omp parallel for
        for(int b = 0; b < chunk_size; b++){
            process_block(data, I_chunk[b], J_chunk[b], nx, ny, Nx, Ny);
        }
    }

    // Leftovers from parallel region
    int pre_spin_down_block_num = end_parallel_region - (N - Nt);

    int* I_chunk_2 = new int[pre_spin_down_block_num];
    int* J_chunk_2 = new int[pre_spin_down_block_num];

    // determine I's and J's
    for(int b = 0; b < pre_spin_down_block_num; b++) {
        I_chunk_2[b] = counter_I;
        J_chunk_2[b] = counter_J;

        counter_I++;
        counter_J--;

        if(counter_J < 0 || counter_I > min(wave_num, Nx)) {
            wave_num++;
            counter_J = min(wave_num - counter_I, Ny - counter_I);
            counter_I = 0;
        }
    }

    #pragma omp parallel for
    for(int b = 0; b < pre_spin_down_block_num; b++) {
        process_block(data, I_chunk_2[b], J_chunk_2[b], nx, ny, Nx, Ny);
    }

    // spin-down phase
    for(int w = spin_up_num_waves; w >= 0; w--) {
        #pragma omp parallel for
        for(int J = Nx; J >= Nx - w; J--) {
            int current_wave = (Nx + Ny - 1) - w;
            int I = current_wave - J;

            process_block(data, I, J, nx, ny, Nx, Ny);
        }
    }




    // for(int w = spin_up_num_waves; w >= 0; w--){
    //     #pragma omp parallel for
    //     for(int I = max(0,w-Nx+1); I >=min(w,Nx); I++) {
    //         int J;

    //         process_block(data, J, I, nx, ny);
    //     }
    // }
}