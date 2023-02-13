#include "wavefront.h"
#include "omp.h"
#include <algorithm> //ask about this 
#define _USE_MATH_DEFINES
using std::min;
using std::max;
#include <cstdio>
#include <cmath>
#include <math.h>

void process_block(double* data, int I, int J, int nx, int ny, int Nx, int Ny){
    int countery=0;
    int counterx=0;
    if(I<Nx-1){
        counterx=(nx+1)/(Nx-1); //1280
    }
    else{
        counterx=(nx+1)%(Nx-1); //1
    }
    if(J<Ny-1){
        countery=(ny+1)/(Ny-1); //1280
    }
    else{
        countery=(ny+1)%(Ny-1); //1
    }
    for(int i=I*counterx;i<I*counterx+counterx;i++){
        for(int j=J*countery;j<J*countery+countery;j++){
           if(i==0){
            data[cartesian2flat(i,j,ny)]=sin(2*M_PI*j/ny);
            printf("2");
           }
           else if(j==0){
            data[cartesian2flat(i,j,ny)]=sin(2*M_PI*i/nx);
            printf("4");
           }
           else{
             data[cartesian2flat(i,j,ny)]=Cx*data[cartesian2flat(i-1,j,ny)]+Cy*data[cartesian2flat(i,j-1,ny)];
             printf("3");
           }
        }
    }
}

void wavefront520(double* data, int nx, int ny, int Nx, int Ny){
    int Nw=Nx+Ny-1;//10+0=11 waves
    int J; int counter;
    for(int w=0;w<Nw;w++){
        counter=0;
        #pragma omp parallel for
        for(int I=std::max(0,w-Nx+1);I<std::min(w,Nx);I++){
            J=min(w-I,Ny-1)-counter;
            printf("inside for parallel 1");
            process_block(data, I, J,nx, ny, Nx, Ny);
            counter=counter+1;
        }
    }
}

// void wavefront520(double* data, int nx, int ny, int Nx, int Ny){
//     int Nw=Nx+Ny-1; //number of waves
//     int nb_of_threads=omp_get_num_threads();

//     int J;
//     //Spin-up phase
//     for (int w=0;w<nb_of_threads-1;w++){
//         #pragma omp parallel for
//         {
//         for(int I=0;I<=w;I++){
//             int tid=omp_get_thread_num();
//             J=w-I;
//             #pragma omp critical
//             { //only one thread to execute it at a time
//                 process_block(data,I,J,nx,ny,Nx,Ny);
//             }   
//         }     
//         }
//     }
//     //Fully spun-up phase
//     int fully_spun_up_blocks=Nx*Ny-nb_of_threads*(nb_of_threads-1);//Nx*Ny-spin_up_blocks-spin_down_blocks. Spin_up_blocks=(nbofthreads)*(nbofthreads-1)/2 (Sumtorial formula)
//     int nb_of_thread_cycles=fully_spun_up_blocks/nb_of_threads;
//     int remainder=fully_spun_up_blocks%nb_of_threads;
// }
