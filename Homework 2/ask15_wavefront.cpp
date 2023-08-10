#include <algorithm> 
#include <cmath>
#include <cstdio>

#include "wavefront.h"
#include "omp.h"

using std::max;
using std::min;
using std::printf;

void process_block(double* data, int I, int J, int nx, int ny, int Nx, int Ny) {
    int counterx; 
    int countery;

    //number of nodes in a specific block in x dimension
    if((nx+1)%Nx==0) {
        counterx = (nx+1)/Nx;
    }
    else 
        counterx = (nx+1)/Nx+1;
    //number of nodes in a specific block in y dimension
    if((ny+1)%Ny==0) {
        countery = (ny+1)/Ny;
    }
    else 
        countery = (ny+1)/Ny+1;
    //backward finite difference method
    for(int i=I*counterx; i<min((I + 1)*counterx, nx+1); i++) {
        for(int j=J*countery; j<min((J + 1)*countery, ny+1); j++) {
            if(i == 0) {
                data[cartesian2flat(i,j,ny+1)] = sin(2*M_PI*j/ny);
            }
            else if(j == 0) {
                data[cartesian2flat(i,j,ny+1)] = sin(2*M_PI*i/nx);
            }
            else {
                data[cartesian2flat(i,j,ny+1)] = Cx*data[cartesian2flat(i-1,j,ny+1)] + Cy*data[cartesian2flat(i,j-1,ny+1)];
            }
        }
    }
}
void wavefront520(double* data, int nx, int ny, int Nx, int Ny) {
    int Nw = Nx+Ny-1; //nb of waves
    int* I_tracker = new int[Nx * Ny];
    int* J_tracker = new int[Nx * Ny];
    int count = 0;
    int J; 
    //storing the i indices and j indices while traversing them in a diagonal-wave format
    for(int w=0; w<Nw; w++) {
        for(int I=std::max(0,1+w-Nx); I<=std::min(Nx,w); I++) {
            J = w-I;
            if (I<=Nx-1 && J<=Ny-1) {//some cases went over Ny-1/Nx-1
                I_tracker[count] = I;
                J_tracker[count] = J;
                count++;
            }
        }
    } 
    int Nt;//nb of threads
    #pragma omp parallel 
    {
        Nt = omp_get_num_threads();
    }
    //spin-up phase
    for (int w=0; w<Nt-1;w++) {
        #pragma omp parallel for
        for (int I=0; I<=w;I++){
            J = w-I;
            process_block(data, I, J, nx, ny, Nx, Ny);
        }
    }

    int Nb_fullyparallel = Nx*Ny-Nt*(Nt-1); //nb of blocks in fully parallelized region
    int Nc; //nb of cycles: each cycle represents a parallel work of the threads i.e nb of threads
    int remainder = 0; //this will be the last part before spin-down: some final work between fully parallel and spin-down
    if (Nb_fullyparallel%Nt==0)
        Nc = Nb_fullyparallel/Nt;
    else {
        Nc = Nb_fullyparallel/Nt;
        remainder = Nb_fullyparallel%Nt;
    }
    //fFully spun-up phase
    for (int k=0;k<Nc;k++) {
        #pragma omp parallel for
        for (int i=0;i<Nt;i++) {
            process_block(data, I_tracker[Nt * (Nt - 1) / 2 + k * Nt + i], J_tracker[Nt * (Nt - 1) / 2 + k * Nt + i], nx, ny, Nx, Ny);
        }
    } 
    //Remainder part of the fully parallel region
    if (remainder != 0) {
        #pragma omp parallel for
        for (int i=0; i<remainder; i++) {
            process_block(data, I_tracker[Nt * (Nt - 1) / 2 + (Nc - 1) * Nt + (Nt - 1) + i + 1], J_tracker[Nt * (Nt - 1) / 2 + (Nc - 1) * Nt + (Nt - 1) + i + 1], nx, ny, Nx, Ny);
        }
    }
    //spin-down phase
    int w;
    int Wave_nb;
    for (w =Nt-2; w>=0;w--) {
        Wave_nb=Nx+Ny-2-w;
        #pragma omp parallel for
        for (int I=Nx-1-w; I<=Nx-1;I++){
            J=Wave_nb-I;
            process_block(data, I, J, nx, ny, Nx, Ny);
        }
    }
    delete[] I_tracker;
    delete[] J_tracker;
}
void wavefront420(double* data, int nx, int ny, int Nx) {
    
}
