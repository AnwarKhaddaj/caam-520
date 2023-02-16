#include <algorithm> 
#include <cmath>
#include <cstdio>

#include "wavefront.h"
#include "omp.h"

using std::min;
using std::max;
using std::printf;

const double PI = M_PI;

int elements(int wave, int I, int J, int waves);//helper method to calculate the number of elements in a specific wave.

void process_block(double* data, int I, int J, int nx, int ny, int Nx, int Ny) {
    int x; 
    int y;
    //determine number of nodes in this specific block in x direction
    if((nx + 1) % Nx == 0) {
        x = (nx + 1) / Nx;
    }
    else 
        x = (nx+1)/Nx+1;
    //determine number of nodes in this specific block in y direction
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
            }
            else if(j == 0) {
                data[cartesian2flat(i, j, ny + 1)] = sin(2 * M_PI * i / nx);
            }
            else {
                data[cartesian2flat(i, j, ny + 1)] = Cx * data[cartesian2flat(i - 1, j, ny + 1)] + Cy * data[cartesian2flat(i, j - 1, ny + 1)];
            }
        }
    }
}
void wavefront520(double* data, int nx, int ny, int Nx, int Ny) {
    int Nw = Nx + Ny - 1;//total number of waves
    int J; 
    int* I_track = new int[Nx*Ny];
    int* J_track = new int[Nx*Ny];
    int count = 0;
    //storing the I indices and J indices while traversing them in a diagonal-wave format
    for(int w=0; w<Nw; w++) {
        for(int I=std::max(0, 1 + w - Nx); I<=std::min(Nx, w); I++) {
            J = w-I;
            I_track[count] = I;
            J_track[count] = J;
            count++;
            }
        }
    }
    int Nt; //number of threads
    #pragma omp parallel 
    {
    Nt = omp_get_num_threads();
    }
    //spin-up phase
    for (int w=0; w<Nt-1; w++) {
        #pragma omp parallel for
        for (int I=0; I<=w; I++){
            J = w-I;
            process_block(data, I, J, nx, ny, Nx, Ny);
        }
    }
    int Nb_fullyparallel = Nx*Ny-Nt*(Nt-1); //number of blocks in fully parallelized region
    int Nc;//number of cycles: number of parallel regions
    int leftover=0;
    if (Nb_fullyparallel % Nt == 0)
        Nc = Nb_fullyparallel/Nt;
    else {
        Nc = Nb_fullyparallel/Nt;
        leftover = Nb_fullyparallel%Nt;
    }
    //fully parallel
    for (int k = 0; k < Nc; k++) {
    #pragma omp parallel for
    for (int i = 0; i < Nt; i++) {
        process_block(data, I_track[Nt*(Nt-1)/2+k*Nt+i],J_track[Nt*(Nt-1)/2+k*Nt+i],nx,ny,Nx,Ny);
        }
    }
    if (leftover != 0) {
        #pragma omp parallel for
        for (int i = 0; i < leftover; i++) {
            process_block(data, I_track[Nt*(Nt-1)/2+(Nc-1)*Nt+(Nt-1)+i+1],J_track[Nt*(Nt-1)/2+(Nc-1)*Nt+(Nt-1)+i+1],nx,ny,Nx,Ny);
        }
    }
    //spin-down phase
    int current_wave;
        int w;
        int wave;
    for ( w = Nt-1-1; w >= 0; w--) {
        wave = Nx+Ny-2-w;
        #pragma omp parallel for
        for (int I = Nx-1-w; I <= Nx-1; I++ ){
            J = wave-I;
            process_block(data, I, J, nx, ny, Nx, Ny);
        }
    }
}

int elements(int wave, int I, int J, int Nw) { //how many elements in a specific wave
    if (wave <= min(I,J)) 
        return wave+1;
    else if (wave > min(I,J) && wave <= Nw - min(I,J))
        return min(I,J);
    else
        return Nw-wave;
}

void wavefront420(double* data, int nx, int ny, int Nx) {
    
}
