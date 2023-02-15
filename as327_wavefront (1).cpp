#include <algorithm> 
#include <cmath>
#include <cstdio>

#include "wavefront.h"
#include "omp.h"

using std::min;
using std::max;
using std::printf;

int elements(int wave, int I, int J, int waves);

void process_block(double* data, int I, int J, int nx, int ny, int Nx, int Ny) {
    int x; 
    int y;

    //determine number of nodes in this specific block in x direction
    //base case of 1 block in x direction
    if((nx + 1) % Nx == 0) {
        x = (nx + 1) / Nx;
    }
    else {
        if(Nx == 1) {
            x = nx + 1;   
        }
        //more than one block in x direction
        else {
            if(I < Nx - 1) {
                x = (nx + 1)/(Nx - 1); 
            }
            else {
                x = (nx + 1) % (Nx - 1); 
            }  
        }
    }

    //determine number of nodes in this specific block in y direction
    //base case of 1 block in y direction
    if((ny + 1) % Ny == 0) {
        y = (ny + 1) / Ny;
    }
    else {
        if(Ny == 1) {
            y = ny + 1;   
        }
        //more than one block in y direction
        else {
            if(J < Ny - 1) {
                y = (ny + 1)/(Ny - 1); 
            }
            else {
                y = (ny + 1) % (Ny - 1); 
            }  
        }
    }
    //calculate nodes using backward finite difference method
    for(int i = I * x; i < (I + 1) * x; i++) {
        for(int j = J * y; j < (J + 1) * y; j++) {
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
    int waves = Nx + Ny - 1;
    int counter;
    int J; 
    for(int w = 0; w < waves; w++) {
        counter = 0;
        //parallelize using OpenMP
        //#pragma omp parallel for
        //loop over the blocks in a wavefront diagonal fashion
        for(int I = std::max(0, 1 + w - Nx); I <= std::min(Nx, w); I++) {
            //printf("I = %d\n",1+w-Nx);
            for(int J = min(Ny - 1, w - I); J >= min(Ny - 1, w - I)-elements(w,I,J,waves); J--) {
            //J = min(Ny - 1, w - I) - counter;
            //process the domain in blocks
            if (I+J == w && I <= Nx-1 && J <= Ny -1) {
            process_block(data, I, J, nx, ny, Nx, Ny);
            printf("%d %d %d\n",w,I, J);
            }
            //counter++;
            }
        }
        
        printf("\n");

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
