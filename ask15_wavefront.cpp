#include <algorithm> 
#include <cmath>

#include "wavefront.h"
#include "omp.h"

using std::min;
using std::max;

void process_block(double* data, int I, int J, int nx, int ny, int Nx, int Ny){
    int counterx; 
    int countery;
    if((nx+1)%Nx==0){
        counterx=(nx+1)/Nx;
    }
    else{
        //Finding the number of nodes in a selected block
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

void wavefront520(double* data, int nx, int ny, int Nx, int Ny){
    int Nw=Nx+Ny-1;//number of waves
    int J;
    int counter;
    for(int w=0;w<Nw;w++){//looping over the waves
        counter=0;
        #pragma omp parallel for
        for(int I=std::max(0,w-Nx+1);I<=std::min(w,Nx);I++){//wavefront parallelization
            J=min(w-I,Ny-1)-counter;
            //J=w-I;
            process_block(data, I, J, nx, ny, Nx, Ny);
            counter=counter+1;
        }
    }
}

void wavefront420(double* data, int nx, int ny, int Nx){
    
}
