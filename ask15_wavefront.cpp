#include "wavefront.h"
#include "omp.h"
// #include <algorithm> //ask about this 
#include <math.h>
//#define _USE_MATH_DEFINES
using std::min;
using std::max;
#include <cmath>


void process_block(double* data, int I, int J, int nx, int ny, int Nx, int Ny){
    int counterx; int countery;
    if(Nx==1){
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
    
    if(Ny==1){
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
    
    for(int i=I*counterx;i<I*counterx+counterx;i++){
        for(int j=J*countery;j<J*countery+countery;j++){
           if(i==0){
            data[cartesian2flat(i,j,ny+1)]=sin(2*M_PI*j/ny);
           }
           else if(j==0){
            data[cartesian2flat(i,j,ny+1)]=sin(2*M_PI*i/nx);
           }
           else{
             data[cartesian2flat(i,j,ny+1)]=Cx*data[cartesian2flat(i-1,j,ny+1)]+Cy*data[cartesian2flat(i,j-1,ny+1)];
           }
        }
    }
}

void wavefront520(double* data, int nx, int ny, int Nx, int Ny){
    int Nw=Nx+Ny-1;
    int J; int counter;
    for(int w=0;w<Nw;w++){
        counter=0;
        #pragma omp parallel for
        for(int I=std::max(0,w-Nx+1);I<=std::min(w,Nx);I++){
            J=min(w-I,Ny-1)-counter;
            process_block(data, I, J,nx, ny, Nx, Ny);
            counter=counter+1;
        }
    }
}


void wavefront420(double* data, int nx, int ny, int Nx){
    
}
