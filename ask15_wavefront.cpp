#include "wavefront.h"
#include "omp.h"
#include <algorithm> //ask about this 
#define _USE_MATH_DEFINES
using std::min;
using std::max;
#include <cmath>
#include <math.h>

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
    
//     if(Nx==1 && Ny==1){
//       counterx=nx+1;
//       countery=ny+1;
//     }
//     else{  
//     if(I<Nx-1){
//         counterx=(nx+1)/(Nx-1); //1280
//     }
//     else{
//         counterx=(nx+1)%(Nx-1); //1
//     }
//     if(J<Ny-1){
//         countery=(ny+1)/(Ny-1); //1280
//     }
//     else{
//         countery=(ny+1)%(Ny-1); //1
//     }    
//     }
    for(int i=I*counterx;i<I*counterx+counterx;i++){
        for(int j=J*countery;j<J*countery+countery;j++){
           if(i==0){
            data[cartesian2flat(i,j,ny+1)]=sin(2*M_PI*j/ny);
            //printf("2 %d %d \t",i,j);
           }
           else if(j==0){
            data[cartesian2flat(i,j,ny+1)]=sin(2*M_PI*i/nx);
            //printf("2");
           }
           else{
             data[cartesian2flat(i,j,ny+1)]=Cx*data[cartesian2flat(i-1,j,ny+1)]+Cy*data[cartesian2flat(i,j-1,ny+1)];
             //printf("3");
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
        for(int I=std::max(0,w-Nx+1);I<=std::min(w,Nx);I++){
            J=min(w-I,Ny-1)-counter;
            //printf("1");
            //printf("thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads());
            process_block(data, I, J,nx, ny, Nx, Ny);
            counter=counter+1;
        }
    }
}


void wavefront420(double* data, int nx, int ny, int Nx){
    
}
