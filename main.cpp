#include <chrono>
#include <cstdio>
#include <string> 

#include <omp.h>

#include "wavefront.h"

using std::printf;
using std::scanf;
using std::chrono::steady_clock;
using std::micro;
using std::stoi;

const int CAAM_420 = 420;
const int CAAM_520 = 520;

void time520(double* data, int nx, int ny, int Nx, int Ny);

int main(int argc, char* argv[]) {
    // Parse arguments
    if (argc < 3 ) {
        printf("Not enough arguments supplied; exiting program...\n");
        return -1;
    } 
    
    int choice = stoi(argv[1], NULL);
    int Nx = 0;
    int Ny = 0;
    if (argc == 3 && choice == CAAM_420) {
        Nx = stoi(argv[2], NULL);
    } else if (argc == 4 && choice == CAAM_520) {
        Nx = stoi(argv[2], NULL);
        Ny = stoi(argv[3], NULL);
    } else {
        printf("Incorrect number of arguments supplied for option %d; exiting program...\n", choice);
        return -1;
    }

    // Allocate memory
    double* data = NULL;
    data = new double[(ny+1)*(nx+1)];
    if (data == NULL) {
        printf("Could not allocate enough memory for data; exiting program... \n");
        return -2;
    }

    // Call the appropriate timing function
    if (choice == CAAM_420){
       time420(data, nx, ny, Nx);
    }
    else if (choice == CAAM_520){
        time520(data, nx, ny, Nx, Ny);
    }
    else {
        printf("Invalid choice entered; exiting program...\n");
        return -3;
    }

    return 0;
}


void time420(double* data, int nx, int ny, int Nx) {
     steady_clock::time_point start = steady_clock::now();
     wavefront420(data, nx, ny, Nx);
     steady_clock::time_point end = steady_clock::now();
     auto diff = end - start;

     printf("Total time (420): \n    Nx = %d\n    T (microseconds): %.5lf\n", Nx, 
         std::chrono::duration<double, micro>(diff).count() );
 }


void time520(double* data, int nx, int ny, int Nx, int Ny) {
    steady_clock::time_point start = steady_clock::now();
    wavefront520(data, nx, ny, Nx, Ny);
    steady_clock::time_point end = steady_clock::now();
    auto diff = end - start;

    printf("Total time (520): \n    Nx = %d, Ny = %d\n    T (microseconds): %.4lf\n", Nx, Ny,
        std::chrono::duration<double, micro>(diff).count() );

}
