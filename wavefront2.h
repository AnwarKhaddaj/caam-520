#ifndef WAVEFRONT_H
#define WAVEFRONT_H

const int nx = 6400;
const int ny = 6400;

const double dx = 1.0 / nx;
const double dy = 1.0 / ny;

const double Cx = dy / (dx + dy);
const double Cy = dx / (dx + dy);

inline int cartesian2flat(int i, int j, int ny) {
    return i*ny + j;
}

void process_block(double* data, int I, int J, int nx, int ny, int Nx, int Ny);

void wavefront420(double* data, int nx, int ny, int Nx); // Ny will = the number of threads

void wavefront520(double* data, int nx, int ny, int Nx, int Ny);

#endif