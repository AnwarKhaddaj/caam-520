const int n_rows = 1280;
const int n_cols = 320;
const int N_cols = n_cols / 32;
const int N = 64;

__global__ void matrix_scale_wrapped(float* A, float* B, float* D, const int n_rows, const int n_cols, const int N);

__device__ int cartesian_2_flat(int i, int j, int stride) { //helper function
    return i*stride+j;
}

void launch_gpu_matrix_scale(float* A_cpu, float* B_cpu, float* D_cpu) {
    // Pointers for memory on the GPU
    float* A_gpu;
    float* D_gpu;
    float* B_gpu;

    // Allocate memory on the GPU
    cudaMalloc(&A_gpu,n_rows*n_cols*sizeof(float));
    cudaMalloc(&D_gpu,n_cols*sizeof(float));
    cudaMalloc(&B_gpu,n_rows*n_cols*sizeof(float));

    // Copy memory to the GPU
    cudaMemcpy(A_gpu,A_cpu,n_rows*n_cols*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(D_gpu,D_cpu,n_cols*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(B_gpu,B_cpu,n_rows*n_cols*sizeof(float),cudaMemcpyHostToDevice);

    // Compute how much shared memory is needed
    int shared_mem = 32*sizeof(float);

    // Compute the size of the CUDA grid (i.e. the number of blocks in each dim)
    int threads_per_block = 32;
    dim3 number_of_blocks(n_rows/N,N_cols);

    // Invoke the CUDA kernel
    matrix_scale_wrapped<<<number_of_blocks,threads_per_block,shared_mem>>>(A_gpu, B_gpu, D_gpu, n_rows, n_cols, N);

    // Copy memory back to the CPU
    cudaMemcpy(A_cpu,A_gpu,n_rows*n_cols*sizeof(float),cudaMemcpyDeviceToHost);
    cudaMemcpy(B_cpu,B_gpu,n_rows*n_cols*sizeof(float),cudaMemcpyDeviceToHost);
    cudaMemcpy(D_cpu,D_gpu,n_cols*sizeof(float),cudaMemcpyDeviceToHost);

    // Free memory on the GPU
    cudaFree(A_gpu);
    cudaFree(B_gpu);
    cudaFree(D_gpu);
}

__global__ void matrix_scale_wrapped(float* A, float* B, float* D, const int n_rows, const int n_cols, const int N) {
    
    int i = threadIdx.x;
    int I = blockIdx.x;
    int J = blockIdx.y;
    
    int Dy = blockDim.y;

    __shared__ float shared_mem[32]; //shared memory
    shared_mem[i] = D[i+J*Dy]; //copying memory from a chunk of D to shared memory

    int el_ind;  //index of the element in A that the thread needs to access
    for (int warp_jump = 0; warp_jump < N; warp_jump++) {
        el_ind = cartesian_2_flat(I*N+warp_jump,J*Dy+i,n_cols); //jumping by one warp as threads access vertical elements
        B[el_ind] = A[el_ind]*shared_mem[warp_jump];
    }
}

int main() {
    // Pointers for the matrices/arrays
    float* A;
    float* B;
    float* D;

    // Allocate the arrays
    A = new float[n_rows*n_cols];
    B = new float[n_rows*n_cols];
    D = new float[n_cols];

    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
            A[j + i*n_cols] = 1;
        }
    }

    for(int j = 0; j < n_cols; j++)
        D[j] = j;

    launch_gpu_matrix_scale(A, B, D);

    // Check B
    int correct = 0;
    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
            if (B[j + i *n_cols] != j)
                correct = -1;
        }
    }

    delete[] A;
    delete[] B;
    delete[] D;
    
    return correct;
}