// i.e. A[i, j] is stored in i * ncols + j element of the vector.
//

#include <iostream>

#include <cuda_runtime.h>  // For CUDA runtime API
#include <helper_cuda.h>  // For checkCudaError macro
#include <helper_timer.h>  // For CUDA SDK timers


// TODO What is a good initial guess for XBD and YBD (both
// greater than 1) ?
// After you get the code to work, experiment with different sizes
// to find the best possible performance
// Note: For meaningful time measurements you need sufficiently large matrices.
#define XBD 16
#define YBD 16

const int TILE_WIDTH = 32;

const dim3 BLOCK_DIM(TILE_WIDTH, TILE_WIDTH, 1);




void printMatrix(const float *h_A,int ROWS,int COLS){
  std::cout << "______________________________________________________________________________________________________________________________________________________\n";
  for (int row = 0; row < ROWS; ++row) {
    std::cout << "| ";
    for (int col = 0; col < COLS; ++col) {
      int idx = row * COLS + col;
      std::cout << h_A[idx] << " ";
    }
    std::cout << "\t\t|\n";
  }
  std::cout << "______________________________________________________________________________________________________________________________________________________\n\n";
}



// Simple CPU implementation of matrix addition.
void CpuMatrixVector(float *matrix_A, float *matrix_B, float *matrix_C, size_t N, size_t K, size_t M) {
  for (size_t i=0;i<N;i++) 
        for(size_t j=0;j<M;j++) {
            for(size_t k=0;k<K;k++){
                matrix_C[i * M + j]=matrix_C[i * M + j]+matrix_A[i * K + k]*matrix_B[k * M + j];
			}
    }
}

__global__ void gpuMatrixVector(float* A_d, float* B_d, float* C_d, int m, int k, int n)
{
    __shared__ float ds_A[TILE_WIDTH][TILE_WIDTH];
    __shared__ float ds_B[TILE_WIDTH][TILE_WIDTH];
    int col = blockIdx.x*blockDim.x + threadIdx.x;
    int row = blockIdx.y*blockDim.y + threadIdx.y;
    int tx = threadIdx.x;
    int ty = threadIdx.y;
    float sum = 0;

    for(int t=0; t<(n-1)/TILE_WIDTH+1; t++)
    {
        if(row<m && t*TILE_WIDTH+tx<n)
            ds_A[ty][tx] = A_d[row*n + t*TILE_WIDTH+tx];
        else
            ds_A[ty][tx] = 0.0;
        if(t*TILE_WIDTH+ty<n && col<k)
            ds_B[ty][tx] = B_d[(t*TILE_WIDTH+ty)*k + col];
        else
            ds_B[ty][tx] = 0.0;
        __syncthreads();
        for(int i=0; i<TILE_WIDTH; i++)
            sum += ds_A[ty][i] * ds_B[i][tx];
        __syncthreads();
    }
    if(row<m && col<k)
        C_d[col+row*k] += sum;
}


int main(int argc, char** argv) {

  if (argc != 4) {
    fprintf(stderr,"[ERROR] Correct usage: %s n k m\nYou should insert matrices dimensons!\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  size_t n = (size_t) atoi(argv[1]);
  size_t k = (size_t) atoi(argv[2]);
  size_t m = (size_t) atoi(argv[3]);

  
  
  // ----------------------- Host memory initialisation ----------------------- //

  float* h_A = new float[n * k];
  float* h_B = new float[k * m];
  float* h_C = new float[n * m];

  float* h_C_result_host = new float[n * m];
  float* h_C_result_device = new float[n * m];

  srand(123456);
  for (size_t row = 0; row < n; ++row) {
    for (size_t col = 0; col < m; ++col) {
      size_t idx = row * m + col;
      h_A[idx] = 100.0f * static_cast<float>(rand()) / RAND_MAX;
      h_B[idx] = 100.0f * static_cast<float>(rand()) / RAND_MAX;
      h_C[idx] = 100.0f * static_cast<float>(rand()) / RAND_MAX;
    }
  }

  memcpy(h_C_result_host,h_C,n * m * sizeof(float));

/*
  puts("\t[MATRIX] Matrix A");
  printMatrix(h_A,n,k);
  puts("\t[MATRIX] Matrix B");
  printMatrix(h_B,k,m);
  puts("\t[MATRIX] Matrix C");
  printMatrix(h_C,n,m);
  puts("\t[MATRIX] Matrix h_C_result_host");
  printMatrix(h_C_result_host,n,m);*/

  std::cout << "Test case: " << n  << " x " << m << std::endl;
// ---------------------- Device memory initialisation ---------------------- //

  float *d_A, *d_B, *d_C;

  checkCudaErrors(cudaMalloc((void**) &d_A, n * k * sizeof(float)));
  checkCudaErrors(cudaMalloc((void**) &d_B, k * m * sizeof(float)));
  checkCudaErrors(cudaMalloc((void**) &d_C, n * m * sizeof(float)));

  // Copy matrices from the host (CPU) to the device (GPU).
  checkCudaErrors(cudaMemcpy(d_A, h_A, n * k * sizeof(float), cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_B, h_B, k * m * sizeof(float), cudaMemcpyHostToDevice));

  // ------------------------ Calculations on the CPU ------------------------- //
  float flopcnt=2.e-6*n*m;
  
  // Create the CUDA SDK timer.
  StopWatchInterface* timer = 0;
  sdkCreateTimer(&timer);

  timer->start();
  CpuMatrixVector(h_A, h_B, h_C_result_host, n, k ,m);

  timer->stop();
  float cpuflops=flopcnt/ timer->getTime();
  std::cout << "  CPU time: " << timer->getTime() << " ms." << " GFLOPS " << cpuflops << std::endl;
/*
  puts("\t[RESULT] Sequential result");
  printMatrix(h_C_result_host,n,m);*/


// ------------------------ Calculations on the GPU ------------------------- //

  // TODO Calculate the dimension of the grid of blocks. A 1D grid suffices.
  
  dim3 GRID_DIM((k-1)/TILE_WIDTH+1, (m-1)/TILE_WIDTH+1, 1);

  timer->reset();
  timer->start();
  gpuMatrixVector<<<GRID_DIM, BLOCK_DIM >>>(d_A, d_B, d_C, n, m, k);
  checkCudaErrors(cudaDeviceSynchronize());

  timer->stop();
  float gpuflops=flopcnt/ timer->getTime();
  std::cout << "  GPU time: " << timer->getTime() << " ms." << " GFLOPS " << gpuflops<<std::endl;

  // Download the resulting vector d_y from the device and store it in h_y_d.
  checkCudaErrors(cudaMemcpy(h_C_result_device, d_C, n * m * sizeof(float),cudaMemcpyDeviceToHost));
/*
  puts("\t[RESULT] CUDA result");
  printMatrix(h_C_result_device,n,m);*/

  // Now let's check if the results are the same.
  float reldiff = 0.0f;
  float diff = 0.0f;
  
  for (int row = 0; row < n; ++row) {
    for(int col = 0; col < m; ++col){
      float maxabs = std::max(std::abs(h_C_result_device[row * m + col]),std::abs(h_C_result_host[row * m + col]));
      if (maxabs == 0.0){
        maxabs=1.0;
      } 
      reldiff = std::max(reldiff, std::abs(h_C_result_device[row * m + col] - h_C_result_host[row * m + col])/maxabs);
      diff = std::max(diff, std::abs(h_C_result_device[row * m + col] - h_C_result_host[row * m + col]));
    }
  }
  std::cout << "Max diff = " << diff << "  Max rel diff = " << reldiff << std::endl;
  // Rel diff should be as close as possible to unit roundoff; float
  // corresponds to IEEE single precision, so unit roundoff is
  // 1.19e-07
  // 

// ------------------------------- Cleaning up ------------------------------ //

  delete timer;

  checkCudaErrors(cudaFree(d_A));
  checkCudaErrors(cudaFree(d_B));
  checkCudaErrors(cudaFree(d_C));

  //delete[] h_A;
  //delete[] h_B;
  //delete[] h_C;
  //delete[] h_C_result_host;
  //delete[] h_C_result_device;
  return 0;
}
