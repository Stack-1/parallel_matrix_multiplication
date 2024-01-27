// i.e. A[i, j] is stored in i * ncols + j element of the vector.
//

#include <iostream>

#include <cuda_runtime.h>  // For CUDA runtime API
#include <helper_cuda.h>  // For checkCudaError macro
#include <helper_timer.h>  // For CUDA SDK timers

#define BLOCK_SIZE 32

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
void CpuMatrixVector(float *matrix_A, float *matrix_B, float *matrix_C, int N, int K, int M) {
    for (int i = 0; i < N; ++i) 
    {
        for (int h = 0; h < K; ++h)
        {
            for (int j = 0; j < M; ++j) 
            {
                matrix_C[i * M + j] += matrix_A[i * K + h] * matrix_B[h * M + j];
            }
        }
    }

}

__global__ void gpuMatrixVector(float* A, float* B, float* C, int n, int k, int m)
{
  int col = blockIdx.x * blockDim.x + threadIdx.x;
  int row = blockIdx.y * blockDim.y + threadIdx.y;
  
  if (row < n && col < m) {
    for (int i = 0; i < k; ++i) {
      C[row * m + col]  += A[row * k + i] * B[i * m + col];
      }
  }

}


int main(int argc, char** argv) {

  if (argc != 4) {
    fprintf(stderr,"[ERROR] Correct usage: %s n k m\nYou should insert matrices dimensons!\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  int n = (int) atoi(argv[1]);
  int k = (int) atoi(argv[2]);
  int m = (int) atoi(argv[3]);

  
  
  // ----------------------- Host memory initialisation ----------------------- //

  float* h_A = new float[n * k];
  float* h_B = new float[k * m];
  float* h_C = new float[n * m];

  float* h_C_result_host = new float[n * m];
  float* h_C_result_device = new float[n * m];

  float time;
  cudaEvent_t start, stop;

  dim3 blockDim(BLOCK_SIZE,BLOCK_SIZE); // Don't need to write z = 1, max 1024
  int gx = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;//(m % blockDim.x==0) ? m / blockDim.x : m / blockDim.x + 1;
  int gy = (m + BLOCK_SIZE - 1) / BLOCK_SIZE;//(n % blockDim.y==0) ? n / blockDim.y : n / blockDim.y + 1;
  dim3 gridDim(gx, gy);
  printf("%d %d\n",gx,gy);

  srand(123456);
  for (int row = 0; row < n; ++row) {
    for (int col = 0; col < m; ++col) {
      int idx = row * m + col;
      h_C[idx] = 100.0f * static_cast<float>(rand()) / RAND_MAX;
    }
  }


  for (int row = 0; row < n; ++row) {
    for (int col = 0; col < k; ++col) {
      int idx = row * k + col;
      h_A[idx] = 100.0f * static_cast<float>(rand()) / RAND_MAX;
    }
  }

  for (int row = 0; row < k; ++row) {
    for (int col = 0; col < m; ++col) {
      int idx = row * m + col;
      h_B[idx] = 100.0f * static_cast<float>(rand()) / RAND_MAX;
    }
  }

  std::cout << "[INFO] Matrix initialization done.\n";

  memcpy(h_C_result_host,h_C,n * m * sizeof(float));
#ifdef DEBUG
  puts("\t[MATRIX] Matrix A");
  printMatrix(h_A,n,k);
  puts("\t[MATRIX] Matrix B");
  printMatrix(h_B,k,m);
  puts("\t[MATRIX] Matrix C");
  printMatrix(h_C,n,m);
#endif
  std::cout << "Test case: " << n  << " x " << m << std::endl;
// ---------------------- Device memory initialisation ---------------------- //

  float *d_A, *d_B, *d_C;

  checkCudaErrors(cudaMalloc((void**) &d_A, n * k * sizeof(float)));
  checkCudaErrors(cudaMalloc((void**) &d_B, k * m * sizeof(float)));
  checkCudaErrors(cudaMalloc((void**) &d_C, n * m * sizeof(float)));

  // Copy matrices from the host (CPU) to the device (GPU).
  checkCudaErrors(cudaMemcpy(d_A, h_A, n * k * sizeof(float), cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_B, h_B, k * m * sizeof(float), cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_C, h_C, n * m * sizeof(float), cudaMemcpyHostToDevice));

  // ------------------------ Calculations on the CPU ------------------------- //
  float flopcnt=2.e-3*n*m;
  
  // Create the CUDA SDK timer.
  StopWatchInterface* timer = 0;
  sdkCreateTimer(&timer);

  timer->start();
  
  CpuMatrixVector(h_A, h_B, h_C_result_host, n, k ,m);

  timer->stop();
  float cpuflops=flopcnt/ timer->getTime();
  std::cout << "  CPU time: " << timer->getTime() << " ms." << " MFLOPS " << cpuflops << std::endl;

#ifdef DEBUG
  puts("\t[RESULT] Sequential result");
  printMatrix(h_C_result_host,n,m);
#endif

// ------------------------ Calculations on the GPU ------------------------- //

  // TODO Calculate the dimension of the grid of blocks. A 1D grid suffices.
  
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);


  gpuMatrixVector<<<gridDim, blockDim >>>(d_A, d_B, d_C, n, k, m);
  checkCudaErrors(cudaDeviceSynchronize());

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time, start, stop);

  std::cout << "Time to generate: " << time << " ms\n";


  float gpuflops=flopcnt/ time;
  std::cout << "  GPU time: " << time << " ms." << " MFLOPS " << gpuflops<<std::endl;

  // Download the resulting vector d_y from the device and store it in h_y_d.
  checkCudaErrors(cudaMemcpy(h_C_result_device, d_C, n * m * sizeof(float),cudaMemcpyDeviceToHost));

#ifdef DEBUG
  puts("\t[RESULT] CUDA result");
  printMatrix(h_C_result_device,n,m);
#endif
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
  std::cout << "Max diff =  " << diff << "  Max rel diff = " << reldiff << std::endl;
  fflush(stdout);
  // Rel diff should be as close as possible to unit roundoff; float
  // corresponds to IEEE single precision, so unit roundoff is
  // 1.19e-07
  // 

// ------------------------------- Cleaning up ------------------------------ //

  delete timer;

  checkCudaErrors(cudaFree(d_A));
  checkCudaErrors(cudaFree(d_B));
  checkCudaErrors(cudaFree(d_C));

  delete[] h_A;
  delete[] h_B;
  delete[] h_C;
  delete[] h_C_result_host;
  delete[] h_C_result_device;
  return 0;
}
