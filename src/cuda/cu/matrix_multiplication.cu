#include <iostream>
#include <string>
#include <iomanip>

#include <string.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdlib.h>
#include <unistd.h>

#include <cuda_runtime.h>  // For CUDA runtime API
#include <helper_cuda.h>  // For checkCudaError macro
#include <helper_timer.h>  // For CUDA SDK timers


#define BLOCK_SIZE 32


/**************************** Logger ********************************/
#define LOG_MESSAGE_SIZE 256

void logger(const char* tag, const char* message) {
   time_t now;
   time(&now);
   printf("%s [%s]: %s\n", ctime(&now), tag, message);
}


void logger_info(const char* message) {
    logger((char *)"INFO",message);
}


void logger_error(const char* message) {
    logger((char *)"ERROR",message);
}


void logger_debug(const char* message) {
    logger((char *)"DEBUG",message);
}

/**********************************************************************/

/**************************** Formatter ********************************/
void getFormattedTime(double seconds, char *formatted_string)
{
    double s(fabs(seconds));
    int h(s/3600);
    int min(s/60 - h*60);
    double sec(s - (h*60 + min)*60);
    std::ostringstream oss;
    oss<<std::setfill('0')<<std::setw(2)<<fabs(seconds)/seconds*h<<":"<<std::setw(2)<<min<<":";
    if (sec/10<1)
        oss<<"0";
    oss<<sec;
    strcpy(formatted_string,oss.str().c_str());
}
/**********************************************************************/


/**************************** FILE I/O ********************************/
#ifndef DATA_DIR
#define DATA_DIR "../../data/"
#endif





static int create_dir(char *rel_path,mode_t mode) {
  return mkdir(rel_path, mode);
}

// Return bytes written in file
int write_matrix_to_file(float *matrix, int rows, int cols, char *matrix_name){
    int ret = 0;
    struct stat st;
    //mode_t mode = st.st_mode & (S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP);
    FILE *matrix_file;
    char file_name[63];
    char dir_name[32]; 
    char log_string[LOG_MESSAGE_SIZE];

    sprintf(dir_name,DATA_DIR);

    if ((ret = stat(dir_name, &st)) == -1) { // If directory not found
        create_dir(dir_name,0770); // TODO: Check for the right parametric bitmask
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Created data directory, with permissions: -rw-r--r-- %s\n",dir_name);
        logger_info(log_string);
    }

    sprintf(&(dir_name)[strlen(dir_name)],"square/matrix%dx%dx%d/",rows,cols,cols);

    if ((ret = stat(dir_name, &st)) == -1) { // If directory not found
        create_dir(dir_name,0770); // TODO: Check for the right parametric bitmask
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Created data directory, with permissions: -rw-r--r-- %s\n",dir_name);
        logger_info(log_string);
    }

    sprintf(file_name,dir_name);

    sprintf(&(file_name)[strlen(dir_name)],"matrix_%s_%dx%d.bin",matrix_name,rows,cols);

    // Open or create file
	if((matrix_file=fopen(file_name, "w+"))==NULL) {
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error opening file %s\n",file_name);
		logger_info(log_string);
		exit(EXIT_FAILURE);
	}

    // Write rows number as first element of the binary file
    ret = fwrite(&rows,sizeof(int),1,matrix_file);
    if(ret == 0){
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error writing bytes on file %s\n",file_name);
		logger_info(log_string);
        exit(EXIT_FAILURE);
    }
    
    // Write cols number as second element of the binary file
    ret = fwrite(&cols,sizeof(int),1,matrix_file);
    if(ret == 0){
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error opening file %s\n",file_name);
		logger_info(log_string);
        exit(EXIT_FAILURE);
    }

    ret = fwrite(matrix,sizeof(float),rows*cols,matrix_file);
    if(ret == 0){
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error opening file %s\n",file_name);
		logger_info(log_string);
        exit(EXIT_FAILURE);
    }
    

    fclose(matrix_file);


    memset(log_string,0,LOG_MESSAGE_SIZE);
    sprintf(log_string,"File %s populated correctly!\n",file_name);
	logger_info(log_string);

    return sizeof(float)*cols*rows + sizeof(int) + sizeof(int);
}


void read_matrix_from_file(float *matrix,int rows_expected, int cols_expected, char *matrix_name){
    FILE *matrix_file;
    char file_name[64];
    char log_string[LOG_MESSAGE_SIZE];
    int rows;
    int cols;
    int ret = 0;

    sprintf(file_name,DATA_DIR);
    sprintf(&(file_name)[strlen(file_name)],"square/matrix%dx%dx%d/matrix_%s_%dx%d.bin",rows_expected,rows_expected,cols_expected,matrix_name,rows_expected,cols_expected);

    // Open or create file
	  if((matrix_file=fopen(file_name, "r"))==NULL) {
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error opening file %s\n",file_name);
		  logger_error(log_string);
		  exit(EXIT_FAILURE);
	  }

    // Read rows number as first element of the binary file
    ret = fread(&rows,sizeof(int),1,matrix_file);
    if(ret == 0){
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error writing bytes on file %s\n",file_name);
        logger_error(log_string); //TODO: Check with errno
        exit(EXIT_FAILURE);
    }else if(rows != rows_expected){    
        logger_error("Matrix passed to the function must be choerent with the size of the matrix in the file!\n");
        exit(EXIT_FAILURE);
    }

    // Read cols number as first element of the binary file
    ret = fread(&cols,sizeof(int),1,matrix_file);
    if(ret == 0){
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error writing bytes on file %s\n",file_name);
        logger_error(log_string); //TODO: Check with errno
        exit(EXIT_FAILURE);
    }else if(cols != cols_expected){
        logger_error("Matrix passed to the function must be choerent with the size of the matrix in the file!\n");
        exit(EXIT_FAILURE);
    }


    ret = fread(matrix,sizeof(float),rows*cols,matrix_file);
    if(ret == 0 && errno == EOF){
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error reading bytes from file %s to matrix! Error returned: %s\n",file_name, strerror(errno));
        logger_error(log_string); 
        exit(EXIT_FAILURE);
    }
    


    fclose(matrix_file);
}


void write_cuda_stats(int N, int K, int M,double total_time,float max_diff,float max_rel_diff,float gflops){
    int ret = 0;
    struct stat st;
    FILE *matrix_file;
    char *file_name = (char *)"data/stats.csv";
    char *dir_name = (char *)"data/"; 
    char log_string[LOG_MESSAGE_SIZE];

    // Create directory if needed
    if ((ret = stat(dir_name, &st)) == -1) { // If directory not found
        create_dir(dir_name,0770); 
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Created data directory, with permissions: -rw-r--r-- %s\n",dir_name);
        logger_info(log_string);
    }


    // Open or create file
	  if((matrix_file=fopen(file_name, "a"))==NULL) {
      memset(log_string,0,LOG_MESSAGE_SIZE);
      sprintf(log_string,"Error opening file %s\n",file_name);
		  logger_info(log_string);
		  exit(EXIT_FAILURE);
	  } 

    fprintf(matrix_file,"%d,%d,%d,%f,%f,%f,%f\n",N,K,M,total_time,max_diff,max_rel_diff,gflops);

    fclose(matrix_file);
}


#ifdef DEBUG


/**********************************************************************/

/**
 * @brief
 * @param matrix
 * @param rows 
 * @param cols
*/
void printMatrix(const float *matrix,int rows,int cols){
  std::cout << "______________________________________________________________________________________________________________________________________________________\n";
  for (int row = 0; row < rows; ++row) {
    std::cout << "| ";
    for (int col = 0; col < cols; ++col) {
      int idx = row * cols + col;
      std::cout << matrix[idx] << " ";
    }
    std::cout << "\t\t|\n";
  }
  std::cout << "______________________________________________________________________________________________________________________________________________________\n\n";
}


/**
 * @brief An implementation of the sequential matrix by matrix multiplication, taking in count memory acces
 * and compiled with -O3. It's used only for debugging and for the developement phase of the project.
 * @param matrix_A The memory in which are the values stored in the matrix A N x K
 * @param matrix_B The memory in which are the values stored in the matrix B K x M
 * @param matrix_C The memory in which are the values stored in the matrix C N x M
 * @param N Rows dimwnsion of the matrix A and C
 * @param K Colums dimension of the matrix A and rows dimension of the matrix B
 * @param M Columns dimension of the matrix B and C
*/
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
#endif


__global__ void gpuMatrixVectorSharedMemory(float* A, float* B, float* C, int N, int K, int M)
{
    __shared__ float matrix_A_shared[BLOCK_SIZE][BLOCK_SIZE]; // Memory in which is stored the matrix A shared between all the threads in the tile
    __shared__ float matrix_B_shared[BLOCK_SIZE][BLOCK_SIZE];
    int col = blockIdx.x*blockDim.x + threadIdx.x;
    int row = blockIdx.y*blockDim.y + threadIdx.y;
    int tx = threadIdx.x;
    int ty = threadIdx.y;
    float sum = 0;
    int num_blocks = (int) (N/BLOCK_SIZE)

    for(int i=0; i<num_blocks; ++i)
    {
        // Load data in the tile of shared memory
        if(row<M && t*BLOCK_SIZE+tx<N)
            matrix_A_shared[ty][tx] = A[row*N + i*BLOCK_SIZE+tx];
        else
            matrix_A_shared[ty][tx] = 0.0;
        if(i*BLOCK_SIZE+ty<N && col<K)
            matrix_B_shared[ty][tx] = B[(t*BLOCK_SIZE+ty)*K + col];
        else
            matrix_B_shared[ty][tx] = 0.0;
        __syncthreads();

        // Actual computation over the single tile
        for(int j=0; j<BLOCK_SIZE; ++j)
            sum += matrix_A_shared[ty][j] * matrix_B_shared[j][tx];
        __syncthreads();
    }
    if(row<M && col<K)
        C[col+row*K] += sum;

}



__global__ void gpuMatrixVector(float* A, float* B, float* C, int N, int K, int M)
{
  int col = blockIdx.x * blockDim.x + threadIdx.x;
  int row = blockIdx.y * blockDim.y + threadIdx.y;
  
  if (row < N && col < M) {
    for (int i = 0; i < K; ++i) {
      C[row * M + col]  += A[row * K + i] * B[i * M + col];
      }
  }

}

int main(int argc, char** argv) {
  int N = 0;
  int K = 0;
  int M = 0;

	if(argc != 4){
        logger_error("Program should be called like ./<elf file name> <N> <K> <M>");
        exit(EXIT_FAILURE);
    }

    if (sscanf (argv[1], "%d", &N) != 1) {
        fprintf(stderr, "[ERROR] - not an integer");
        logger_error("Program should be called like ./<elf file name> <N> <K> <M>");
        exit(EXIT_FAILURE);
    }

    if (sscanf (argv[2], "%d", &K) != 1) {
        fprintf(stderr, "[ERROR] - not an integer");
        logger_error("Program should be called like ./<elf file name> <N> <K> <M>");
        exit(EXIT_FAILURE);
    }

	if (sscanf (argv[3], "%d", &M) != 1) {
        fprintf(stderr, "[ERROR] - not an integer");
        logger_error("Program should be called like ./<elf file name> <N> <K> <M>");
        exit(EXIT_FAILURE);
    }


  
  // ----------------------- Host memory initialisation ----------------------- //
  float* h_A = new float[N * K];
  float* h_B = new float[K * M];
  float* h_C = new float[N * M];

  float* h_C_result_host = new float[N * M];
  float* h_C_result_device = new float[N * M];

  float time;
  char formatted_time_string[64];
  cudaEvent_t start, stop;

  dim3 blockDim(BLOCK_SIZE,BLOCK_SIZE); // Don't need to write z = 1, max 1024
  int gx = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;//(m % blockDim.x==0) ? m / blockDim.x : m / blockDim.x + 1;
  int gy = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;//(n % blockDim.y==0) ? n / blockDim.y : n / blockDim.y + 1;
  dim3 gridDim(gx, gy);
  printf("%d %d\n",gx,gy);

  read_matrix_from_file(h_A,N,K,(char *)"A");
  read_matrix_from_file(h_B,K,M,(char *)"B");
  read_matrix_from_file(h_C,N,M,(char *)"C");
  read_matrix_from_file(h_C_result_host,N,M,(char *)"sequential_C");

  std::cout << "[INFO] Matrix initialization done.\n";

#ifdef DEBUG
  puts("\t[MATRIX] Matrix A");
  printMatrix(h_A,N,K);
  puts("\t[MATRIX] Matrix B");
  printMatrix(h_B,K,M);
  puts("\t[MATRIX] Matrix C");
  printMatrix(h_C,N,M);
#endif
  std::cout << "Test case: " << N  << " x " << M << std::endl;
// ---------------------- Device memory initialisation ---------------------- //

  float *d_A, *d_B, *d_C;

  checkCudaErrors(cudaMalloc((void**) &d_A, N * K * sizeof(float)));
  checkCudaErrors(cudaMalloc((void**) &d_B, K * M * sizeof(float)));
  checkCudaErrors(cudaMalloc((void**) &d_C, N * M * sizeof(float)));

  // Copy matrices from the host (CPU) to the device (GPU).
  checkCudaErrors(cudaMemcpy(d_A, h_A, N * K * sizeof(float), cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_B, h_B, K * M * sizeof(float), cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_C, h_C, N * M * sizeof(float), cudaMemcpyHostToDevice));

#ifdef DEBUG
  // ------------------------ Calculations on the CPU ------------------------- //

  // Create the CUDA SDK timer.
  StopWatchInterface* timer = 0;
  sdkCreateTimer(&timer);

  timer->start();
  
  CpuMatrixVector(h_A, h_B, h_C_result_host, N, K ,M);

  timer->stop();
  double cpuflops=flopcnt/ timer->getTime();
  std::cout << "  CPU time: " << timer->getTime() << " ms." << " MFLOPS " << cpuflops << std::endl;

  puts("\t[RESULT] Sequential result");
  printMatrix(h_C_result_host,N,M);
#endif

// ------------------------ Calculations on the GPU ------------------------- //
  
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);


  gpuMatrixVector<<<gridDim, blockDim >>>(d_A, d_B, d_C, N, K, M);
  checkCudaErrors(cudaDeviceSynchronize());

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time, start, stop);

  std::cout << "Time to generate: " << time << " ms\n";


  double gpuflops= 2.0*N*K*M / time;
  memset(formatted_time_string, 0, strlen(formatted_time_string));
  getFormattedTime(time/1000.0f,formatted_time_string);
  std::cout << "  GPU time: " << formatted_time_string << " MFLOPS " << gpuflops<<std::endl;

  // Download the resulting vector d_y from the device and store it in h_y_d.
  checkCudaErrors(cudaMemcpy(h_C_result_device, d_C, N * M * sizeof(float),cudaMemcpyDeviceToHost));

#ifdef DEBUG
  puts("\t[RESULT] CUDA result");
  printMatrix(h_C_result_device,N,M);
#endif
  // Now let's check if the results are the same.
  float reldiff = 0.0f;
  float diff = 0.0f;
  
  for (int row = 0; row < N; ++row) {
    for(int col = 0; col < M; ++col){
      float maxabs = std::max(std::abs(h_C_result_device[row * M + col]),std::abs(h_C_result_host[row * M + col]));
      if (maxabs == 0.0){
        maxabs=1.0;
      } 
      reldiff = std::max(reldiff, std::abs(h_C_result_device[row * M + col] - h_C_result_host[row * M + col])/maxabs);
      diff = std::max(diff, std::abs(h_C_result_device[row * M + col] - h_C_result_host[row * M + col]));
    }
  }
  std::cout << "Max diff =  " << diff << "  Max rel diff = " << reldiff << std::endl;
  fflush(stdout);
  write_cuda_stats(N, K, M,time,diff,reldiff,gflops);


// ------------------------------- Cleaning up ------------------------------ //

  write_matrix_to_file(h_C_result_device,N,M,(char *)"cuda_C");
#ifdef DEBUG
  delete timer;
#endif
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
