#ifndef __MPI__COMPUTATION__
#define __MPI__COMPUTATION__

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include "../../utils/logger/logger.h"
#include "../../utils/string_formatter/formatter.h"
#include "../../utils/file_io/stackio.h"

#define BLOCK_ROWS 30
#define BLOCK_COLS 30
#define AUDIT if(my_rank == 0)

#define max(a,b) \
    ({ __typeof__ (a) _a = (a); \
        __typeof__ (b) _b = (b); \
        _a > _b ? _a : _b; })

typedef struct computed_dimensions{
    int rows;
    int cols;
}computed_dimensions;

#ifndef FORMATTED_STRING_SIZE
#define FORMATTED_STRING_SIZE 256
#endif

#define LOG_MESSAGE_SIZE 512


void mpi_matrix_multiplication(float *matrix_A, float *matrix_B, float **matrix_C, int N, int K, int M, int my_rank);
void reduce_partial_results(float *matrix, float **return_matrix, int rows, int cols, int process_grid_cols,int my_rank, MPI_Comm comm);
computed_dimensions *compute_submatrix_dimension(int matrix_rows, int matrix_cols, int block_rows, int block_cols, int process_gird_rows, int process_gird_cols, int my_rank);
void generate_block_cyclic_distribution_matrix_A(char *matrix_file_name, float **matrix_A,int N, int K, int N_subarray,int K_subarray,int block_rows, int block_cols, int process_gird_rows, int process_gird_cols, int my_rank,MPI_Comm comm);
void generate_rows_distribution(char *matrix_file_name, float **matrix,int K, int M, int K_subarray,int M_subarray,int block_rows, int process_gird_cols, int my_rank, MPI_Comm comm);
void generate_rows_distribution_C(char *matrix_file_name, float **matrix,int N, int M, int N_subarray,int M_subarray,int block_rows, int process_gird_rows,int process_gird_cols, int my_rank, MPI_Comm comm);
void write_result_to_file(char *filename, float *result, int sub_matrix_size, int matrix_rows, int matrix_cols, int block_rows, int process_gird_rows, int process_gird_cols, int my_rank,MPI_Comm comm);


#endif