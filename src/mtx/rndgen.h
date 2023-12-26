#ifndef __MTX__
#define __MTX__


void print_matrix(float **matrix,size_t rows, size_t cols);
void generate_real_value_matrix(float **matrix,size_t rows, size_t cols);
void generate_zero_value_matrix(float **matrix,size_t rows, size_t cols);
void compute_sequential_matrix_by_matrix_multiplication(float **matrix_A,float **matrix_B,float **matrix_C,size_t rows,size_t cols);

#endif

