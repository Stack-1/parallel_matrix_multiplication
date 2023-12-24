#ifndef STACK_IO
#define STACK_IO

void print_matrix(float **matrix,size_t rows, size_t cols);
size_t write_matrix_to_file(float **matrix, size_t rows, size_t cols);
void read_matrix_from_file(char *file_name,float **matrix,size_t rows_expected, size_t cols_expected);
void test();

#endif