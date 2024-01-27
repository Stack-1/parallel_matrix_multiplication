#ifndef STACK_IO
#define STACK_IO


void print_matrix(float *matrix,int rows, int cols);
int write_matrix_to_file(float *matrix, int rows, int cols, char *matrix_name);
void read_matrix_from_file(float *matrix,int rows_expected, int cols_expected, char *matrix_name);

#endif