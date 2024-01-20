#ifndef STACK_IO
#define STACK_IO

void print_matrix(float **matrix,size_t rows, size_t cols);
size_t write_matrix_to_file(float **matrix, size_t rows, size_t cols, char *matrix_name);
void read_matrix_from_file(float **matrix,size_t rows_expected, size_t cols_expected, char *matrix_name);
void test();
void getFormattedTime(double time_in_seconds, char *formatted_String);

#endif