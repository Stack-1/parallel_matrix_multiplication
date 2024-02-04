#ifndef STACK_IO
#define STACK_IO

#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <stdbool.h>
#include "../logger/logger.h"

void print_matrix(float *matrix,int rows, int cols);
int write_matrix_to_file(float *matrix, int rows, int cols, char *matrix_name,bool is_square);
void read_matrix_from_file(float *matrix,int rows_expected, int cols_expected, char *matrix_name,bool is_square);
void read_matrix_from_file_mpi(float *matrix,int rows_expected, int cols_expected, char *matrix_name,bool is_square);
void write_times_to_txt_file(char *matrix_name, int rows, int cols, double time);
void write_sequential_computation_csv(int N,int K,int M,double read_time, double comptation_time,double write_time,double flops);

#endif