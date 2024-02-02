#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "../..//utils/file_io/stackio.h"
#include "../../utils/logger/logger.h"
#include "../../utils/string_formatter/formatter.h"

#ifndef LOG_MESSAGE_SIZE
#define LOG_MESSAGE_SIZE 256
#endif

void compute_sequential_matrix_by_matrix_multiplication(double *matrix_A,double *matrix_B,double *matrix_C,int N,int K, int M){
	for (int i=0; i<N; ++i){ 
        for(int k=0; k<K; ++k) {
            for(int j=0; j<M; ++j){    
                matrix_C[i * M + j] = matrix_C[i * M + j] + matrix_A[i * K + k] * matrix_B[k * M + j];
			}
    	}
	}
}

int main(int argc, char *argv[]){
	int N = 0;
    int K = 0;
	int M = 0;
    clock_t read_start, read_end, computation_start, computation_end;
	double cpu_time_used;
	double total_time = 0.0;
	char logger_message[LOG_MESSAGE_SIZE];
	char formatted_string[64];

    
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

	double *matrix_A;
	double *matrix_B;
	double *matrix_C;

    puts("");
	puts("---------------------------------------------------------------------------[START]---------------------------------------------------------------------------");


	memset(logger_message,0,LOG_MESSAGE_SIZE);
	sprintf(logger_message,"Starting sequential computation for matrix %d x %d",N,M);
	logger_info(logger_message);

	matrix_A = (double *)malloc(N*K*sizeof(double));
	matrix_B = (double *)malloc(K*M*sizeof(double));
	matrix_C = (double *)malloc(N*M*sizeof(double));
	

	logger_info("Matrix memory correctly allocated on HOST");


	// Getting data from previously generated file
	read_start = clock();
	read_matrix_from_file(matrix_A,N,K,(char *)"A");
	read_matrix_from_file(matrix_B,K,M,(char *)"B");
	read_matrix_from_file(matrix_C,N,M,(char *)"C");
	read_end = clock();
	cpu_time_used = ((double) (read_end - read_start)) / CLOCKS_PER_SEC;
	total_time += cpu_time_used;

	getFormattedTime(cpu_time_used,(char *)formatted_string);
	memset(logger_message,0,LOG_MESSAGE_SIZE);
	sprintf(logger_message,"Reading matrix A,B and C from memory ended in %s",formatted_string);
	logger_info(logger_message);

	logger_info("Matrix values acquired correctly on HOST");

#ifdef DEBUG
		puts("Matrix A:");
		print_matrix(matrix_A,N,K);
	
		puts("Matrix B:");
		print_matrix(matrix_B,K,M);
	
		puts("Matrix C:");
		print_matrix(matrix_C,N,M);
#endif

	
	computation_start = clock();
	compute_sequential_matrix_by_matrix_multiplication(matrix_A,matrix_B,matrix_C,N,K,M);
	computation_end = clock();
	cpu_time_used = ((double) (computation_end - computation_start)) / CLOCKS_PER_SEC;
	total_time += cpu_time_used;

	memset(formatted_string,0,64);
	getFormattedTime(cpu_time_used,(char *)formatted_string);

	memset(logger_message,0,LOG_MESSAGE_SIZE);
	sprintf(logger_message,"Sequential computation ended in %s",formatted_string);
	logger_info(logger_message);
	
#ifdef DEBUG
		puts("New matrix C:");
		print_matrix(matrix_C,N,M);
#endif

	write_matrix_to_file(matrix_C,N,M,(char *)"sequential_C");

	memset(logger_message,0,LOG_MESSAGE_SIZE);
	sprintf(logger_message,"Computation + I/O time: %f\n",total_time);
	logger_info(logger_message);
	puts("----------------------------------------------------------------------------[END]----------------------------------------------------------------------------");
	puts("");

	return 0;
}