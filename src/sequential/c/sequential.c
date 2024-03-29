#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>
#include <stdbool.h>
#include "../..//utils/file_io/stackio.h"
#include "../../utils/logger/logger.h"
#include "../../utils/string_formatter/formatter.h"

#ifndef LOG_MESSAGE_SIZE
#define LOG_MESSAGE_SIZE 256
#endif

#ifndef FORMATTED_STRING_SIZE
#define FORMATTED_STRING_SIZE 64
#endif

void compute_sequential_matrix_by_matrix_multiplication(float *matrix_A,float *matrix_B,float *matrix_C,int N,int K, int M){  
	for(int i=0; i<N; ++i){
		for(int k=0; k<K; ++k){
			for(int j=0; j<M; ++j){ 
                matrix_C[i * M + j] += matrix_A[i * K + k] * matrix_B[k * M + j];
			}
    	}
	}
}

int main(int argc, char *argv[]){
	int N = 0;
    int K = 0;
	int M = 0;

	struct timeval stop, start;	
	double total_time = 0.0;
	double read_from_file_time = 0.0;
	double computation_time = 0.0;
	double write_to_file_time = 0.0;


	char logger_message[LOG_MESSAGE_SIZE];
	char formatted_string[FORMATTED_STRING_SIZE];
	bool is_square;
	double gflops;
    
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

	is_square = (N==K && K==M) ? true : false;

	float *matrix_A;
	float *matrix_B;
	float *matrix_C;

    puts("");
	puts("---------------------------------------------------------------------------[START]---------------------------------------------------------------------------");


	memset(logger_message,0,LOG_MESSAGE_SIZE);
	sprintf(logger_message,"Starting sequential computation for matrix %d x %d",N,M);
	logger_info(logger_message);

	matrix_A = (float *)malloc(N*K*sizeof(float));
	matrix_B = (float *)malloc(K*M*sizeof(float));
	matrix_C = (float *)malloc(N*M*sizeof(float));
	

	logger_info("Matrix memory correctly allocated on HOST");


	// Getting data from previously generated file
	gettimeofday(&start, NULL);

	read_matrix_from_file(matrix_A,N,K,M,N,K,(char *)"A",is_square);
	read_matrix_from_file(matrix_B,N,K,M,K,M,(char *)"B",is_square);
	read_matrix_from_file(matrix_C,N,K,M,N,M,(char *)"C",is_square);
	
	gettimeofday(&stop, NULL);

	read_from_file_time = (double)(stop.tv_usec - start.tv_usec) / 1000000 + (double)(stop.tv_sec - start.tv_sec);
	total_time += read_from_file_time;

	getFormattedTime(read_from_file_time,(char *)formatted_string);
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

	
	gettimeofday(&start, NULL);

	compute_sequential_matrix_by_matrix_multiplication(matrix_A,matrix_B,matrix_C,N,K,M);
	
	gettimeofday(&stop, NULL);

	computation_time = (double)(stop.tv_usec - start.tv_usec) / 1000000 + (double)(stop.tv_sec - start.tv_sec);
	total_time += computation_time;

	memset(formatted_string,0,FORMATTED_STRING_SIZE);
	getFormattedTime(computation_time,(char *)formatted_string);

	memset(logger_message,0,LOG_MESSAGE_SIZE);
	sprintf(logger_message,"Sequential computation ended in %s",formatted_string);
	logger_info(logger_message);
	
#ifdef DEBUG
		puts("New matrix C:");
		print_matrix(matrix_C,N,M);
#endif

	gettimeofday(&start, NULL);

	write_matrix_to_file(matrix_C,N,K,M,N,M,(char *)"sequential_C",is_square);

	gettimeofday(&stop, NULL);

	write_to_file_time = (double)(stop.tv_usec - start.tv_usec) / 1000000 + (double)(stop.tv_sec - start.tv_sec);
	total_time += write_to_file_time;

	float partial = (float)(2.0*N*K*M)/total_time;
	gflops = (double)( partial / 1000000000);

	write_sequential_computation_csv(N,K,M,read_from_file_time,computation_time,write_to_file_time,gflops);


	memset(logger_message,0,LOG_MESSAGE_SIZE);
	sprintf(logger_message,"Computation + I/O time: %f\n",total_time);
	logger_info(logger_message);
	puts("----------------------------------------------------------------------------[END]----------------------------------------------------------------------------");
	puts("");

	return 0;
}